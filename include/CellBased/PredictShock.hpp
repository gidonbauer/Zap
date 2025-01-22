#ifndef ZAP_CELL_BASED_PREDICT_SHOCK_HPP_
#define ZAP_CELL_BASED_PREDICT_SHOCK_HPP_

#include "CellBased/Definitions.hpp"
#include "CellBased/Geometry.hpp"
#include "CellBased/Quadrature.hpp"

#include "Igor/Logging.hpp"

namespace Zap::CellBased {

// -------------------------------------------------------------------------------------------------
template <typename PassiveFloat>
[[nodiscard]] constexpr auto
eval_piecewise_linear(PassiveFloat t, const std::vector<SimCoord<PassiveFloat>>& points) noexcept
    -> std::pair<SimCoord<PassiveFloat>, SimCoord<PassiveFloat>> {
  IGOR_ASSERT(t >= 0 && t <= 1, "t must be in [0,1] but is {}", t);
  const auto segment = t * static_cast<PassiveFloat>(points.size() - 1);
  auto first_idx     = static_cast<size_t>(std::floor(segment));
  auto last_idx      = static_cast<size_t>(std::floor(segment + 1));
  if (last_idx == points.size()) {
    first_idx -= 1;
    last_idx -= 1;
  }

  IGOR_ASSERT(first_idx < points.size(), "Index {} is out of bounds.", first_idx);
  IGOR_ASSERT(last_idx < points.size(), "Index {} is out of bounds.", last_idx);
  const auto t_prime = segment - static_cast<PassiveFloat>(first_idx);
  IGOR_ASSERT(t_prime >= 0 && t_prime <= 1, "t_prime must be in [0,1] but is {}", t_prime);

  const auto value = points[first_idx] + t_prime * (points[last_idx] - points[first_idx]);

  if (first_idx != 0UZ && last_idx != points.size() - 1UZ &&
      (approx_eq(t_prime, PassiveFloat{0}) || approx_eq(t_prime, PassiveFloat{1}))) {
    Igor::Todo("Handle normal at kinks.");
  }
  const auto tangent = points[last_idx] - points[first_idx];
  const auto normal  = SimCoord<PassiveFloat>{tangent.y, -tangent.x};

  return std::make_pair(value, normal);
}

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, typename PassiveFloat>
[[nodiscard]] constexpr auto
predict_shock_front(const std::vector<SimCoord<ActiveFloat>>& original_shock,
                    PassiveFloat eps) noexcept -> std::vector<SimCoord<PassiveFloat>> {
  IGOR_ASSERT(original_shock.size() >= 2UZ,
              "Require at least two points on the shock curve but got only {}",
              original_shock.size());
  const auto N = original_shock.size();

  std::vector<SimCoord<PassiveFloat>> pred_shock(N);

  // First point
  {
    const SimCoord<PassiveFloat> p1 = {ad::value(original_shock[0].x),
                                       ad::value(original_shock[0].y)};
    const SimCoord<PassiveFloat> p2 = {ad::value(original_shock[1].x),
                                       ad::value(original_shock[1].y)};

    const SimCoord<PassiveFloat> t = p2 - p1;
    const SimCoord<PassiveFloat> n = SimCoord<PassiveFloat>{t.y, -t.x}.normalized();

    const SimCoord<PassiveFloat> d = {ad::derivative(original_shock[0].x),
                                      ad::derivative(original_shock[0].y)};

    const SimCoord<PassiveFloat> projected_d = d.dot(n) * n;

    pred_shock[0] = p1 + eps * projected_d;
  }

  // Intermediate points
  for (size_t i = 1; i < N - 1; ++i) {
    const SimCoord<PassiveFloat> p  = {ad::value(original_shock[i].x),
                                       ad::value(original_shock[i].y)};
    const SimCoord<PassiveFloat> pp = {ad::value(original_shock[i - 1].x),
                                       ad::value(original_shock[i - 1].y)};
    const SimCoord<PassiveFloat> pn = {ad::value(original_shock[i + 1].x),
                                       ad::value(original_shock[i + 1].y)};

    const SimCoord<PassiveFloat> tp = p - pp;
    const SimCoord<PassiveFloat> np = SimCoord<PassiveFloat>{tp.y, -tp.x}.normalized();

    const SimCoord<PassiveFloat> tn = pn - p;
    const SimCoord<PassiveFloat> nn = SimCoord<PassiveFloat>{tn.y, -tn.x}.normalized();

    const SimCoord<PassiveFloat> n = ((np + nn) / 2).normalized();

    const SimCoord<PassiveFloat> d           = {ad::derivative(original_shock[i].x),
                                                ad::derivative(original_shock[i].y)};
    const SimCoord<PassiveFloat> projected_d = d.dot(n) * n;

    pred_shock[i] = p + eps * projected_d;
  }

  // Last point
  {
    const SimCoord<PassiveFloat> p1 = {ad::value(original_shock[N - 2].x),
                                       ad::value(original_shock[N - 2].y)};
    const SimCoord<PassiveFloat> p2 = {ad::value(original_shock[N - 1].x),
                                       ad::value(original_shock[N - 1].y)};

    const SimCoord<PassiveFloat> t = p2 - p1;
    const SimCoord<PassiveFloat> n = SimCoord<PassiveFloat>{t.y, -t.x}.normalized();

    const SimCoord<PassiveFloat> d = {ad::derivative(original_shock[N - 1].x),
                                      ad::derivative(original_shock[N - 1].y)};

    const SimCoord<PassiveFloat> projected_d = d.dot(n) * n;

    pred_shock[N - 1] = p1 + eps * projected_d;
  }

  return pred_shock;
}

// -------------------------------------------------------------------------------------------------
template <size_t N = 15UZ, bool simple_dif_func = false, typename PassiveFloat>
[[nodiscard]] constexpr auto compare_curves(const std::vector<SimCoord<PassiveFloat>>& real_curve,
                                            const std::vector<SimCoord<PassiveFloat>>& pred_curve,
                                            size_t num_slices = 1) noexcept -> PassiveFloat {
  IGOR_ASSERT(num_slices > 0,
              "Number of slices for the calculation of the integral must be larger than zero.");

  auto dif_func = [&](PassiveFloat t) {
    if constexpr (simple_dif_func) {
      return (std::get<0>(eval_piecewise_linear(t, pred_curve)) -
              std::get<0>(eval_piecewise_linear(t, real_curve)))
          .norm();
    } else {
      const auto [p, n] = eval_piecewise_linear(t, pred_curve);

      auto min_dist = std::numeric_limits<PassiveFloat>::max();
      for (size_t i = 0; i < real_curve.size() - 1; ++i) {
        const auto intersect =
            Geometry::line_intersect<true>(p, p + n, real_curve[i], real_curve[i + 1]);
        if (intersect.has_value()) {
          const auto dist = (p - *intersect).norm();
          min_dist        = std::min(min_dist, dist);
        }
      }

      if (min_dist == std::numeric_limits<PassiveFloat>::max()) {
        const auto other = std::get<0>(eval_piecewise_linear(t, real_curve));
        min_dist         = (p - other).norm();
      }

      return min_dist;
    }
  };

  PassiveFloat integral = 0;
  for (size_t i = 0; i < num_slices; ++i) {
    const PassiveFloat begin = static_cast<PassiveFloat>(i) / static_cast<PassiveFloat>(num_slices);
    const PassiveFloat end =
        static_cast<PassiveFloat>(i + 1) / static_cast<PassiveFloat>(num_slices);
    integral += quadrature<N>(dif_func, begin, end);
  }
  return integral;
}

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_PREDICT_SHOCK_HPP_
