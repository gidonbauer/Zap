#ifndef ZAP_CELL_BASED_RECONSTRUCT_SHOCK_HPP_
#define ZAP_CELL_BASED_RECONSTRUCT_SHOCK_HPP_

#include "CellBased/Definitions.hpp"

#include "Igor/Logging.hpp"

namespace Zap::CellBased {

// TODO: Think about spline reconstruction.

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, typename PassiveFloat>
[[nodiscard]] constexpr auto
piecewise_linear(PassiveFloat t, const std::vector<SimCoord<ActiveFloat>>& points) noexcept
    -> SimCoord<ActiveFloat> {
  assert(t >= 0 && t <= 1);
  const auto segment = t * static_cast<PassiveFloat>(points.size() - 1);
  auto first_idx     = static_cast<size_t>(std::floor(segment));
  auto last_idx      = static_cast<size_t>(std::floor(segment + 1));
  if (last_idx == points.size()) {
    first_idx -= 1;
    last_idx -= 1;
  }

  assert(first_idx < points.size());
  assert(last_idx < points.size());
  const auto t_prime = segment - static_cast<PassiveFloat>(first_idx);
  assert(t_prime >= 0 && t_prime <= 1);

  return points[first_idx] + t_prime * (points[last_idx] - points[first_idx]);
}

template <typename ActiveFloat, typename PassiveFloat>
[[nodiscard]] constexpr auto
piecewise_linear(const std::vector<PassiveFloat>& ts,
                 const std::vector<SimCoord<ActiveFloat>>& points) noexcept
    -> std::vector<SimCoord<ActiveFloat>> {
  std::vector<SimCoord<ActiveFloat>> curve_points(ts.size());
  std::transform(std::cbegin(ts),
                 std::cend(ts),
                 std::begin(curve_points),
                 [&points](PassiveFloat t) { return piecewise_linear(t, points); });
  return curve_points;
}

// -------------------------------------------------------------------------------------------------
namespace detail {

template <typename T>
struct Spline {
  T a, b, c, d, x;
};

// Algorithm taken from: https://en.wikipedia.org/wiki/Spline_(mathematics)
template <typename ActiveFloat, typename PassiveFloat>
[[nodiscard]] constexpr auto natural_cubic_spline(const std::vector<PassiveFloat>& x,
                                                  const std::vector<ActiveFloat>& y) noexcept
    -> std::vector<Spline<ActiveFloat>> {
  IGOR_ASSERT(x.size() == y.size(),
              "Expected x and y to have the same size, but x.size() = {} and y.size() = {}.",
              x.size(),
              y.size());
  IGOR_ASSERT(x.size() > 1, "Require at least two points, but got only {}.", x.size());
  const auto n = x.size() - 1;

  std::vector<ActiveFloat> a(n + 1);
  for (size_t i = 0; i < n + 1; ++i) {
    a[i] = y[i];
  }

  std::vector<PassiveFloat> h(n);
  for (size_t i = 0; i < n; ++i) {
    h[i] = x[i + 1] - x[i];
  }

  std::vector<ActiveFloat> alpha(n);
  alpha[0] = PassiveFloat{0.0};
  for (size_t i = 1; i < n; ++i) {
    alpha[i] = PassiveFloat{3.0} * (a[i + 1] - a[i]) / h[i] -
               PassiveFloat{3.0} * (a[i] - a[i - 1]) / h[i - 1];
  }

  std::vector<PassiveFloat> l(n + 1);
  std::vector<PassiveFloat> mu(n + 1);
  std::vector<ActiveFloat> z(n + 1);

  l[0]  = PassiveFloat{1.0};
  mu[0] = PassiveFloat{0.0};
  z[0]  = PassiveFloat{0.0};

  for (size_t i = 1; i < n; ++i) {
    l[i]  = PassiveFloat{2.0} * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
    mu[i] = h[i] / l[i];
    z[i]  = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
  }

  l[n] = PassiveFloat{1.0};
  z[n] = PassiveFloat{0.0};

  std::vector<ActiveFloat> b(n);
  std::vector<ActiveFloat> c(n + 1);
  c[n] = PassiveFloat{0.0};
  std::vector<ActiveFloat> d(n);

  for (int j_ = static_cast<int>(n - 1); j_ >= 0; --j_) {
    const auto j = static_cast<size_t>(j_);

    c[j] = z[j] - mu[j] * c[j + 1];
    b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / PassiveFloat{3.0};
    d[j] = (c[j + 1] - c[j]) / (PassiveFloat{3.0} * h[j]);
  }

  std::vector<Spline<ActiveFloat>> splines(n);
  for (size_t i = 0; i < n; ++i) {
    auto& s = splines[i];
    s.a     = a[i];
    s.b     = b[i];
    s.c     = c[i];
    s.d     = d[i];
    s.x     = x[i];
  }

  return splines;
}

template <typename ActiveFloat, typename PassiveFloat>
[[nodiscard]] constexpr auto
eval_natural_cubic_spline(PassiveFloat t,
                          const std::vector<Spline<ActiveFloat>>& x_spline,
                          const std::vector<Spline<ActiveFloat>>& y_spline) noexcept
    -> SimCoord<ActiveFloat> {
  // Assumes that x-component of the splines are eqaully spaces

  IGOR_ASSERT(x_spline.size() == y_spline.size(),
              "Expected x_spline and y_spline to have the same size, but x_spline.size() = {} and "
              "y_spline.size() = {}.",
              x_spline.size(),
              y_spline.size());

  auto idx = static_cast<size_t>(t * static_cast<PassiveFloat>(x_spline.size() - 1));
  if (idx == x_spline.size()) { idx -= 1; }
  IGOR_ASSERT(idx < x_spline.size(), "idx {} is out of bounds.", idx);

  const auto& sx = x_spline[idx];
  const auto& sy = y_spline[idx];

  auto eval = [t](const Spline<ActiveFloat>& s) {
    return s.a + s.b * (t - s.x) + s.c * (t - s.x) * (t - s.x) +
           s.d * (t - s.x) * (t - s.x) * (t - s.x);
  };

  return SimCoord<ActiveFloat>{eval(sx), eval(sy)};
}

}  // namespace detail

template <typename ActiveFloat, typename PassiveFloat>
[[nodiscard]] constexpr auto
natural_cubic_spline(const std::vector<PassiveFloat>& ts,
                     const std::vector<SimCoord<ActiveFloat>>& points) noexcept
    -> std::vector<SimCoord<ActiveFloat>> {

  std::vector<PassiveFloat> rs(points.size());
  std::generate(std::begin(rs), std::end(rs), [i = 0UZ, n = rs.size()]() mutable {
    return static_cast<PassiveFloat>(i++) / static_cast<PassiveFloat>(n - 1);
  });
  std::vector<ActiveFloat> xs(points.size());
  std::transform(std::cbegin(points),
                 std::cend(points),
                 std::begin(xs),
                 [](const SimCoord<ActiveFloat>& p) { return p.x; });
  std::vector<ActiveFloat> ys(points.size());
  std::transform(std::cbegin(points),
                 std::cend(points),
                 std::begin(ys),
                 [](const SimCoord<ActiveFloat>& p) { return p.y; });

  const auto x_spline = detail::natural_cubic_spline(rs, xs);
  const auto y_spline = detail::natural_cubic_spline(rs, ys);

  std::vector<SimCoord<ActiveFloat>> curve_points(ts.size());
  std::transform(std::cbegin(ts),
                 std::cend(ts),
                 std::begin(curve_points),
                 [&x_spline, &y_spline](PassiveFloat t) {
                   return detail::eval_natural_cubic_spline(t, x_spline, y_spline);
                 });
  return curve_points;
}

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, typename PassiveFloat>
[[nodiscard]] constexpr auto smoothstep(PassiveFloat t,
                                        const std::vector<SimCoord<ActiveFloat>>& points) noexcept
    -> SimCoord<ActiveFloat> {
  assert(t >= 0 && t <= 1);
  const auto segment = t * static_cast<PassiveFloat>(points.size() - 1);
  auto first_idx     = static_cast<size_t>(std::floor(segment));
  auto last_idx      = static_cast<size_t>(std::floor(segment + 1));
  if (last_idx == points.size()) {
    first_idx -= 1;
    last_idx -= 1;
  }

  assert(first_idx < points.size());
  assert(last_idx < points.size());
  auto t_prime = segment - static_cast<PassiveFloat>(first_idx);
  t_prime      = 3 * t_prime * t_prime - 2 * t_prime * t_prime * t_prime;

  return points[first_idx] + t_prime * (points[last_idx] - points[first_idx]);
}

template <typename ActiveFloat, typename PassiveFloat>
[[nodiscard]] constexpr auto smoothstep(const std::vector<PassiveFloat>& ts,
                                        const std::vector<SimCoord<ActiveFloat>>& points) noexcept
    -> std::vector<SimCoord<ActiveFloat>> {
  std::vector<SimCoord<ActiveFloat>> curve_points(ts.size());
  std::transform(std::cbegin(ts),
                 std::cend(ts),
                 std::begin(curve_points),
                 [&points](PassiveFloat t) { return smoothstep(t, points); });
  return curve_points;
}

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_RECONSTRUCT_SHOCK_HPP_
