#ifndef ZAP_CELL_BASED_QUADRATURE_HPP_
#define ZAP_CELL_BASED_QUADRATURE_HPP_

#include <cmath>
#include <numeric>

#include <AD/ad.hpp>

#include "CellBased/Geometry.hpp"
#include "CellBased/QuadratureTables.hpp"

namespace Zap::CellBased {

namespace detail {

// -------------------------------------------------------------------------------------------------
template <size_t N = 15UZ, typename FUNC, Point2D_c PointType>
[[nodiscard]] constexpr auto
quadrature_four_corners(FUNC f, const Geometry::Polygon<PointType>& domain) noexcept
    -> decltype(std::declval<PointType>().x) {
  static_assert(N > 0UZ && N <= detail::MAX_QUAD_N);
  static_assert(
      std::is_same_v<decltype(std::declval<PointType>().x), decltype(std::declval<PointType>().y)>,
      "Expect x- and y-component of PointType to have the same type.");

  using ActiveFloat  = decltype(std::declval<PointType>().x);
  using PassiveFloat = std::conditional_t<ad::mode<ActiveFloat>::is_ad_type,
                                          typename ad::mode<ActiveFloat>::passive_t,
                                          ActiveFloat>;

  IGOR_ASSERT(domain.size() == 4UZ, "Domain has to have exactly four corners.");

  constexpr auto& gauss_points  = detail::gauss_points_table<PassiveFloat>[N - 1UZ];
  constexpr auto& gauss_weights = detail::gauss_weights_table<PassiveFloat>[N - 1UZ];
  static_assert(gauss_points.size() == gauss_weights.size(),
                "Weights and points must have the same size.");
  static_assert(gauss_points.size() == N, "Number of weights and points must be equal to N.");
  static_assert(approx_eq(std::reduce(gauss_weights.cbegin(), gauss_weights.cend()),
                          static_cast<PassiveFloat>(2)),
                "Weights must add up to 2.");

  // Jacobian matrix of coordinate transform
  // J11
  auto dxdxi = [&domain](PassiveFloat eta) -> ActiveFloat {
    return 1.0 / 4.0 *
           ((domain[1].x - domain[0].x) * (1 - eta) + (domain[2].x - domain[3].x) * (1 + eta));
  };
  // J12
  auto dydxi = [&domain](PassiveFloat eta) -> ActiveFloat {
    return 1.0 / 4.0 *
           ((domain[1].y - domain[0].y) * (1 - eta) + (domain[2].y - domain[3].y) * (1 + eta));
  };
  // J21
  auto dxdeta = [&domain](PassiveFloat xi) -> ActiveFloat {
    return 1.0 / 4.0 *
           ((domain[3].x - domain[0].x) * (1 - xi) + (domain[2].x - domain[1].x) * (1 + xi));
  };
  // J22
  auto dydeta = [&domain](PassiveFloat xi) -> ActiveFloat {
    return 1.0 / 4.0 *
           ((domain[3].y - domain[0].y) * (1 - xi) + (domain[2].y - domain[1].y) * (1 + xi));
  };

  auto abs_det_J = [&](PassiveFloat xi, PassiveFloat eta) -> ActiveFloat {
    return std::abs(dxdxi(eta) * dydeta(xi) - dydxi(eta) * dxdeta(xi));
  };

  // Shape functions
  constexpr std::array<PassiveFloat (*)(PassiveFloat, PassiveFloat), 4UZ> psi = {
      [](PassiveFloat xi, PassiveFloat eta) -> PassiveFloat {
        return (1.0 / 4.0) * (1.0 - xi) * (1.0 - eta);
      },
      [](PassiveFloat xi, PassiveFloat eta) -> PassiveFloat {
        return (1.0 / 4.0) * (1.0 + xi) * (1.0 - eta);
      },
      [](PassiveFloat xi, PassiveFloat eta) -> PassiveFloat {
        return (1.0 / 4.0) * (1.0 + xi) * (1.0 + eta);
      },
      [](PassiveFloat xi, PassiveFloat eta) -> PassiveFloat {
        return (1.0 / 4.0) * (1.0 - xi) * (1.0 + eta);
      },
  };

  auto integral = static_cast<ActiveFloat>(0);
  for (size_t xidx = 0; xidx < gauss_points.size(); ++xidx) {
    for (size_t yidx = 0; yidx < gauss_points.size(); ++yidx) {
      const auto xi  = gauss_points[xidx];
      const auto wx  = gauss_weights[xidx];
      const auto eta = gauss_points[yidx];
      const auto wy  = gauss_weights[yidx];

      auto x = static_cast<ActiveFloat>(0);
      auto y = static_cast<ActiveFloat>(0);
      for (size_t i = 0; i < psi.size(); ++i) {
        x += psi[i](xi, eta) * domain[i].x;
        y += psi[i](xi, eta) * domain[i].y;
      }

      integral += wx * wy * f(x, y) * abs_det_J(xi, eta);
    }
  }
  return integral;
}

}  // namespace detail

// -------------------------------------------------------------------------------------------------
template <size_t N = 15UZ, typename FUNC, Point2D_c PointType>
[[nodiscard]] constexpr auto quadrature(FUNC f, Geometry::Polygon<PointType> domain) noexcept
    -> decltype(std::declval<PointType>().x) {
  static_assert(N > 0UZ && N <= detail::MAX_QUAD_N);
  static_assert(
      std::is_same_v<decltype(std::declval<PointType>().x), decltype(std::declval<PointType>().y)>,
      "Expect x- and y-component of PointType to have the same type.");

  if (domain.size() == 3UZ) {
    // TODO: Do something smarter for triangles.
    domain.add_point((domain[0] + domain[1]) / 2);
    return detail::quadrature_four_corners<N>(f, domain);
  } else if (domain.size() == 4UZ) {
    return detail::quadrature_four_corners<N>(f, domain);
  } else if (domain.size() == 5UZ) {
    const Geometry::Polygon<PointType> domain1({
        domain[0],
        domain[1],
        domain[2],
        domain[3],
    });

    // TODO: Do something smarter for the resulting triangle.
    const Geometry::Polygon<PointType> domain2(
        {domain[3], domain[4], domain[0], (domain[3] + domain[0]) / 2});

    return detail::quadrature_four_corners<N>(f, domain1) +
           detail::quadrature_four_corners<N>(f, domain2);
  } else {
    Igor::Todo("{} corners is not implemented yet.", domain.size());
    std::unreachable();
  }
}

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_QUADRATURE_HPP_
