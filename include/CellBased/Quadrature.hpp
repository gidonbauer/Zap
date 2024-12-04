#ifndef ZAP_CELL_BASED_QUADRATURE_HPP_
#define ZAP_CELL_BASED_QUADRATURE_HPP_

#include <cmath>
#ifdef ZAP_QUAD_2
#include <numbers>
#endif  // ZAP_QUAD_2

#include "Geometry.hpp"

#include "Igor/Math.hpp"

namespace Zap::CellBased {

namespace detail {

// -------------------------------------------------------------------------------------------------
template <typename FUNC, Point2D_c PointType>
[[nodiscard]] constexpr auto
quadrature_four_corners(FUNC f, const Geometry::Polygon<PointType>& domain) noexcept
    -> decltype(std::declval<PointType>().x) {
  static_assert(
      std::is_same_v<decltype(std::declval<PointType>().x), decltype(std::declval<PointType>().y)>,
      "Expect x- and y-component of PointType to have the same type.");

  using Float = decltype(std::declval<PointType>().x);

  IGOR_ASSERT(domain.size() == 4UZ, "Domain has to have exactly four corners.");

#ifdef ZAP_QUAD_64
  constexpr std::array<Float, 64UZ> xs{
      -0.0243502926634244, 0.0243502926634244, -0.0729931217877990, 0.0729931217877990,
      -0.1214628192961206, 0.1214628192961206, -0.1696444204239928, 0.1696444204239928,
      -0.2174236437400071, 0.2174236437400071, -0.2646871622087674, 0.2646871622087674,
      -0.3113228719902110, 0.3113228719902110, -0.3572201583376681, 0.3572201583376681,
      -0.4022701579639916, 0.4022701579639916, -0.4463660172534641, 0.4463660172534641,
      -0.4894031457070530, 0.4894031457070530, -0.5312794640198946, 0.5312794640198946,
      -0.5718956462026340, 0.5718956462026340, -0.6111553551723933, 0.6111553551723933,
      -0.6489654712546573, 0.6489654712546573, -0.6852363130542333, 0.6852363130542333,
      -0.7198818501716109, 0.7198818501716109, -0.7528199072605319, 0.7528199072605319,
      -0.7839723589433414, 0.7839723589433414, -0.8132653151227975, 0.8132653151227975,
      -0.8406292962525803, 0.8406292962525803, -0.8659993981540928, 0.8659993981540928,
      -0.8893154459951141, 0.8893154459951141, -0.9105221370785028, 0.9105221370785028,
      -0.9295691721319396, 0.9295691721319396, -0.9464113748584028, 0.9464113748584028,
      -0.9610087996520538, 0.9610087996520538, -0.9733268277899110, 0.9733268277899110,
      -0.9833362538846260, 0.9833362538846260, -0.9910133714767443, 0.9910133714767443,
      -0.9963401167719553, 0.9963401167719553, -0.9993050417357722, 0.9993050417357722,
  };
  constexpr std::array<Float, 64UZ> ws{
      0.0486909570091397, 0.0486909570091397, 0.0485754674415034, 0.0485754674415034,
      0.0483447622348030, 0.0483447622348030, 0.0479993885964583, 0.0479993885964583,
      0.0475401657148303, 0.0475401657148303, 0.0469681828162100, 0.0469681828162100,
      0.0462847965813144, 0.0462847965813144, 0.0454916279274181, 0.0454916279274181,
      0.0445905581637566, 0.0445905581637566, 0.0435837245293235, 0.0435837245293235,
      0.0424735151236536, 0.0424735151236536, 0.0412625632426235, 0.0412625632426235,
      0.0399537411327203, 0.0399537411327203, 0.0385501531786156, 0.0385501531786156,
      0.0370551285402400, 0.0370551285402400, 0.0354722132568824, 0.0354722132568824,
      0.0338051618371416, 0.0338051618371416, 0.0320579283548516, 0.0320579283548516,
      0.0302346570724025, 0.0302346570724025, 0.0283396726142595, 0.0283396726142595,
      0.0263774697150547, 0.0263774697150547, 0.0243527025687109, 0.0243527025687109,
      0.0222701738083833, 0.0222701738083833, 0.0201348231535302, 0.0201348231535302,
      0.0179517157756973, 0.0179517157756973, 0.0157260304760247, 0.0157260304760247,
      0.0134630478967186, 0.0134630478967186, 0.0111681394601311, 0.0111681394601311,
      0.0088467598263639, 0.0088467598263639, 0.0065044579689784, 0.0065044579689784,
      0.0041470332605625, 0.0041470332605625, 0.0017832807216964, 0.0017832807216964,
  };
#elifdef ZAP_QUAD_5
  constexpr std::array<Float, 5UZ> xs{
      0.0,
      1.0 / 3.0 * Igor::constexpr_sqrt(5.0 - 2.0 * Igor::constexpr_sqrt(10.0 / 7.0)),
      -1.0 / 3.0 * Igor::constexpr_sqrt(5.0 - 2.0 * Igor::constexpr_sqrt(10.0 / 7.0)),
      1.0 / 3.0 * Igor::constexpr_sqrt(5.0 + 2.0 * Igor::constexpr_sqrt(10.0 / 7.0)),
      -1.0 / 3.0 * Igor::constexpr_sqrt(5.0 + 2.0 * Igor::constexpr_sqrt(10.0 / 7.0)),
  };
  constexpr std::array<Float, 5UZ> ws{
      128.0 / 255.0,
      (322.0 + 13.0 * Igor::constexpr_sqrt(70.0)) / 900.0,
      (322.0 + 13.0 * Igor::constexpr_sqrt(70.0)) / 900.0,
      (322.0 - 13.0 * Igor::constexpr_sqrt(70.0)) / 900.0,
      (322.0 - 13.0 * Igor::constexpr_sqrt(70.0)) / 900.0,
  };
#elifdef ZAP_QUAD_2
  constexpr std::array<Float, 2UZ> xs{
      1.0 / std::numbers::sqrt3_v<Float>,
      -1.0 / std::numbers::sqrt3_v<Float>,
  };
  constexpr std::array<Float, 2UZ> ws{
      1.0,
      1.0,
  };
#else
  constexpr std::array<Float, 15UZ> xs{
      0.0000000000000000,
      -0.2011940939974345,
      0.2011940939974345,
      -0.3941513470775634,
      0.3941513470775634,
      -0.5709721726085388,
      0.5709721726085388,
      -0.7244177313601701,
      0.7244177313601701,
      -0.8482065834104272,
      0.8482065834104272,
      -0.9372733924007060,
      0.9372733924007060,
      -0.9879925180204854,
      0.9879925180204854,
  };
  constexpr std::array<Float, 15UZ> ws{
      0.2025782419255613,
      0.1984314853271116,
      0.1984314853271116,
      0.1861610000155622,
      0.1861610000155622,
      0.1662692058169939,
      0.1662692058169939,
      0.1395706779261543,
      0.1395706779261543,
      0.1071592204671719,
      0.1071592204671719,
      0.0703660474881081,
      0.0703660474881081,
      0.0307532419961173,
      0.0307532419961173,
  };
#endif
  static_assert(xs.size() == ws.size());

  // Jacobian matrix of coordinate transform
  const Float J11     = (domain[1].x + domain[2].x) - (domain[0].x + domain[3].x);
  const Float J12     = (domain[2].x + domain[3].x) - (domain[0].x + domain[1].x);
  const Float J21     = (domain[1].y + domain[2].y) - (domain[0].y + domain[3].y);
  const Float J22     = (domain[2].y + domain[3].y) - (domain[0].y + domain[1].y);
  const Float abs_det = (1.0 / 16.0) * std::abs(J11 * J22 - J12 * J21);

  constexpr std::array<Float (*)(Float, Float), 4UZ> psi = {
      [](Float xi, Float eta) -> Float { return (1.0 / 4.0) * (1.0 - xi) * (1.0 - eta); },
      [](Float xi, Float eta) -> Float { return (1.0 / 4.0) * (1.0 + xi) * (1.0 - eta); },
      [](Float xi, Float eta) -> Float { return (1.0 / 4.0) * (1.0 + xi) * (1.0 + eta); },
      [](Float xi, Float eta) -> Float { return (1.0 / 4.0) * (1.0 - xi) * (1.0 + eta); },
  };

  auto integral = static_cast<Float>(0);
  for (size_t xidx = 0; xidx < xs.size(); ++xidx) {
    for (size_t yidx = 0; yidx < xs.size(); ++yidx) {
      const auto xi  = xs[xidx];
      const auto wx  = ws[xidx];
      const auto eta = xs[yidx];
      const auto wy  = ws[yidx];

      auto x = static_cast<Float>(0);
      auto y = static_cast<Float>(0);
      for (size_t i = 0; i < psi.size(); ++i) {
        x += psi[i](xi, eta) * domain[i].x;
        y += psi[i](xi, eta) * domain[i].y;
      }

      // Igor::Debug("xi = {}, eta = {} => x = {}, y = {}", xi, eta, x, y);

      integral += wx * wy * f(x, y);
    }
  }
  return abs_det * integral;
}

}  // namespace detail

// -------------------------------------------------------------------------------------------------
template <typename FUNC, Point2D_c PointType>
[[nodiscard]] constexpr auto quadrature(FUNC f, Geometry::Polygon<PointType> domain) noexcept
    -> decltype(std::declval<PointType>().x) {
  static_assert(
      std::is_same_v<decltype(std::declval<PointType>().x), decltype(std::declval<PointType>().y)>,
      "Expect x- and y-component of PointType to have the same type.");

  if (domain.size() == 3UZ) {
    domain.add_point((domain[0] + domain[1]) / 2);
    return detail::quadrature_four_corners(f, domain);
  } else if (domain.size() == 4UZ) {
    return detail::quadrature_four_corners(f, domain);
  } else if (domain.size() == 5UZ) {
    const Geometry::Polygon<PointType> domain1({
        domain[0],
        domain[1],
        domain[2],
        domain[3],
    });

    const Geometry::Polygon<PointType> domain2(
        {domain[3], domain[4], domain[0], (domain[3] + domain[0]) / 2});

    return detail::quadrature_four_corners(f, domain1) +
           detail::quadrature_four_corners(f, domain2);
  } else {
    Igor::Todo("{} corners is not implemented yet.", domain.size());
    std::unreachable();
  }
}

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_QUADRATURE_HPP_
