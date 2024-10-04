#ifndef ZAP_CELL_BASED_GRID_HELPER_HPP_
#define ZAP_CELL_BASED_GRID_HELPER_HPP_

#include <Eigen/Dense>

#include "CellBased/Cell.hpp"
#include "CellBased/Definitions.hpp"

namespace Zap::CellBased {

template <typename Float>
[[nodiscard]] constexpr auto approx_eq(Float a, Float b, Float tol = EPS<Float>) noexcept -> bool {
  return std::abs(a - b) <= tol;
};

template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto point_on_boundary_grid_coord(const Point<Float>& p,
                                                          const Cell<Float, DIM>& cell) -> bool {
  // TODO: Do we need to add tolerance to the second check that makes sure that we are within the
  // cell? (line 23 and 26)
  return ((approx_eq(p(X), static_cast<Float>(cell.x_idx)) ||
           approx_eq(p(X), static_cast<Float>(cell.x_idx + 1))) &&
          (p(Y) - static_cast<Float>(cell.y_idx) >= -EPS<Float> &&
           p(Y) - static_cast<Float>(cell.y_idx + 1) <= EPS<Float>)) ||
         ((approx_eq(p(Y), static_cast<Float>(cell.y_idx)) ||
           approx_eq(p(Y), static_cast<Float>(cell.y_idx + 1))) &&
          (p(X) - static_cast<Float>(cell.x_idx) >= -EPS<Float> &&
           p(X) - static_cast<Float>(cell.x_idx + 1) <= EPS<Float>));
};

template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto point_in_cell_grid_coord(const Point<Float>& p,
                                                      const Cell<Float, DIM>& cell) {
  return p(X) - static_cast<Float>(cell.x_idx) >= -EPS<Float> &&     //
         p(X) - static_cast<Float>(cell.x_idx + 1) <= EPS<Float> &&  //
         p(Y) - static_cast<Float>(cell.y_idx) >= -EPS<Float> &&     //
         p(Y) - static_cast<Float>(cell.y_idx + 1) <= EPS<Float>;
};

// template <typename Float, size_t DIM>
// [[nodiscard]] constexpr auto
// piecewise_linear_curve(Float t, const std::vector<Point<Float>>& points) noexcept -> Point<Float>
// {
//   assert(t >= 0 && t <= 1);
//
//   const auto segment = t * static_cast<Float>(points.size() - 1);
//   auto first_idx     = static_cast<size_t>(std::floor(segment));
//   auto last_idx      = static_cast<size_t>(std::floor(segment + 1));
//   if (last_idx == points.size()) {
//     first_idx -= 1;
//     last_idx -= 1;
//   }
//
//   // IGOR_DEBUG_PRINT(t);
//   // IGOR_DEBUG_PRINT(segment);
//   // IGOR_DEBUG_PRINT(first_idx);
//   // IGOR_DEBUG_PRINT(last_idx);
//   // IGOR_DEBUG_PRINT(points.size());
//   assert(first_idx < points.size());
//   assert(last_idx < points.size());
//
//   const auto t_prime = segment - static_cast<Float>(first_idx);
//
//   // IGOR_DEBUG_PRINT(t_prime);
//   // std::cerr << "----------------------------------------\n";
//
//   return points[first_idx] + t_prime * (points[last_idx] - points[first_idx]);
// }

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_GRID_HELPER_HPP_
