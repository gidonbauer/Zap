#ifndef ZAP_CELL_BASED_GRID_HELPER_HPP_
#define ZAP_CELL_BASED_GRID_HELPER_HPP_

#include <vector>

#include <Eigen/Dense>

#include "CellBased/Cell.hpp"

namespace Zap::CellBased {

enum : size_t { X, Y, POINT_SIZE };

template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto point_on_boundary(const Eigen::Vector<Float, POINT_SIZE>& p,
                                               const Cell<Float, DIM>& cell) -> bool {

  const auto approx_eq = [](Float a, Float b) -> bool {
    constexpr Float eps = 1e-8;
    return std::abs(a - b) <= eps;
  };

  return ((approx_eq(p(X), cell.x_min) || approx_eq(p(X), cell.x_min + cell.dx)) &&
          (p(Y) >= cell.y_min && p(Y) <= cell.y_min + cell.dy)) ||
         ((approx_eq(p(Y), cell.y_min) || approx_eq(p(Y), cell.y_min + cell.dy)) &&
          (p(X) >= cell.x_min && p(X) <= cell.x_min + cell.dx));
};

template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto point_in_cell(const Eigen::Vector<Float, POINT_SIZE>& p,
                                           const Cell<Float, DIM>& cell) {
  constexpr Float eps = 1e-8;

  return p(X) - cell.x_min >= -eps &&             //
         p(X) - (cell.x_min + cell.dx) <= eps &&  //
         p(Y) - cell.y_min >= -eps &&             //
         p(Y) - (cell.y_min + cell.dy) <= eps;
};

template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto
piecewise_linear_curve(Float t,
                       const std::vector<Eigen::Vector<Float, POINT_SIZE>>& points) noexcept
    -> Eigen::Vector<Float, POINT_SIZE> {
  assert(t >= 0 && t <= 1);

  const auto segment = t * static_cast<Float>(points.size() - 1);
  auto first_idx     = static_cast<size_t>(std::floor(segment));
  auto last_idx      = static_cast<size_t>(std::floor(segment + 1));
  if (last_idx == points.size()) {
    first_idx -= 1;
    last_idx -= 1;
  }

  // IGOR_DEBUG_PRINT(t);
  // IGOR_DEBUG_PRINT(segment);
  // IGOR_DEBUG_PRINT(first_idx);
  // IGOR_DEBUG_PRINT(last_idx);
  // IGOR_DEBUG_PRINT(points.size());
  assert(first_idx < points.size());
  assert(last_idx < points.size());

  const auto t_prime = segment - static_cast<Float>(first_idx);

  // IGOR_DEBUG_PRINT(t_prime);
  // std::cerr << "----------------------------------------\n";

  return points[first_idx] + t_prime * (points[last_idx] - points[first_idx]);
}

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_GRID_HELPER_HPP_
