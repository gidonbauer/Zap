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

template <typename ActiveFloat, typename PassiveFloat, size_t DIM, Point2D_c PointType>
[[nodiscard]] constexpr auto point_on_boundary(const PointType& p,
                                               const Cell<ActiveFloat, PassiveFloat, DIM>& cell)
    -> bool {
  constexpr CoordType coord_type = PointType2CoordType<PointType>;

  return ((approx_eq(p.x, cell.template x_min<coord_type>()) ||
           approx_eq(p.x, cell.template x_min<coord_type>() + cell.template dx<coord_type>())) &&
          (p.y - cell.template y_min<coord_type>() >= -EPS<PassiveFloat> &&
           p.y - (cell.template y_min<coord_type>() + cell.template dy<coord_type>()) <=
               EPS<PassiveFloat>)) ||
         ((approx_eq(p.y, cell.template y_min<coord_type>()) ||
           approx_eq(p.y, cell.template y_min<coord_type>() + cell.template dy<coord_type>())) &&
          (p.x - cell.template x_min<coord_type>() >= -EPS<PassiveFloat> &&
           p.x - (cell.template x_min<coord_type>() + cell.template dx<coord_type>()) <=
               EPS<PassiveFloat>));
};

template <typename ActiveFloat, typename PassiveFloat, size_t DIM, Point2D_c PointType>
[[nodiscard]] constexpr auto point_in_cell(const PointType& p,
                                           const Cell<ActiveFloat, PassiveFloat, DIM>& cell) {
  constexpr CoordType coord_type = PointType2CoordType<PointType>;

  return p.x - cell.template x_min<coord_type>() >= -EPS<PassiveFloat> &&
         p.x - (cell.template x_min<coord_type>() + cell.template dx<coord_type>()) <=
             EPS<PassiveFloat> &&
         p.y - cell.template y_min<coord_type>() >= -EPS<PassiveFloat> &&
         p.y - (cell.template y_min<coord_type>() + cell.template dy<coord_type>()) <=
             EPS<PassiveFloat>;
};

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_GRID_HELPER_HPP_
