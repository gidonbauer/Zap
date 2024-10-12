#ifndef ZAP_CELL_BASED_GRID_HPP_
#define ZAP_CELL_BASED_GRID_HPP_

#include <numeric>

#include "CellBased/Cell.hpp"
#include "CellBased/Definitions.hpp"
#include "CellBased/GridHelper.hpp"
#include "IO/Fwd.hpp"

#include "Igor/TypeName.hpp"

namespace Zap::CellBased {

// IDEA: Introduce scaled grid
//        - each cell is 1x1; global scaling factor for x and y
//        - each cell only holds its indices (and value)
//        - cut are always in [0, 1]
//        - "grid coordinates" vs "simulation coordinates"
//        - handle waves over periodic boundary conditions via grid coordinates

template <typename Float, size_t DIM>
class UniformGrid {
  enum { ON_MIN = -1, NOT_ON = 0, ON_MAX = 1 };

  std::vector<Cell<Float, DIM>> m_cells;
  size_t m_nx;
  size_t m_ny;

  Float m_x_min;
  Float m_x_max;
  Float m_y_min;
  Float m_y_max;
  Float m_dx;
  Float m_dy;

  std::vector<size_t> m_cut_cell_idxs;

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto to_vec_idx(size_t xi, size_t yi) const noexcept -> size_t {
    assert(xi < m_nx);
    assert(yi < m_ny);
    return yi * m_nx + xi;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto to_mat_idx(size_t i) const noexcept -> std::pair<size_t, size_t> {
    assert(i < m_cells.size());
    return std::make_pair(i % m_nx, i / m_nx);
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto to_grid_coord(const SimCoord<Float>& simulation_coord) const noexcept
      -> GridCoord<Float> {
    return GridCoord<Float>{
        .x = (simulation_coord.x - m_x_min) / m_dx,
        .y = (simulation_coord.y - m_y_min) / m_dy,
    };
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto to_sim_coord(const GridCoord<Float>& grid_coord) const noexcept
      -> SimCoord<Float> {
    return SimCoord<Float>{
        .x = grid_coord.x * m_dx + m_x_min,
        .y = grid_coord.y * m_dy + m_y_min,
    };
  }

 public:
  // -----------------------------------------------------------------------------------------------
  constexpr UniformGrid(
      Float x_min, Float x_max, size_t nx, Float y_min, Float y_max, size_t ny) noexcept
      : m_cells(nx * ny),
        m_nx(nx),
        m_ny(ny),
        m_x_min(x_min),
        m_x_max(x_max),
        m_y_min(y_min),
        m_y_max(y_max),
        m_dx((x_max - x_min) / static_cast<Float>(nx)),
        m_dy((y_max - y_min) / static_cast<Float>(ny)) {
    assert(nx > 0);
    assert(ny > 0);

    for (size_t yi = 0; yi < ny; ++yi) {
      for (size_t xi = 0; xi < nx; ++xi) {
        // clang-format off
        m_cells[to_vec_idx(xi, yi)] = Cell<Float, DIM>{
            .m_x_min    = m_x_min,
            .m_y_min    = m_y_min,
            .m_dx       = m_dx,
            .m_dy       = m_dy,
            .x_idx      = xi,
            .y_idx      = yi,
            .left_idx   = xi == 0        ? NULL_INDEX : to_vec_idx(xi - 1, yi),
            .right_idx  = xi == (nx - 1) ? NULL_INDEX : to_vec_idx(xi + 1, yi),
            .bottom_idx = yi == 0        ? NULL_INDEX : to_vec_idx(xi, yi - 1),
            .top_idx    = yi == (ny - 1) ? NULL_INDEX : to_vec_idx(xi, yi + 1),
        };
        // clang-format on
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void periodic_boundary(int sides = ALL) noexcept {
    // Left side
    if ((sides & LEFT) != 0) {
      for (size_t yi = 0; yi < m_ny; ++yi) {
        const auto left_idx        = to_vec_idx(0, yi);
        const auto right_idx       = to_vec_idx(m_nx - 1, yi);
        m_cells[left_idx].left_idx = right_idx;
      }
    }

    // Right side
    if ((sides & RIGHT) != 0) {
      for (size_t yi = 0; yi < m_ny; ++yi) {
        const auto left_idx          = to_vec_idx(0, yi);
        const auto right_idx         = to_vec_idx(m_nx - 1, yi);
        m_cells[right_idx].right_idx = left_idx;
      }
    }

    // Bottom side
    if ((sides & BOTTOM) != 0) {
      for (size_t xi = 0; xi < m_nx; ++xi) {
        const auto bottom_idx          = to_vec_idx(xi, 0);
        const auto top_idx             = to_vec_idx(xi, m_ny - 1);
        m_cells[bottom_idx].bottom_idx = top_idx;
      }
    }

    // Top side
    if ((sides & TOP) != 0) {
      for (size_t xi = 0; xi < m_nx; ++xi) {
        const auto bottom_idx    = to_vec_idx(xi, 0);
        const auto top_idx       = to_vec_idx(xi, m_ny - 1);
        m_cells[top_idx].top_idx = bottom_idx;
      }
    }
  }

 private:
  // -----------------------------------------------------------------------------------------------
  constexpr void set_boundary_impl(int sides, size_t idx) noexcept {
    // Left side
    if ((sides & LEFT) != 0) {
      for (size_t yi = 0; yi < m_ny; ++yi) {
        const auto left_idx        = to_vec_idx(0, yi);
        m_cells[left_idx].left_idx = idx;
      }
    }

    // Right side
    if ((sides & RIGHT) != 0) {
      for (size_t yi = 0; yi < m_ny; ++yi) {
        const auto right_idx         = to_vec_idx(m_nx - 1, yi);
        m_cells[right_idx].right_idx = idx;
      }
    }

    // Bottom side
    if ((sides & BOTTOM) != 0) {
      for (size_t xi = 0; xi < m_nx; ++xi) {
        const auto bottom_idx          = to_vec_idx(xi, 0);
        m_cells[bottom_idx].bottom_idx = idx;
      }
    }

    // Top side
    if ((sides & TOP) != 0) {
      for (size_t xi = 0; xi < m_nx; ++xi) {
        const auto top_idx       = to_vec_idx(xi, m_ny - 1);
        m_cells[top_idx].top_idx = idx;
      }
    }
  }

 public:
  // -----------------------------------------------------------------------------------------------
  constexpr void zero_flux_boundary(int sides = ALL) noexcept {
    set_boundary_impl(sides, ZERO_FLUX_INDEX);
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void same_value_boundary(int sides = ALL) noexcept {
    set_boundary_impl(sides, SAME_VALUE_INDEX);
  }

  // -----------------------------------------------------------------------------------------------
  template <typename FUNC>
  constexpr void fill_center(FUNC f) noexcept {
    const auto get_centroid = [](const auto& points) -> Eigen::Vector<Float, 2> {
      Eigen::Vector<Float, 2> centroid = Eigen::Vector<Float, 2>::Zero();
      for (const auto& p : points) {
        centroid += p;
      }
      return centroid / points.size();
    };

    for (auto& cell : m_cells) {
      if (cell.is_cartesian()) {
        auto& value = cell.get_cartesian().value;
        if constexpr (DIM == 1) {
          value(0) =
              f(cell.template x_min<SIM_C>() + static_cast<Float>(0.5) * cell.template dx<SIM_C>(),
                cell.template y_min<SIM_C>() + static_cast<Float>(0.5) * cell.template dy<SIM_C>());
        } else {
          value =
              f(cell.template x_min<SIM_C>() + static_cast<Float>(0.5) * cell.template dx<SIM_C>(),
                cell.template y_min<SIM_C>() + static_cast<Float>(0.5) * cell.template dy<SIM_C>());
        }
      } else if (cell.is_cut()) {
        Eigen::Vector<Float, 2> left_centroid =
            get_centroid(get_left_points<Cell<Float, DIM>, Float>(cell));
        Eigen::Vector<Float, 2> right_centroid =
            get_centroid(get_right_points<Cell<Float, DIM>, Float>(cell));

        auto& cell_value = cell.get_cut();
        if constexpr (DIM == 1) {
          cell_value.left_value(0)  = f(left_centroid(X), left_centroid(Y));
          cell_value.right_value(0) = f(right_centroid(X), right_centroid(Y));
        } else {
          cell_value.left_value  = f(left_centroid(X), left_centroid(Y));
          cell_value.right_value = f(right_centroid(X), right_centroid(Y));
        }
      } else {
        Igor::Panic("Unknown cell type with variant index {}.", cell.value.index());
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <typename FUNC>
  constexpr void fill_four_point(FUNC f) noexcept {
    for (auto& cell : m_cells) {
      if (cell.is_cartesian()) {
        auto& value = cell.get_cartesian().value;
        if constexpr (DIM == 1) {
          value(0) = (f(cell.template x_min<SIM_C>() + cell.template dx<SIM_C>() / 4,
                        cell.template y_min<SIM_C>() + cell.template dy<SIM_C>() / 4) +
                      f(cell.template x_min<SIM_C>() + 3 * cell.template dx<SIM_C>() / 4,
                        cell.template y_min<SIM_C>() + cell.template dy<SIM_C>() / 4) +
                      f(cell.template x_min<SIM_C>() + cell.template dx<SIM_C>() / 4,
                        cell.template y_min<SIM_C>() + 3 * cell.template dy<SIM_C>() / 4) +
                      f(cell.template x_min<SIM_C>() + 3 * cell.template dx<SIM_C>() / 4,
                        cell.template y_min<SIM_C>() + 3 * cell.template dy<SIM_C>() / 4)) /
                     4;
        } else {
          value = (f(cell.template x_min<SIM_C>() + cell.template dx<SIM_C>() / 4,
                     cell.template y_min<SIM_C>() + cell.template dy<SIM_C>() / 4) +
                   f(cell.template x_min<SIM_C>() + 3 * cell.template dx<SIM_C>() / 4,
                     cell.template y_min<SIM_C>() + cell.template dy<SIM_C>() / 4) +
                   f(cell.template x_min<SIM_C>() + cell.template dx<SIM_C>() / 4,
                     cell.template y_min<SIM_C>() + 3 * cell.template dy<SIM_C>() / 4) +
                   f(cell.template x_min<SIM_C>() + 3 * cell.template dx<SIM_C>() / 4,
                     cell.template y_min<SIM_C>() + 3 * cell.template dy<SIM_C>() / 4)) /
                  4;
        }
      } else if (cell.is_cut()) {
        auto& cut_value = cell.get_cut();

        // Left side
        {
          const auto corners = cell.template get_left_points<SIM_C>();
          for (size_t j = 0; j < corners.size(); ++j) {
            auto point = SimCoord<Float>::Zero();
            for (size_t i = 0; i < corners.size(); ++i) {
              point += (1 + (i == j)) * corners[i] / static_cast<Float>(corners.size() + 1);
            }

            if constexpr (DIM == 1) {
              cut_value.left_value(0) += f(point.x, point.y);
            } else {
              cut_value.left_value += f(point.x, point.y);
            }
          }
          cut_value.left_value /= static_cast<Float>(corners.size());
        }

        // Right side
        {
          const auto corners = cell.template get_right_points<SIM_C>();
          for (size_t j = 0; j < corners.size(); ++j) {
            auto point = SimCoord<Float>::Zero();
            for (size_t i = 0; i < corners.size(); ++i) {
              point += (1 + (i == j)) * corners[i] / static_cast<Float>(corners.size() + 1);
            }

            if constexpr (DIM == 1) {
              cut_value.right_value(0) += f(point.x, point.y);
            } else {
              cut_value.right_value += f(point.x, point.y);
            }
          }
          cut_value.right_value /= static_cast<Float>(corners.size());
        }
      } else {
        Igor::Panic("Unknown cell type with variant index {}.", cell.value.index());
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <Point2D_c PointType>
  [[nodiscard]] constexpr auto classify_cut(const Cell<Float, DIM>& cell,
                                            PointType& cut1_point,
                                            PointType& cut2_point) const noexcept -> CutType {
    constexpr CoordType coord_type = PointType2CoordType<PointType>;
    auto on_single_side            = [](Side side) -> bool {
      return static_cast<int>((side & LEFT) > 0) + static_cast<int>((side & BOTTOM) > 0) +
                 static_cast<int>((side & RIGHT) > 0) + static_cast<int>((side & TOP) > 0) ==
             1;
    };
    auto remove_common = [](Side& side1, Side& side2) {
      const Side common = static_cast<Side>(side1 & side2);
      side1             = static_cast<Side>(side1 ^ common);
      side2             = static_cast<Side>(side2 ^ common);
    };

    const int cut1_on_x =
        approx_eq(cell.template x_min<coord_type>(), cut1_point.x) * ON_MIN +
        approx_eq(cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                  cut1_point.x) *
            ON_MAX;
    const int cut1_on_y =
        approx_eq(cell.template y_min<coord_type>(), cut1_point.y) * ON_MIN +
        approx_eq(cell.template y_min<coord_type>() + cell.template dy<coord_type>(),
                  cut1_point.y) *
            ON_MAX;
    assert(cut1_on_x != NOT_ON || cut1_on_y != NOT_ON);
    Side cut1_loc = static_cast<Side>((cut1_on_x == ON_MIN) * LEFT | (cut1_on_x == ON_MAX) * RIGHT |
                                      (cut1_on_y == ON_MIN) * BOTTOM | (cut1_on_y == ON_MAX) * TOP);

    const int cut2_on_x =
        approx_eq(cell.template x_min<coord_type>(), cut2_point.x) * ON_MIN +
        approx_eq(cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                  cut2_point.x) *
            ON_MAX;
    const int cut2_on_y =
        approx_eq(cell.template y_min<coord_type>(), cut2_point.y) * ON_MIN +
        approx_eq(cell.template y_min<coord_type>() + cell.template dy<coord_type>(),
                  cut2_point.y) *
            ON_MAX;
    assert(cut2_on_x != NOT_ON || cut2_on_y != NOT_ON);
    Side cut2_loc = static_cast<Side>((cut2_on_x == ON_MIN) * LEFT | (cut2_on_x == ON_MAX) * RIGHT |
                                      (cut2_on_y == ON_MIN) * BOTTOM | (cut2_on_y == ON_MAX) * TOP);

    IGOR_ASSERT(cut1_loc != cut2_loc,
                "Expected the cuts to be at different sides of the cell, but both cuts are on "
                "side {}. cut1_point = {}, cut2_point = {}. PointType = {}.",
                cut1_loc,
                cut1_point,
                cut2_point,
                Igor::type_name<PointType>());

    remove_common(cut1_loc, cut2_loc);

    IGOR_ASSERT(on_single_side(cut1_loc) && on_single_side(cut2_loc),
                "Expected cut1 on cut2 to be on only one side after removing common sides, but "
                "cut1_loc = {} and cut2_loc = {}.",
                cut1_loc,
                cut2_loc);

    IGOR_ASSERT(cut1_loc != cut2_loc,
                "Expected the cuts to be at different sides of the cell, but both cuts are on "
                "side {} (== {}).",
                cut1_loc,
                cut2_loc);

    if (cut1_loc > cut2_loc) {
      std::swap(cut1_loc, cut2_loc);
      std::swap(cut1_point, cut2_point);
    }

    CutType type;
    switch (cut1_loc | cut2_loc) {
      case LEFT | BOTTOM:  type = CutType::BOTTOM_LEFT; break;
      case BOTTOM | RIGHT: type = CutType::BOTTOM_RIGHT; break;
      case RIGHT | TOP:    type = CutType::TOP_RIGHT; break;
      case TOP | LEFT:     type = CutType::TOP_LEFT; break;
      case LEFT | RIGHT:   type = CutType::MIDDLE_HORI; break;
      case BOTTOM | TOP:   type = CutType::MIDDLE_VERT; break;
      default:
        Igor::Panic("Invalid combination cut1_loc = {} and cut2_loc = {}",
                    static_cast<std::underlying_type_t<Side>>(cut1_loc),
                    static_cast<std::underlying_type_t<Side>>(cut2_loc));
        std::unreachable();
    }

    // Make cut relative to cell size, i.e. cut in [0, 1]
    cut1_point.x =
        (cut1_point.x - cell.template x_min<coord_type>()) / cell.template dx<coord_type>();
    cut1_point.y =
        (cut1_point.y - cell.template y_min<coord_type>()) / cell.template dy<coord_type>();
    cut2_point.x =
        (cut2_point.x - cell.template x_min<coord_type>()) / cell.template dx<coord_type>();
    cut2_point.y =
        (cut2_point.y - cell.template y_min<coord_type>()) / cell.template dy<coord_type>();

    IGOR_ASSERT(cut1_point.x >= -EPS<Float> && cut1_point.x - 1 <= EPS<Float>,
                "Cut1-x {} is not in [0, 1].",
                cut1_point.x);
    IGOR_ASSERT(cut1_point.y >= -EPS<Float> && cut1_point.y - 1 <= EPS<Float>,
                "Cut1-y {} is not in [0, 1].",
                cut1_point.y);
    IGOR_ASSERT(cut2_point.x >= -EPS<Float> && cut2_point.x - 1 <= EPS<Float>,
                "Cut2-x {} is not in [0, 1].",
                cut2_point.x);
    IGOR_ASSERT(cut2_point.y >= -EPS<Float> && cut2_point.y - 1 <= EPS<Float>,
                "Cut2-y {} is not in [0, 1].",
                cut2_point.y);

    return type;
  }

  // -----------------------------------------------------------------------------------------------
  template <Point2D_c PointType>
  [[nodiscard]] constexpr auto find_next_cell_to_cut(const Cell<Float, DIM>& cell,
                                                     const PointType& exit_point) const noexcept
      -> size_t {
    IGOR_ASSERT(point_in_cell(exit_point, cell), "Point {} is not in cell {}", exit_point, cell);

    constexpr CoordType coord_type = PointType2CoordType<PointType>;

    const int on_x = approx_eq(cell.template x_min<coord_type>(), exit_point.x) * ON_MIN +
                     approx_eq(cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                               exit_point.x) *
                         ON_MAX;
    const int on_y = approx_eq(cell.template y_min<coord_type>(), exit_point.y) * ON_MIN +
                     approx_eq(cell.template y_min<coord_type>() + cell.template dy<coord_type>(),
                               exit_point.y) *
                         ON_MAX;
    assert(on_x != NOT_ON || on_y != NOT_ON);

    if (on_x != NOT_ON && on_y != NOT_ON) {
      // Bottom left corner
      if (on_x == ON_MIN && on_y == ON_MIN) {
        assert(is_cell(cell.bottom_idx));
        assert(is_cell(m_cells[cell.bottom_idx].left_idx));

        // TODO: Add check that this is actually the correct next cell
        Igor::Todo("Exit on bottom left corner.");
        return m_cells[cell.bottom_idx].left_idx;
      }
      // Top left corner
      else if (on_x == ON_MIN && on_y == ON_MAX) {
        assert(is_cell(cell.top_idx));
        assert(is_cell(m_cells[cell.top_idx].left_idx));

        // TODO: Add check that this is actually the correct next cell
        Igor::Todo("Exit on top left corner.");
        return m_cells[cell.top_idx].left_idx;
      }
      // Bottom right corner
      else if (on_x == ON_MAX && on_y == ON_MIN) {
        assert(is_cell(cell.bottom_idx));
        assert(is_cell(m_cells[cell.bottom_idx].right_idx));

        // TODO: Add check that this is actually the correct next cell
        Igor::Todo("Exit on bottom right corner.");
        return m_cells[cell.bottom_idx].right_idx;
      }
      // Top right corner
      else {
        assert(is_cell(cell.top_idx));
        assert(is_cell(m_cells[cell.top_idx].right_idx));

        // TODO: Add check that this is actually the correct next cell
        // Igor::Todo("Exit on top right corner.");
        // return m_cells[cell.top_idx].right_idx;

        return cell.top_idx;
      }
    } else {
      // Left side
      if (on_x == ON_MIN) {
        assert(is_cell(cell.left_idx));
        return cell.left_idx;
      }
      // Right side
      else if (on_x == ON_MAX) {
        assert(is_cell(cell.right_idx));
        return cell.right_idx;
      }
      // Bottom side
      else if (on_y == ON_MIN) {
        assert(is_cell(cell.bottom_idx));
        return cell.bottom_idx;
      }
      // Top side
      else {
        assert(is_cell(cell.top_idx));
        return cell.top_idx;
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <typename ParamCurve>
  [[nodiscard]] constexpr auto cut_curve(ParamCurve curve) noexcept -> bool {
    enum { T_BELOW_EXIT = -1, T_IS_EXIT, T_ABOVE_EXIT };

    using PointType = GridCoord<Float>;

    const auto t_location = [this,
                             curve](Float t, const Cell<Float, DIM>& cell, Float t_entry) -> int {
      const auto pos = to_grid_coord(curve(t));
      if (t <= t_entry || (point_in_cell(pos, cell) && !point_on_boundary(pos, cell))) {
        return T_BELOW_EXIT;
      } else if (point_on_boundary(pos, cell)) {
        return T_IS_EXIT;
      }
      return T_ABOVE_EXIT;
    };

    Float t             = 0;
    const auto init_pos = to_grid_coord(curve(t));
    auto cell_idx       = find_cell(init_pos);
    if (cell_idx == NULL_INDEX) {
      Igor::Warn("Point {} on parametrized curve is not in the grid.", curve(t));
      return false;
    }

    if (!m_cut_cell_idxs.empty()) {
      Igor::Warn("Grid is already cut, this method expects the grid to no be cut.");
      return false;
    }

    m_cut_cell_idxs.push_back(cell_idx);
    std::vector<PointType> entry_points{};
    entry_points.push_back(init_pos);

    while (t < 1) {
      const auto& cell = m_cells[cell_idx];
      Float t_lower    = t;
      Float t_upper    = std::min(Float{1}, t + Float{0.5});

      bool found_exit = false;
      while (!found_exit) {
        int t_lower_loc = t_location(t_lower, cell, t);
        int t_upper_loc = t_location(t_upper, cell, t);
        if (t_lower_loc == T_IS_EXIT && t_upper_loc == T_IS_EXIT) {
          Igor::Warn("lower t ({}) and upper t ({}) are both exit points.", t_lower, t_upper);
        }
        assert(!(t_lower_loc == T_IS_EXIT && t_upper_loc == T_IS_EXIT));

        if (!(t_lower_loc == T_BELOW_EXIT || t_lower_loc == T_IS_EXIT)) {
          Igor::Warn("lower t ({}) is outside of cell.", t_lower);
        }
        assert(t_lower_loc == T_BELOW_EXIT || t_lower_loc == T_IS_EXIT);

        if (!(t_upper_loc == T_ABOVE_EXIT || t_upper_loc == T_IS_EXIT)) {
          IGOR_DEBUG_PRINT(t);
          IGOR_DEBUG_PRINT(t_lower);
          IGOR_DEBUG_PRINT(t_upper);
          IGOR_DEBUG_PRINT(curve(t_upper));
          IGOR_DEBUG_PRINT(to_grid_coord(curve(t_upper)));
          IGOR_DEBUG_PRINT(point_on_boundary(to_grid_coord(curve(t_upper)), cell));
          Igor::Warn("upper t ({}) is inside of cell.", t_upper);
        }
        assert(t_upper_loc == T_ABOVE_EXIT || t_upper_loc == T_IS_EXIT);

        if (t_lower_loc == T_IS_EXIT) {
          t          = t_lower;
          found_exit = true;
        } else if (t_upper_loc == T_IS_EXIT) {
          t          = t_upper;
          found_exit = true;
        } else {
          const auto t_mid = (t_upper + t_lower) / 2;
          switch (t_location(t_mid, cell, t)) {
            case T_BELOW_EXIT: t_lower = t_mid; break;
            case T_ABOVE_EXIT: t_upper = t_mid; break;
            case T_IS_EXIT:
              t          = t_mid;
              found_exit = true;
              break;
          }
        }
      }
      const auto next_pos = to_grid_coord(curve(t));
      entry_points.push_back(next_pos);

      if (approx_eq(t, Float{1})) { break; }

      cell_idx = find_next_cell_to_cut(cell, next_pos);
      m_cut_cell_idxs.push_back(cell_idx);
    }

    assert(m_cut_cell_idxs.size() + 1 == entry_points.size());
    for (size_t i = 0; i < m_cut_cell_idxs.size(); ++i) {
      auto& cell_to_cut = m_cells[m_cut_cell_idxs[i]];
      auto cut1_point   = entry_points[i];
      auto cut2_point   = entry_points[i + 1];

      const auto type = classify_cut(cell_to_cut, cut1_point, cut2_point);

      Eigen::Vector<Float, DIM> old_value = Eigen::Vector<Float, DIM>::Zero();
      if (cell_to_cut.is_cartesian()) { old_value = cell_to_cut.get_cartesian().value; }

      cell_to_cut.value = CutValue<Float, DIM>{
          .left_value  = old_value,
          .right_value = old_value,
          .type        = type,
          .rel_cut1    = {cut1_point.x, cut1_point.y},
          .rel_cut2    = {cut2_point.x, cut2_point.y},
      };
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  template <ExtendType extend_type = ExtendType::MAX>
  [[nodiscard]] constexpr auto
  cut_piecewise_linear(const std::vector<SimCoord<Float>>& points) noexcept -> bool {
    std::vector<GridCoord<Float>> points_grid_coord(points.size());
    std::transform(std::cbegin(points),
                   std::cend(points),
                   std::begin(points_grid_coord),
                   [this](const SimCoord<Float>& p) { return to_grid_coord(p); });
    return cut_piecewise_linear<extend_type>(points_grid_coord);
  }

  // -----------------------------------------------------------------------------------------------
  template <ExtendType extend_type = ExtendType::MAX>
  [[nodiscard]] constexpr auto cut_piecewise_linear(std::vector<GridCoord<Float>> points) noexcept
      -> bool {
    if (!m_cut_cell_idxs.empty()) {
      Igor::Warn("Grid is already cut, this method expects the grid to no be cut.");
      return false;
    }

    auto get_whole_min = [](Float x1, Float x2) -> int {
      return static_cast<int>(std::ceil(std::min(x1, x2)));
    };
    auto get_whole_max = [](Float x1, Float x2) -> int {
      return static_cast<int>(std::floor(std::max(x1, x2)));
    };

    // Remove points outside of grid
    std::erase_if(points, [this](const GridCoord<Float>& p) { return !point_in_grid(p); });
    IGOR_ASSERT(points.size() >= 2,
                "Only {} points inside of the grid, require at least 2.",
                points.size());

    // TODO: Add option to extend to nearest intersection point
    // - Extend end points -------------------------------------------------------------------------
    if constexpr (extend_type != ExtendType::NONE) {
      // - Handle first point; extend curve backwards to cut cell containing the first point -------
      {
        auto& p0       = points[0];
        const auto& p1 = points[1];

        const GridCoord<Float> s = p1 - p0;
        const std::array<Float, 4> rs{
            (std::floor(p0.x) - p0.x) / s.x,
            (std::ceil(p0.x) - p0.x) / s.x,
            (std::floor(p0.y) - p0.y) / s.y,
            (std::ceil(p0.y) - p0.y) / s.y,
        };

        Float r_entry = -std::numeric_limits<Float>::max();
        for (Float r : rs) {
          if constexpr (extend_type == ExtendType::MAX) {
            if (r < EPS<Float> && r > r_entry) { r_entry = r; }
          } else {
            if (std::abs(r) < std::abs(r_entry)) { r_entry = r; }
          }
        }
        IGOR_ASSERT(r_entry != -std::numeric_limits<Float>::max(),
                    "Could not find the extended intersection point for p0={} and p1={}. Potential "
                    "r-values are rs={}",
                    p0,
                    p1,
                    rs);
        p0 += r_entry * s;
      }

      // - Handle last point; extend curve forwards to cut cell containing the last point ----------
      {
        const auto& p0 = points[points.size() - 2];
        auto& p1       = points[points.size() - 1];

        const GridCoord<Float> s = p1 - p0;
        const std::array<Float, 4> rs{
            (std::floor(p1.x) - p1.x) / s.x,
            (std::ceil(p1.x) - p1.x) / s.x,
            (std::floor(p1.y) - p1.y) / s.y,
            (std::ceil(p1.y) - p1.y) / s.y,
        };

        Float r_exit = std::numeric_limits<Float>::max();
        for (Float r : rs) {
          if constexpr (extend_type == ExtendType::MAX) {
            if (r > -EPS<Float> && r < r_exit) { r_exit = r; }
          } else {
            if (std::abs(r) < std::abs(r_exit)) { r_exit = r; }
          }
        }
        IGOR_ASSERT(r_exit != std::numeric_limits<Float>::max(),
                    "Could not find the extended intersection point for p0={} and p1={}. Potential "
                    "r-values are rs={}",
                    p0,
                    p1,
                    rs);

        p1 += r_exit * s;
      }
    }
    // - Extend end points -------------------------------------------------------------------------

    std::vector<GridCoord<Float>> intersect_points{};
    ptrdiff_t begin_points_idx = 0;
    for (size_t i = 0; i < points.size() - 1; ++i) {
      const auto& p0 = points[i];
      const auto& p1 = points[i + 1];

      const auto slope = p1 - p0;

      const bool search_whole_x = !approx_eq(slope.x, static_cast<Float>(0));
      const bool search_whole_y = !approx_eq(slope.y, static_cast<Float>(0));
      if (!(search_whole_x || search_whole_y)) {
        Igor::Warn("p0 = {} and p1 = {} are the same point.", p0, p1);
        continue;
      }

      // Intersections where x in N
      if (search_whole_x) {
        for (int x = get_whole_min(p0.x, p1.x); x <= get_whole_max(p0.x, p1.x); ++x) {
          const Float r = (static_cast<Float>(x) - p0.x) / slope.x;
          const Float y = p0.y + r * slope.y;
          intersect_points.emplace_back(static_cast<Float>(x), y);
        }
      }

      // Intersections where y in N
      if (search_whole_y) {
        for (int y = get_whole_min(p0.y, p1.y); y <= get_whole_max(p0.y, p1.y); ++y) {
          const Float r = (static_cast<Float>(y) - p0.y) / slope.y;
          const Float x = p0.x + r * slope.x;

          // If x is a whole number, the point was already added in the previous loop, can ignore it
          if (search_whole_x && approx_eq(x, std::round(x))) { continue; }

          intersect_points.emplace_back(x, static_cast<Float>(y));
        }
      }

      std::sort(std::next(std::begin(intersect_points), begin_points_idx),
                std::end(intersect_points),
                [&p0](const GridCoord<Float>& lhs, const GridCoord<Float>& rhs) {
                  const auto lhs_dist = (lhs - p0).norm();
                  const auto rhs_dist = (rhs - p0).norm();
                  return lhs_dist < rhs_dist;
                });
      begin_points_idx = static_cast<ptrdiff_t>(intersect_points.size());
    }

    // Remove duplicate points
    // TODO: Can we do this on the fly while filling intersect_points?
    intersect_points.erase(
        std::unique(intersect_points.begin(),
                    intersect_points.end(),
                    [](const GridCoord<Float>& lhs, const GridCoord<Float>& rhs) {
                      return (lhs - rhs).norm() < EPS<Float>;
                    }),
        intersect_points.end());

    m_cut_cell_idxs.resize(intersect_points.size() - 1);
    for (size_t i = 0; i < intersect_points.size() - 1; ++i) {
      const auto entry_point = intersect_points[i];
      const auto exit_point  = intersect_points[i + 1];
      if (approx_eq((entry_point - exit_point).norm(), static_cast<Float>(0))) {
        Igor::Warn("entry_point {} and exit_point {} are the same point, ignore.",
                   entry_point,
                   exit_point);
        continue;
      }

      const auto mid_point = (entry_point + exit_point) / 2;
      const auto cell_idx  = find_cell(mid_point);
      IGOR_ASSERT(cell_idx != NULL_INDEX,
                  "Expected mid_point {} to be in grid [{}, {}]x[{}, {}], but is not.",
                  mid_point,
                  0,
                  m_nx,
                  0,
                  m_ny);
      m_cut_cell_idxs[i] = cell_idx;
    }

    IGOR_ASSERT(
        m_cut_cell_idxs.size() + 1 == intersect_points.size(),
        "Expected to find exactly one point more than cells, but found {} points and {} cells",
        intersect_points.size(),
        m_cut_cell_idxs.size());
    for (size_t i = 0; i < m_cut_cell_idxs.size(); ++i) {
      auto& cell_to_cut = m_cells[m_cut_cell_idxs[i]];
      auto cut1_point   = intersect_points[i];
      auto cut2_point   = intersect_points[i + 1];

      const auto type = classify_cut(cell_to_cut, cut1_point, cut2_point);

      Eigen::Vector<Float, DIM> old_value = Eigen::Vector<Float, DIM>::Zero();
      if (cell_to_cut.is_cartesian()) { old_value = cell_to_cut.get_cartesian().value; }

      cell_to_cut.value = CutValue<Float, DIM>{
          .left_value  = old_value,
          .right_value = old_value,
          .type        = type,
          .rel_cut1    = {cut1_point.x, cut1_point.y},
          .rel_cut2    = {cut2_point.x, cut2_point.y},
      };
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void merge_cut_cells() noexcept {
    for (size_t cell_idx : m_cut_cell_idxs) {
      assert(is_cell(cell_idx));
      auto& cell = m_cells[cell_idx];
      assert(cell.is_cut());

      const auto left_value   = cell.get_cut().left_value;
      const auto left_polygon = cell.template get_cut_left_polygon<SIM_C>();

      const auto right_value   = cell.get_cut().right_value;
      const auto right_polygon = cell.template get_cut_right_polygon<SIM_C>();

      const auto cell_polygon = cell.template get_cartesian_polygon<SIM_C>();

      cell.value = CartesianValue<Float, DIM>{
          .value = (left_value * left_polygon.area() + right_value * right_polygon.area()) /
                   cell_polygon.area()};
    }
    m_cut_cell_idxs.clear();
  }

  // -----------------------------------------------------------------------------------------------
  template <typename PointType>
  [[nodiscard]] constexpr auto point_in_grid(const PointType& point) const noexcept -> bool {
    if constexpr (std::is_same_v<std::remove_cvref_t<PointType>, SimCoord<Float>>) {
      return point.x - m_x_min >= -EPS<Float> && m_x_max - point.x >= -EPS<Float> &&  //
             point.y - m_y_min >= -EPS<Float> && m_y_max - point.y >= -EPS<Float>;
    } else {
      return point.x >= -EPS<Float> && static_cast<Float>(m_nx) - point.x >= -EPS<Float> &&  //
             point.y >= -EPS<Float> && static_cast<Float>(m_ny) - point.y >= -EPS<Float>;
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <Point2D_c PointType>
  [[nodiscard]] constexpr auto find_cell(const PointType& pos) const noexcept -> size_t {
    const GridCoord<Float> gpos = [this, &pos] {
      if constexpr (std::is_same_v<std::remove_cvref_t<PointType>, SimCoord<Float>>) {
        return to_grid_coord(pos);
      } else {
        (void)this;
        return pos;
      }
    }();

    if (gpos.x < 0 || gpos.x > static_cast<Float>(m_nx)) {
      Igor::Warn("x-component of position {} is not in grid with nx={}", pos, m_nx);
      return NULL_INDEX;
    }
    if (gpos.y < 0 || gpos.y > static_cast<Float>(m_ny)) {
      Igor::Warn("y-component of position {} is not in grid with ny={}", pos, m_ny);
      return NULL_INDEX;
    }

    const size_t xi = [this, &gpos] {
      if (approx_eq(gpos.x, static_cast<Float>(m_nx))) {
        return m_nx - 1;
      } else {
        return static_cast<size_t>(std::floor(gpos.x));
      }
    }();

    const size_t yi = [this, &gpos] {
      if (approx_eq(gpos.y, static_cast<Float>(m_ny))) {
        return m_ny - 1;
      } else {
        return static_cast<size_t>(std::floor(gpos.y));
      }
    }();

    IGOR_ASSERT(xi < m_nx,
                "gpos.x={:.16f} translates to xi={} which is out of bounds for nx={}.",
                gpos.x,
                xi,
                m_nx);
    IGOR_ASSERT(yi < m_nx,
                "gpos.y={:.16f} translates to yi={} which is out of bounds for ny={}.",
                gpos.y,
                yi,
                m_ny);

    return to_vec_idx(xi, yi);
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto get_shock_curve() const noexcept -> std::vector<SimCoord<Float>> {
    if (m_cut_cell_idxs.empty()) { return {}; }

    std::vector<SimCoord<Float>> curve_points(m_cut_cell_idxs.size() + 1);

    {
      const auto& cell = m_cells[m_cut_cell_idxs.front()];
      curve_points[0]  = cell.template cut1<SIM_C>();
    }

    for (size_t i = 0; i < m_cut_cell_idxs.size(); ++i) {
      const auto& cell    = m_cells[m_cut_cell_idxs[i]];
      curve_points[i + 1] = cell.template cut2<SIM_C>();
    }

    return curve_points;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto get_existing_neighbours(size_t idx) const noexcept
      -> SmallVector<size_t> {
    if (!is_cell(idx)) { return {}; }

    const auto& center_cell = m_cells[idx];
    SmallVector<size_t> neighbours{};

    // Left cell
    if (is_cell(center_cell.left_idx)) { neighbours.push_back(center_cell.left_idx); }
    // Bottom cell
    if (is_cell(center_cell.bottom_idx)) { neighbours.push_back(center_cell.bottom_idx); }
    // Right cell
    if (is_cell(center_cell.right_idx)) { neighbours.push_back(center_cell.right_idx); }
    // Top cell
    if (is_cell(center_cell.top_idx)) { neighbours.push_back(center_cell.top_idx); }

    // Bottom left diagonal cell
    if (is_cell(center_cell.left_idx) && is_cell(center_cell.bottom_idx)) {
      assert(m_cells[center_cell.left_idx].bottom_idx == m_cells[center_cell.bottom_idx].left_idx);
      if (is_cell(m_cells[center_cell.left_idx].bottom_idx)) {
        neighbours.push_back(m_cells[center_cell.left_idx].bottom_idx);
      }
    }
    // Bottom right diagonal cell
    if (is_cell(center_cell.right_idx) && is_cell(center_cell.bottom_idx)) {
      assert(m_cells[center_cell.right_idx].bottom_idx ==
             m_cells[center_cell.bottom_idx].right_idx);
      if (is_cell(m_cells[center_cell.right_idx].bottom_idx)) {
        neighbours.push_back(m_cells[center_cell.right_idx].bottom_idx);
      }
    }
    // Top left diagonal cell
    if (is_cell(center_cell.left_idx) && is_cell(center_cell.top_idx)) {
      assert(m_cells[center_cell.left_idx].top_idx == m_cells[center_cell.top_idx].left_idx);
      if (is_cell(m_cells[center_cell.left_idx].top_idx)) {
        neighbours.push_back(m_cells[center_cell.left_idx].top_idx);
      }
    }
    // Top right diagonal cell
    if (is_cell(center_cell.right_idx) && is_cell(center_cell.top_idx)) {
      assert(m_cells[center_cell.right_idx].top_idx == m_cells[center_cell.top_idx].right_idx);
      if (is_cell(m_cells[center_cell.right_idx].top_idx)) {
        neighbours.push_back(m_cells[center_cell.right_idx].top_idx);
      }
    }

    return neighbours;
  }

  // -----------------------------------------------------------------------------------------------
  template <Point2D_c PointType>
  [[nodiscard]] constexpr auto eval(const PointType& point) const noexcept
      -> Eigen::Vector<Float, DIM> {
    const auto cell_idx = find_cell(point);
    if (cell_idx == NULL_INDEX) {
      Igor::Panic(
          "Point {} is not in grid [{}, {}]x[{}, {}].", point, m_x_min, m_x_max, m_y_min, m_y_max);
    }

    const auto& cell = m_cells[cell_idx];
    if (cell.is_cartesian()) {
      return cell.get_cartesian().value;
    } else if (cell.is_cut()) {
      constexpr CoordType coord_type = PointType2CoordType<PointType>;

      const auto point_in_left =
          cell.template get_cut_left_polygon<coord_type>().point_in_polygon(point);
      const auto point_in_right =
          cell.template get_cut_right_polygon<coord_type>().point_in_polygon(point);

      if (point_in_left) {
        return cell.get_cut().left_value;
      } else if (point_in_right) {
        return cell.get_cut().right_value;
      } else {
#ifdef ZAP_WARN_GRID_EVAL
        if (point_in_left && point_in_right) {
          Igor::Warn("Point {} is in both subcells, take the average.", point);
        } else if (!(point_in_left || point_in_right)) {
          Igor::Warn("Point {} is in neither subcell, probably a rounding error, take the average.",
                     point);
        }
#endif  // ZAP_WARN_GRID_EVAL
        return (cell.get_cut().left_value + cell.get_cut().right_value) / 2;
      }
    } else {
      Igor::Panic("Unknown cell type with variant index {}", cell.value.index());
      std::unreachable();
    }
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto nx() const noexcept -> size_t { return m_nx; }
  [[nodiscard]] constexpr auto ny() const noexcept -> size_t { return m_ny; }
  [[nodiscard]] constexpr auto min_delta() const noexcept -> Float { return std::min(m_dx, m_dy); }
  [[nodiscard]] constexpr auto x_min() const noexcept -> Float { return m_x_min; }
  [[nodiscard]] constexpr auto x_max() const noexcept -> Float { return m_x_max; }
  [[nodiscard]] constexpr auto y_min() const noexcept -> Float { return m_y_min; }
  [[nodiscard]] constexpr auto y_max() const noexcept -> Float { return m_y_max; }
  [[nodiscard]] constexpr auto scale_x() const noexcept -> Float { return 1 / m_dx; }
  [[nodiscard]] constexpr auto scale_y() const noexcept -> Float { return 1 / m_dy; }
  [[nodiscard]] constexpr auto cells() noexcept -> std::vector<Cell<Float, DIM>>& {
    return m_cells;
  }
  [[nodiscard]] constexpr auto cells() const noexcept -> const std::vector<Cell<Float, DIM>>& {
    return m_cells;
  }
  [[nodiscard]] constexpr auto cut_cell_idxs() const noexcept -> const std::vector<size_t>& {
    return m_cut_cell_idxs;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto size() const noexcept -> size_t { return m_cells.size(); }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto operator[](size_t idx) const noexcept -> const Cell<Float, DIM>& {
    assert(idx < m_nx * m_ny);
    return m_cells[idx];
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto operator[](size_t idx) noexcept -> Cell<Float, DIM>& {
    assert(idx < m_nx * m_ny);
    return m_cells[idx];
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto begin() noexcept { return m_cells.begin(); }
  [[nodiscard]] constexpr auto begin() const noexcept { return m_cells.begin(); }
  [[nodiscard]] constexpr auto cbegin() const noexcept { return m_cells.cbegin(); }
  [[nodiscard]] constexpr auto end() noexcept { return m_cells.end(); }
  [[nodiscard]] constexpr auto end() const noexcept { return m_cells.end(); }
  [[nodiscard]] constexpr auto cend() const noexcept { return m_cells.cend(); }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto is_cell(size_t idx) const noexcept -> bool {
    return idx < m_cells.size();
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void dump_cells(std::ostream& out) const noexcept {
    for (size_t i = 0; i < m_cells.size(); ++i) {
      const auto& cell = m_cells[i];
      out << i << ": " << cell;
      if (i + 1 == m_cells.size()) {
        out << '\n';
      } else {
        out << ",\n";
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  friend struct Cell<Float, DIM>;
};

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_GRID_HPP_
