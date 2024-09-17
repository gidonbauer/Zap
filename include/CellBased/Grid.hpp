#ifndef ZAP_CELL_BASED_GRID_HPP_
#define ZAP_CELL_BASED_GRID_HPP_

#include <numeric>

#include "CellBased/Cell.hpp"
#include "CellBased/Definitions.hpp"
#include "CellBased/GridHelper.hpp"
#include "IO/Fwd.hpp"

namespace Zap::CellBased {

template <typename Float, size_t DIM>
class UniformGrid {
  enum { ON_MIN = -1, NOT_ON = 0, ON_MAX = 1 };

  struct m_CornerIndices {
    size_t top_left;
    size_t top_right;
    size_t bottom_left;
    size_t bottom_right;
  };

  using m_Cell = Cell<Float, DIM>;
  std::vector<m_Cell> m_cells;
  size_t m_nx;
  size_t m_ny;
  std::vector<size_t> m_cut_cell_idxs;

  Float m_x_min;
  Float m_x_max;
  Float m_y_min;
  Float m_y_max;
  Float m_min_delta;

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
  [[nodiscard]] constexpr auto to_corner_idx(size_t xi,
                                             size_t yi) const noexcept -> m_CornerIndices {
    assert(xi < m_nx);
    assert(yi < m_ny);
    return {
        .top_left     = (xi + 0) + (m_nx + 1) * (yi + 0),
        .top_right    = (xi + 1) + (m_nx + 1) * (yi + 0),
        .bottom_left  = (xi + 0) + (m_nx + 1) * (yi + 1),
        .bottom_right = (xi + 1) + (m_nx + 1) * (yi + 1),
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
        m_y_max(y_max) {
    assert(nx > 0);
    assert(ny > 0);

    const auto dx = (x_max - x_min) / static_cast<Float>(nx);
    const auto dy = (y_max - y_min) / static_cast<Float>(ny);
    m_min_delta   = std::min(dx, dy);

    for (size_t yi = 0; yi < ny; ++yi) {
      for (size_t xi = 0; xi < nx; ++xi) {
        // clang-format off
        m_cells[to_vec_idx(xi, yi)] = m_Cell{
            .x_min      = x_min + static_cast<Float>(xi) * dx,
            .dx         = dx,
            .y_min      = y_min + static_cast<Float>(yi) * dy,
            .dy         = dy,
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
          value(0) = f(cell.x_min + static_cast<Float>(0.5) * cell.dx,
                       cell.y_min + static_cast<Float>(0.5) * cell.dy);
        } else {
          value = f(cell.x_min + static_cast<Float>(0.5) * cell.dx,
                    cell.y_min + static_cast<Float>(0.5) * cell.dy);
        }
      } else if (cell.is_cut()) {
        Eigen::Vector<Float, 2> left_centroid = get_centroid(get_left_points<m_Cell, Float>(cell));
        Eigen::Vector<Float, 2> right_centroid =
            get_centroid(get_right_points<m_Cell, Float>(cell));

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
          value(0) = (f(cell.x_min + cell.dx / 4, cell.y_min + cell.dy / 4) +
                      f(cell.x_min + 3 * cell.dx / 4, cell.y_min + cell.dy / 4) +
                      f(cell.x_min + cell.dx / 4, cell.y_min + 3 * cell.dy / 4) +
                      f(cell.x_min + 3 * cell.dx / 4, cell.y_min + 3 * cell.dy / 4)) /
                     4;
        } else {
          value = (f(cell.x_min + cell.dx / 4, cell.y_min + cell.dy / 4) +
                   f(cell.x_min + 3 * cell.dx / 4, cell.y_min + cell.dy / 4) +
                   f(cell.x_min + cell.dx / 4, cell.y_min + 3 * cell.dy / 4) +
                   f(cell.x_min + 3 * cell.dx / 4, cell.y_min + 3 * cell.dy / 4)) /
                  4;
        }
      } else if (cell.is_cut()) {
        auto& cut_value = cell.get_cut();

        // Left side
        {
          const auto corners = get_left_points<m_Cell, Float>(cell);

          for (size_t j = 0; j < corners.size(); ++j) {
            Eigen::Vector<Float, 2> point = Eigen::Vector<Float, 2>::Zero();
            for (size_t i = 0; i < corners.size(); ++i) {
              point += (1 + (i == j)) * corners[i] / (corners.size() + 1);
            }

            if constexpr (DIM == 1) {
              cut_value.left_value(0) += f(point(X), point(Y));
            } else {
              cut_value.left_value += f(point(X), point(Y));
            }
          }
          cut_value.left_value /= static_cast<Float>(corners.size());
        }

        // Right side
        {
          const auto corners = get_right_points<m_Cell, Float>(cell);

          for (size_t j = 0; j < corners.size(); ++j) {
            Eigen::Vector<Float, 2> point = Eigen::Vector<Float, 2>::Zero();
            for (size_t i = 0; i < corners.size(); ++i) {
              point += (1 + (i == j)) * corners[i] / (corners.size() + 1);
            }

            if constexpr (DIM == 1) {
              cut_value.right_value(0) += f(point(X), point(Y));
            } else {
              cut_value.right_value += f(point(X), point(Y));
            }
          }
          cut_value.right_value /= static_cast<Float>(corners.size());
        }
      } else {
        Igor::Panic("Unknown cell type with variant index {}.", cell.value.index());
      }
    }
  }

 private:
  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto approx_eq(Float a, Float b) const noexcept -> bool {
    return std::abs(a - b) <= EPS<Float>;
  };

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto classify_cut(const m_Cell& cell,
                                            Point<Float>& cut1_point,
                                            Point<Float>& cut2_point) const noexcept -> CutType {
    const int cut1_on_x = approx_eq(cell.x_min, cut1_point(X)) * ON_MIN +
                          approx_eq(cell.x_min + cell.dx, cut1_point(X)) * ON_MAX;
    const int cut1_on_y = approx_eq(cell.y_min, cut1_point(Y)) * ON_MIN +
                          approx_eq(cell.y_min + cell.dy, cut1_point(Y)) * ON_MAX;
    assert(cut1_on_x != NOT_ON || cut1_on_y != NOT_ON);
    Side cut1_loc = cut1_on_x == ON_MIN   ? LEFT
                    : cut1_on_y == ON_MIN ? BOTTOM
                    : cut1_on_x == ON_MAX ? RIGHT
                                          : TOP;

    const int cut2_on_x = approx_eq(cell.x_min, cut2_point(X)) * ON_MIN +
                          approx_eq(cell.x_min + cell.dx, cut2_point(X)) * ON_MAX;
    const int cut2_on_y = approx_eq(cell.y_min, cut2_point(Y)) * ON_MIN +
                          approx_eq(cell.y_min + cell.dy, cut2_point(Y)) * ON_MAX;
    assert(cut2_on_x != NOT_ON || cut2_on_y != NOT_ON);
    Side cut2_loc = cut2_on_x == ON_MIN   ? LEFT
                    : cut2_on_y == ON_MIN ? BOTTOM
                    : cut2_on_x == ON_MAX ? RIGHT
                                          : TOP;

    if (cut1_loc == cut2_loc) {
      std::stringstream s{};
      s << cell;
      Igor::Debug("Cell to cut: {}", s.str());
      Igor::Debug("cut1_loc = {}", static_cast<std::underlying_type_t<Side>>(cut1_loc));
      Igor::Debug("cut2_loc = {}", static_cast<std::underlying_type_t<Side>>(cut2_loc));
      IGOR_DEBUG_PRINT(cut1_point);
      IGOR_DEBUG_PRINT(cut2_point);
      Igor::Todo("Cut only on one side. I think we can just ignore then.");
    }
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

    return type;
  }

  [[nodiscard]] constexpr auto
  find_next_cell_to_cut(const m_Cell& cell,
                        const Point<Float>& exit_point) const noexcept -> size_t {
    assert(point_in_cell(exit_point, cell));

    const int on_x = approx_eq(cell.x_min, exit_point(X)) * ON_MIN +
                     approx_eq(cell.x_min + cell.dx, exit_point(X)) * ON_MAX;
    const int on_y = approx_eq(cell.y_min, exit_point(Y)) * ON_MIN +
                     approx_eq(cell.y_min + cell.dy, exit_point(Y)) * ON_MAX;
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
        Igor::Panic("Expected to exit on top, not on left.");
        return cell.left_idx;
      }
      // Right side
      else if (on_x == ON_MAX) {
        assert(is_cell(cell.right_idx));
        Igor::Panic("Expected to exit on top, not on right.");
        return cell.right_idx;
      }
      // Bottom side
      else if (on_y == ON_MIN) {
        assert(is_cell(cell.bottom_idx));
        Igor::Panic("Expected to exit on top, not on bottom.");
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
 public:
  template <typename ParamCurve>
  [[nodiscard]] constexpr auto cut_curve(ParamCurve curve) noexcept -> bool {
    enum { T_BELOW_EXIT = -1, T_IS_EXIT, T_ABOVE_EXIT };

    const auto t_location = [curve](Float t, const m_Cell& cell, Float t_entry) -> int {
      const auto pos = curve(t);
      if (t <= t_entry || (point_in_cell(pos, cell) && !point_on_boundary(pos, cell))) {
        return T_BELOW_EXIT;
      } else if (point_on_boundary(pos, cell)) {
        return T_IS_EXIT;
      }
      return T_ABOVE_EXIT;
    };

    Float t             = 0;
    const auto init_pos = curve(t);
    const auto cell_it =
        std::find_if(std::cbegin(m_cells), std::cend(m_cells), [&](const auto& cell) {
          return point_in_cell(init_pos, cell);
        });
    if (cell_it == std::cend(m_cells)) {
      Igor::Warn("Point {} on parametrized curve is not in the grid.", init_pos);
      return false;
    }
    assert(std::distance(std::cbegin(m_cells), cell_it) >= 0);
    auto cell_idx = static_cast<size_t>(std::distance(std::cbegin(m_cells), cell_it));

    if (!m_cut_cell_idxs.empty()) {
      Igor::Warn("Grid is already cut, this method expects the grid to no be cut.");
      return false;
    }

    m_cut_cell_idxs.push_back(static_cast<size_t>(cell_idx));
    std::vector<Point<Float>> entry_points{};
    entry_points.push_back(init_pos);

    while (t < 1) {
      const auto& cell = m_cells[cell_idx];
      Float t_lower    = t;
      Float t_upper    = std::min(Float{1}, t + 0.5);

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
          IGOR_DEBUG_PRINT(point_on_boundary(curve(t_upper), cell));
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
      const auto next_pos = curve(t);
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
          .x1_cut      = cut1_point(X),
          .y1_cut      = cut1_point(Y),
          .x2_cut      = cut2_point(X),
          .y2_cut      = cut2_point(Y),
      };
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto point_in_grid(const Point<Float>& point) const noexcept -> bool {
    return point(X) - m_x_min >= -EPS<Float> && m_x_max - point(X) >= -EPS<Float> &&  //
           point(Y) - m_y_min >= -EPS<Float> && m_y_max - point(Y) >= -EPS<Float>;
  }

  [[nodiscard]] constexpr auto
  cut_piecewise_linear(std::vector<Point<Float>> points) noexcept -> bool {
    // Remove points outside of grid
    std::erase_if(points, [this](const Point<Float>& p) { return !point_in_grid(p); });
    assert(points.size() >= 2);

    const auto cell_it =
        std::find_if(std::cbegin(m_cells), std::cend(m_cells), [&](const auto& cell) {
          return point_in_cell(points[0], cell);
        });
    assert(cell_it != std::cend(m_cells));

    assert(std::distance(std::cbegin(m_cells), cell_it) >= 0);
    auto cell_idx = static_cast<size_t>(std::distance(std::cbegin(m_cells), cell_it));

    if (!m_cut_cell_idxs.empty()) {
      Igor::Warn("Grid is already cut, this method expects the grid to no be cut.");
      return false;
    }
    std::vector<Point<Float>> entry_points{};

    // Handle first point; extend curve backwards to cut cell containing the first point
    {
      const auto& cell = m_cells[cell_idx];
      const auto& p0   = points[0];
      const auto& p1   = points[1];

      const Point<Float> s = p1 - p0;
      const std::array<Float, 4> rs{
          (cell.x_min - p0(X)) / s(X),
          (cell.x_min + cell.dx - p0(X)) / s(X),
          (cell.y_min - p0(Y)) / s(Y),
          (cell.y_min + cell.dy - p0(Y)) / s(Y),
      };

      Float r_entry = -std::numeric_limits<Float>::max();
      for (Float r : rs) {
        if (r < EPS<Float> && r > r_entry) { r_entry = r; }
      }
      if (!(r_entry != -std::numeric_limits<Float>::max())) {
        IGOR_DEBUG_PRINT(s);
        IGOR_DEBUG_PRINT(r_entry);
        IGOR_DEBUG_PRINT(rs);
        IGOR_DEBUG_PRINT(p0);
        IGOR_DEBUG_PRINT(p1);
      }
      assert(r_entry != -std::numeric_limits<Float>::max());

      entry_points.emplace_back(p0 + r_entry * s);
      m_cut_cell_idxs.push_back(cell_idx);
    }

    // Handle all point pairs
    for (size_t point_idx = 0; point_idx < points.size() - 1; ++point_idx) {
      auto p0        = points[point_idx];
      const auto& p1 = points[point_idx + 1];

      const Point<Float> s = (p1 - p0).normalized();

      while (!point_in_cell(p1, m_cells[cell_idx])) {
        const auto& cell = m_cells[cell_idx];
        if (!point_in_cell(p0, cell)) {
          std::cout << "cell = " << cell << '\n';
          IGOR_DEBUG_PRINT(p0);
          IGOR_DEBUG_PRINT(p1);
        }
        assert(point_in_cell(p0, cell));

        const std::array<Float, 4> rs{
            (cell.x_min - p0(X)) / s(X),
            (cell.x_min + cell.dx - p0(X)) / s(X),
            (cell.y_min - p0(Y)) / s(Y),
            (cell.y_min + cell.dy - p0(Y)) / s(Y),
        };
        Float r_exit = std::numeric_limits<Float>::max();
        for (Float r : rs) {
          if (r > -EPS<Float> && r < r_exit) { r_exit = r; }
        }
        if (!(r_exit < std::numeric_limits<Float>::max())) {
          IGOR_DEBUG_PRINT(s);
          IGOR_DEBUG_PRINT(r_exit);
          IGOR_DEBUG_PRINT(rs);
          IGOR_DEBUG_PRINT(p0);
          IGOR_DEBUG_PRINT(p1);
        }
        assert(r_exit < std::numeric_limits<Float>::max());

        const Point<Float> p0_next = p0 + r_exit * s;
        if (!point_in_cell(p0_next, cell)) {
          std::cout << "cell = " << cell << '\n';
          IGOR_DEBUG_PRINT(p0);
          IGOR_DEBUG_PRINT(p0_next);
          IGOR_DEBUG_PRINT(p1);
          IGOR_DEBUG_PRINT(s);
          IGOR_DEBUG_PRINT(r_exit);
          IGOR_DEBUG_PRINT(entry_points);
          IGOR_DEBUG_PRINT(m_cut_cell_idxs);
          IGOR_DEBUG_PRINT(point_idx);
          IGOR_DEBUG_PRINT(points.size());
        }

        p0 = p0_next;
        entry_points.push_back(p0);
        cell_idx = find_next_cell_to_cut(cell, entry_points.back());
        m_cut_cell_idxs.push_back(cell_idx);
      }
    }

    // Handle last point; extend curve backwards to cut cell containing the first point
    {
      const auto& cell = m_cells[cell_idx];
      const auto& p0   = points[points.size() - 2];
      const auto& p1   = points[points.size() - 1];

      const Point<Float> s = p1 - p0;
      const std::array<Float, 4> rs{
          (cell.x_min - p1(X)) / s(X),
          (cell.x_min + cell.dx - p1(X)) / s(X),
          (cell.y_min - p1(Y)) / s(Y),
          (cell.y_min + cell.dy - p1(Y)) / s(Y),
      };

      Float r_exit = std::numeric_limits<Float>::max();
      for (Float r : rs) {
        if (r > -EPS<Float> && r < r_exit) { r_exit = r; }
      }
      assert(r_exit != std::numeric_limits<Float>::max());

      entry_points.emplace_back(p1 + r_exit * s);
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
          .x1_cut      = cut1_point(X),
          .y1_cut      = cut1_point(Y),
          .x2_cut      = cut2_point(X),
          .y2_cut      = cut2_point(Y),
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
      const auto left_polygon = cell.get_cut_left_polygon();

      const auto right_value   = cell.get_cut().right_value;
      const auto right_polygon = cell.get_cut_right_polygon();

      const auto cell_polygon = cell.get_cartesian_polygon();

      cell.value = CartesianValue<Float, DIM>{
          .value = (left_value * left_polygon.area() + right_value * right_polygon.area()) /
                   cell_polygon.area()};
    }
    m_cut_cell_idxs.clear();
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto abs_max_value() const noexcept -> Float {
    if constexpr (DIM == 1) {
      assert(m_cells.size() > 0);
      return std::transform_reduce(
          std::cbegin(m_cells),
          std::cend(m_cells),
          Float{0},
          [](const auto& a, const auto& b) { return std::max(a, b); },
          [](const m_Cell& cell) {
            if (cell.is_cartesian()) {
              return std::abs(cell.get_cartesian().value(0));
            } else if (cell.is_cut()) {
              const auto& value = cell.get_cut();
              return std::max(std::abs(value.left_value(0)), std::abs(value.right_value(0)));
            } else {
              Igor::Panic("Unknown variant entry with index {}", cell.value.index());
              std::unreachable();
            }
          });
    } else {
      assert(m_cells.size() > 0);
      return std::transform_reduce(
          std::cbegin(m_cells),
          std::cend(m_cells),
          Float{0},
          [](const auto& a, const auto& b) { return std::max(a, b); },
          [](const m_Cell& cell) {
            if (cell.is_cartesian()) {
              return std::transform_reduce(
                  std::cbegin(cell.get_cartesian().value),
                  std::cend(cell.get_cartesian().value),
                  Float{0},
                  [](const auto& a, const auto& b) { return std::max(a, b); },
                  [](const auto& a) { return std::abs(a); });
            } else if (cell.is_cut()) {
              return std::max(std::transform_reduce(
                                  std::cbegin(cell.get_cut().left_value),
                                  std::cend(cell.get_cut().left_value),
                                  Float{0},
                                  [](const auto& a, const auto& b) { return std::max(a, b); },
                                  [](const auto& a) { return std::abs(a); }),
                              std::transform_reduce(
                                  std::cbegin(cell.get_cut().right_value),
                                  std::cend(cell.get_cut().right_value),
                                  Float{0},
                                  [](const auto& a, const auto& b) { return std::max(a, b); },
                                  [](const auto& a) { return std::abs(a); }));
            } else {
              Igor::Panic("Unknown variant entry with index {}", cell.value.index());
              std::unreachable();
            }
          });
    }
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto get_shock_curve() const noexcept -> std::vector<Point<Float>> {
    if (m_cut_cell_idxs.empty()) { return {}; }

    std::vector<Point<Float>> curve_points(m_cut_cell_idxs.size() + 1);

    {
      const auto& cell = m_cells[m_cut_cell_idxs.front()];
      curve_points[0]  = {cell.get_cut().x1_cut, cell.get_cut().y1_cut};
    }

    for (size_t i = 0; i < m_cut_cell_idxs.size(); ++i) {
      const auto& cell    = m_cells[m_cut_cell_idxs[i]];
      curve_points[i + 1] = {cell.get_cut().x2_cut, cell.get_cut().y2_cut};
    }

    return curve_points;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto nx() const noexcept -> size_t { return m_nx; }
  [[nodiscard]] constexpr auto ny() const noexcept -> size_t { return m_ny; }
  [[nodiscard]] constexpr auto min_delta() const noexcept -> Float { return m_min_delta; }
  [[nodiscard]] constexpr auto x_min() const noexcept -> Float { return m_x_min; }
  [[nodiscard]] constexpr auto x_max() const noexcept -> Float { return m_x_max; }
  [[nodiscard]] constexpr auto y_min() const noexcept -> Float { return m_y_min; }
  [[nodiscard]] constexpr auto y_max() const noexcept -> Float { return m_y_max; }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto size() const noexcept -> size_t { return m_cells.size(); }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto operator[](size_t idx) const noexcept -> const m_Cell& {
    assert(idx < m_nx * m_ny);
    return m_cells[idx];
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto operator[](size_t idx) noexcept -> m_Cell& {
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
  [[nodiscard]] constexpr auto is_left(size_t idx) const noexcept -> bool {
    const auto [xi, yi] = to_mat_idx(idx);
    return xi == 0;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto is_bottom(size_t idx) const noexcept -> bool {
    const auto [xi, yi] = to_mat_idx(idx);
    return yi == 0;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto
  get_existing_neighbours(size_t idx) const noexcept -> SmallVector<size_t> {
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
  template <typename, typename>
  friend class Solver;

  template <Zap::IO::VTKFormat F>
  friend class Zap::IO::VTKWriter;

  template <typename F, size_t D>
  friend class Zap::IO::IncCellWriter;

  template <typename F, size_t D>
  friend class Zap::IO::IncCellReader;
};

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_GRID_HPP_
