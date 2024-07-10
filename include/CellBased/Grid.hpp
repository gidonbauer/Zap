#ifndef ZAP_CELL_BASED_GRID_HPP_
#define ZAP_CELL_BASED_GRID_HPP_

#include <numeric>
#include <optional>

#include <AD/ad.hpp>

#include "CellBased/Cell.hpp"
#include "IO/Fwd.hpp"

namespace Zap::CellBased {

enum Side : int {
  BOTTOM = 0b0001,
  RIGHT  = 0b0010,
  TOP    = 0b0100,
  LEFT   = 0b1000,
  ALL    = LEFT | RIGHT | BOTTOM | TOP,
};

template <typename Float, size_t DIM>
class Grid {
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

  mutable std::optional<Float> m_min_delta = std::nullopt;
  mutable std::optional<Float> m_x_min     = std::nullopt;
  mutable std::optional<Float> m_x_max     = std::nullopt;
  mutable std::optional<Float> m_y_min     = std::nullopt;
  mutable std::optional<Float> m_y_max     = std::nullopt;

  // -----------------------------------------------------------------------------------------------
  constexpr Grid(size_t nx, size_t ny)
      : m_cells(nx * ny),
        m_nx(nx),
        m_ny(ny) {}

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto to_vec_idx(size_t xi, size_t yi) const noexcept -> size_t {
    assert(xi < m_nx);
    assert(yi < m_ny);
    return yi * m_nx + xi;
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
  [[nodiscard]] static constexpr auto
  Uniform(Float x_min, Float x_max, size_t nx, Float y_min, Float y_max, size_t ny) -> Grid {
    assert(nx > 0);
    assert(ny > 0);
    Grid grid(nx, ny);
    const auto dx    = (x_max - x_min) / static_cast<Float>(nx);
    const auto dy    = (y_max - y_min) / static_cast<Float>(ny);
    grid.m_min_delta = std::min(dx, dy);
    grid.m_x_min     = x_min;
    grid.m_x_max     = x_max;
    grid.m_y_min     = y_min;
    grid.m_y_max     = y_max;

    for (size_t yi = 0; yi < ny; ++yi) {
      for (size_t xi = 0; xi < nx; ++xi) {
        // clang-format off
        grid.m_cells[grid.to_vec_idx(xi, yi)] = m_Cell{
            .x_min      = x_min + static_cast<Float>(xi) * dx,
            .dx         = dx,
            .y_min      = y_min + static_cast<Float>(yi) * dy,
            .dy         = dy,
            .left_idx   = xi == 0        ? NULL_INDEX : grid.to_vec_idx(xi - 1, yi),
            .right_idx  = xi == (nx - 1) ? NULL_INDEX : grid.to_vec_idx(xi + 1, yi),
            .bottom_idx = yi == 0        ? NULL_INDEX : grid.to_vec_idx(xi, yi - 1),
            .top_idx    = yi == (ny - 1) ? NULL_INDEX : grid.to_vec_idx(xi, yi + 1),
        };
        // clang-format on
      }
    }
    return grid;
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
    constexpr Float eps = 1e-8;
    auto approx_eq      = [](Float a, Float b) { return std::abs(a - b) <= eps; };
    enum { ON_MIN = -1, NOT_ON = 0, ON_MAX = 1 };

    size_t counter = 0;
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
        Igor::Debug("Cell index = {}", counter);
        auto& cell_value = cell.get_cut();
        {
          std::stringstream s{};
          s << cell;
          Igor::Debug("cell = {}", s.str());
        }

        const int cut1_on_x = approx_eq(cell.x_min, cell_value.x1_cut) * ON_MIN +
                              approx_eq(cell.x_min + cell.dx, cell_value.x1_cut) * ON_MAX;
        const int cut1_on_y = approx_eq(cell.y_min, cell_value.y1_cut) * ON_MIN +
                              approx_eq(cell.y_min + cell.dy, cell_value.y1_cut) * ON_MAX;

        const int cut2_on_x = approx_eq(cell.x_min, cell_value.x2_cut) * ON_MIN +
                              approx_eq(cell.x_min + cell.dx, cell_value.x2_cut) * ON_MAX;
        const int cut2_on_y = approx_eq(cell.y_min, cell_value.y2_cut) * ON_MIN +
                              approx_eq(cell.y_min + cell.dy, cell_value.y2_cut) * ON_MAX;
        IGOR_DEBUG_PRINT(cut1_on_x);
        IGOR_DEBUG_PRINT(cut1_on_y);
        IGOR_DEBUG_PRINT(cut2_on_x);
        IGOR_DEBUG_PRINT(cut2_on_y);

        assert(cut1_on_x != NOT_ON || cut1_on_y != NOT_ON);
        assert(cut2_on_x != NOT_ON || cut2_on_y != NOT_ON);

        Igor::Todo("Implement fill_center for cut cells.");
      } else {
        Igor::Panic("Unknown cell type with variant index {}.", cell.value.index());
      }
      ++counter;
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <typename FUNC>
  constexpr void fill_four_point(FUNC f) noexcept {
    for (auto& cell : m_cells) {
      if (cell.is_cartesian()) {
        auto& value = cell.get_cartesian();
        if constexpr (DIM == 1) {
          value(0) = (f(cell.x_min + cell.dx / 3, cell.y_min + cell.dy / 3) +
                      f(cell.x_min + 2 * cell.dx / 3, cell.y_min + cell.dy / 3) +
                      f(cell.x_min + cell.dx / 3, cell.y_min + 2 * cell.dy / 3) +
                      f(cell.x_min + 2 * cell.dx / 3, cell.y_min + 2 * cell.dy / 3)) /
                     4;
        } else {
          value = (f(cell.x_min + cell.dx / 3, cell.y_min + cell.dy / 3) +
                   f(cell.x_min + 2 * cell.dx / 3, cell.y_min + cell.dy / 3) +
                   f(cell.x_min + cell.dx / 3, cell.y_min + 2 * cell.dy / 3) +
                   f(cell.x_min + 2 * cell.dx / 3, cell.y_min + 2 * cell.dy / 3)) /
                  4;
        }
      } else {
        Igor::Todo("Only cartesian cells are implemented.");
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <typename ParamShock>
  [[nodiscard]] constexpr auto cut_init_shock(ParamShock init_shock) noexcept -> bool {
    using ad_t          = ad::gt1s<Float>::type;
    constexpr Float eps = 1e-8;
    enum { ON_MIN = -1, NOT_ON = 0, ON_MAX = 1 };
    enum { X, Y, POS_SIZE };

    auto approx_eq       = [](Float a, Float b) -> bool { return std::abs(a - b) <= eps; };
    auto pos_on_boundary = [approx_eq](const Eigen::Vector<Float, POS_SIZE>& pos,
                                       const m_Cell& cell) -> bool {
      return ((approx_eq(pos(X), cell.x_min) || approx_eq(pos(X), cell.x_min + cell.dx)) &&
              (pos(Y) >= cell.y_min && pos(Y) <= cell.y_min + cell.dy)) ||
             ((approx_eq(pos(Y), cell.y_min) || approx_eq(pos(Y), cell.y_min + cell.dy)) &&
              (pos(X) >= cell.x_min && pos(X) <= cell.x_min + cell.dx));
    };
    auto pos_in_cell = [pos_on_boundary](const Eigen::Vector<Float, POS_SIZE>& pos,
                                         const m_Cell& cell) -> bool {
      return ((pos(X) >= cell.x_min) && (pos(X) <= cell.x_min + cell.dx) &&  //
              (pos(Y) >= cell.y_min) && (pos(Y) <= cell.y_min + cell.dy)) ||
             pos_on_boundary(pos, cell);
      // (pos(X) - cell.x_min) <= cell.dx && (pos(1) - cell.y_min) <= cell.dy;
    };
    auto finished_iteration = [approx_eq](Float t) {
      return t > Float{1} || approx_eq(t, Float{1});
    };

    const auto init_pos = init_shock(Float{0});
    const auto cell_it =
        std::find_if(std::cbegin(m_cells), std::cend(m_cells), [&](const auto& cell) {
          return pos_in_cell(init_pos, cell);
        });
    if (cell_it == std::cend(m_cells)) {
      Igor::Warn("Point {} on parametrized initial shock is not in the grid.", init_pos);
      return false;
    }
    assert(std::distance(std::cbegin(m_cells), cell_it) >= 0);
    auto cell_idx = static_cast<size_t>(std::distance(std::cbegin(m_cells), cell_it));

    ad_t t_ad            = 0.0;
    ad::derivative(t_ad) = 1.0;

    std::vector<size_t> cut_cells{};
    cut_cells.push_back(static_cast<size_t>(cell_idx));
    std::vector<Eigen::Vector<Float, POS_SIZE>> entry_points{};
    entry_points.push_back(init_pos);

    for ([[maybe_unused]] size_t iter = 0; true; ++iter) {
      const auto& cell = m_cells.at(cell_idx);

      auto pos_ad = init_shock(t_ad);

      const Eigen::Vector<Float, POS_SIZE> pos{ad::value(pos_ad(X)), ad::value(pos_ad(Y))};
      const Eigen::Vector<Float, POS_SIZE> der{ad::derivative(pos_ad(X)),
                                               ad::derivative(pos_ad(Y))};

      const auto in_cell     = pos_in_cell(pos, cell);
      const auto r_condition = [approx_eq, in_cell](Float r) {
        return !approx_eq(r, Float{0}) &&       // Not 0
               ((in_cell && (r > Float{0})) ||  // If in cell, then r > 0
                (!in_cell && (r < Float{0})));  // If not in cell, then r < 0
      };
      const auto update_r = [r_condition](Float& r, Float prop_r) {
        if (r_condition(prop_r) && std::abs(r) > std::abs(prop_r)) {
          r = prop_r;
        }
      };
      Float r = std::numeric_limits<Float>::infinity();
      if (der(X) != Float{0}) {
        update_r(r, (cell.x_min - pos(X)) / der(X));
        update_r(r, (cell.x_min + cell.dx - pos(X)) / der(X));
      }
      if (der(Y) != Float{0}) {
        update_r(r, (cell.y_min - pos(Y)) / der(Y));
        update_r(r, (cell.y_min + cell.dy - pos(Y)) / der(Y));
      }
      if (std::isinf(r)) {
        std::stringstream s{};
        s << cell;
        Igor::Debug("t = {}", ad::value(t_ad));
        Igor::Debug("pos = {}", pos);
        Igor::Debug("der = {}", der);
        Igor::Warn("Did not find exit point for cell {}.", s.str());
        return false;
      }

      ad::value(t_ad) += r;

      const auto next_pos_ad = init_shock(t_ad);
      Eigen::Vector<Float, POS_SIZE> next_pos{ad::value(next_pos_ad(X)), ad::value(next_pos_ad(Y))};
      const Eigen::Vector<Float, POS_SIZE> next_der{ad::derivative(next_pos_ad(X)),
                                                    ad::derivative(next_pos_ad(Y))};
      if (pos_on_boundary(next_pos, cell)) {
        entry_points.push_back(next_pos);

        if (finished_iteration(ad::value(t_ad))) {
          break;
        }

        const int on_x = approx_eq(cell.x_min, next_pos(X)) * ON_MIN +
                         approx_eq(cell.x_min + cell.dx, next_pos(X)) * ON_MAX;
        const int on_y = approx_eq(cell.y_min, next_pos(Y)) * ON_MIN +
                         approx_eq(cell.y_min + cell.dy, next_pos(Y)) * ON_MAX;
        assert(on_x != NOT_ON || on_y != NOT_ON);

        if (on_x != NOT_ON && on_y != NOT_ON) {
          Igor::Todo("Handle exit on corner, e.g. using `next_der`.");
        }

        if (on_x == ON_MIN) {
          assert(is_cell(cell.left_idx));
          cut_cells.push_back(cell.left_idx);
          cell_idx = cell.left_idx;
        } else if (on_x == ON_MAX) {
          assert(is_cell(cell.right_idx));
          cut_cells.push_back(cell.right_idx);
          cell_idx = cell.right_idx;
        } else if (on_y == ON_MIN) {
          assert(is_cell(cell.bottom_idx));
          cut_cells.push_back(cell.bottom_idx);
          cell_idx = cell.bottom_idx;
        } else if (on_y == ON_MAX) {
          assert(is_cell(cell.top_idx));
          cut_cells.push_back(cell.top_idx);
          cell_idx = cell.top_idx;
        }
      }
    }

    assert(cut_cells.size() + 1 == entry_points.size());
    for (size_t i = 0; i < cut_cells.size(); ++i) {
      auto& cell_to_cut = m_cells[cut_cells[i]];
      auto cut1_point   = entry_points[i];
      auto cut2_point   = entry_points[i + 1];

      const int cut1_on_x = approx_eq(cell_to_cut.x_min, cut1_point(X)) * ON_MIN +
                            approx_eq(cell_to_cut.x_min + cell_to_cut.dx, cut1_point(X)) * ON_MAX;
      const int cut1_on_y = approx_eq(cell_to_cut.y_min, cut1_point(Y)) * ON_MIN +
                            approx_eq(cell_to_cut.y_min + cell_to_cut.dy, cut1_point(Y)) * ON_MAX;
      assert(cut1_on_x != NOT_ON || cut1_on_y != NOT_ON);
      Side cut1_loc = cut1_on_x == ON_MIN   ? LEFT
                      : cut1_on_y == ON_MIN ? BOTTOM
                      : cut1_on_x == ON_MAX ? RIGHT
                                            : TOP;

      const int cut2_on_x = approx_eq(cell_to_cut.x_min, cut2_point(X)) * ON_MIN +
                            approx_eq(cell_to_cut.x_min + cell_to_cut.dx, cut2_point(X)) * ON_MAX;
      const int cut2_on_y = approx_eq(cell_to_cut.y_min, cut2_point(Y)) * ON_MIN +
                            approx_eq(cell_to_cut.y_min + cell_to_cut.dy, cut2_point(Y)) * ON_MAX;
      assert(cut2_on_x != NOT_ON || cut2_on_y != NOT_ON);
      Side cut2_loc = cut2_on_x == ON_MIN   ? LEFT
                      : cut2_on_y == ON_MIN ? BOTTOM
                      : cut2_on_x == ON_MAX ? RIGHT
                                            : TOP;

      if (cut1_loc == cut2_loc) {
        Igor::Todo("Cut only on one side. I think we can just ignore then.");
      }
      if (cut1_loc > cut2_loc) {
        std::swap(cut1_loc, cut2_loc);
        std::swap(cut1_point, cut2_point);
      }

      CutType type;
      switch (cut1_loc | cut2_loc) {
        case LEFT | BOTTOM:
          type = CutType::BOTTOM_LEFT;
          break;
        case BOTTOM | RIGHT:
          type = CutType::BOTTOM_RIGHT;
          break;
        case RIGHT | TOP:
          type = CutType::TOP_RIGHT;
          break;
        case TOP | LEFT:
          type = CutType::TOP_LEFT;
          break;
        case LEFT | RIGHT:
          type = CutType::MIDDLE_HORI;
          break;
        case BOTTOM | TOP:
          type = CutType::MIDDLE_VERT;
          break;
        default:
          Igor::Panic("Invalid combination cut1_loc = {} and cut2_loc = {}",
                      static_cast<std::underlying_type_t<Side>>(cut1_loc),
                      static_cast<std::underlying_type_t<Side>>(cut2_loc));
          std::unreachable();
      }

      cell_to_cut.value = CutValue<Float, DIM>{
          .left_value  = Eigen::Vector<Float, DIM>::Zero(),
          .right_value = Eigen::Vector<Float, DIM>::Zero(),
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
      Igor::Todo("Not implemented for DIM={}", DIM);
      return 0;
    }
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto nx() const noexcept -> size_t { return m_nx; }
  [[nodiscard]] constexpr auto ny() const noexcept -> size_t { return m_ny; }

  [[nodiscard]] constexpr auto min_delta() const noexcept -> Float {
    if (!m_min_delta.has_value()) {
      assert(m_cells.size() > 0);
      m_min_delta = std::transform_reduce(
          std::cbegin(m_cells),
          std::cend(m_cells),
          std::min(m_cells[0].dx, m_cells[0].dy),
          [](const auto& a, const auto& b) { return std::min(a, b); },
          [](const m_Cell& cell) { return std::min(cell.dx, cell.dy); });
    }
    return *m_min_delta;
  }

  [[nodiscard]] constexpr auto x_min() const noexcept -> Float {
    if (!m_x_min.has_value()) {
      Igor::Todo("Implement finding x_min.");
      m_x_min = 0;
    }
    return *m_x_min;
  }

  [[nodiscard]] constexpr auto x_max() const noexcept -> Float {
    if (!m_x_max.has_value()) {
      Igor::Todo("Implement finding x_max.");
      m_x_max = 0;
    }
    return *m_x_max;
  }

  [[nodiscard]] constexpr auto y_min() const noexcept -> Float {
    if (!m_y_min.has_value()) {
      Igor::Todo("Implement finding y_min.");
      m_y_min = 0;
    }
    return *m_y_min;
  }

  [[nodiscard]] constexpr auto y_max() const noexcept -> Float {
    if (!m_y_max.has_value()) {
      Igor::Todo("Implement finding y_max.");
      m_y_max = 0;
    }
    return *m_y_max;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto operator[](size_t idx) const noexcept -> const m_Cell& {
    assert(idx < m_nx * m_ny);
    return m_cells[idx];
  }

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
  template <typename F, typename G>
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
