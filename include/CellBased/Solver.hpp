#ifndef ZAP_CELL_BASED_SOLVER_HPP_
#define ZAP_CELL_BASED_SOLVER_HPP_

#include "CellBased/Grid.hpp"

namespace Zap::CellBased {

template <typename FluxX, typename FluxY>
class Solver {
  FluxX m_numerical_flux_x;
  FluxY m_numerical_flux_y;

  template <typename Float, size_t DIM>
  requires(DIM == 1)
  [[nodiscard]] constexpr auto
  apply_boundary_flux(const Cell<Float, DIM>& cell, Side side, size_t idx) const noexcept -> Float {
    switch (idx) {
      case SAME_VALUE_INDEX:
        {
          assert(cell.is_cartesian() && "Non-cartesian cell not implemented.");
          const auto& cell_value = cell.get_cartesian().value;
          switch (side) {
            case LEFT:
            case RIGHT:  return m_numerical_flux_x(cell_value(0), cell_value(0));

            case BOTTOM:
            case TOP:    return m_numerical_flux_y(cell_value(0), cell_value(0));

            default:     Igor::Panic("Must be BOTTOM, RIGHT, TOP, or LEFT."); std::unreachable();
          }
        }
      case ZERO_FLUX_INDEX:
        {
          return static_cast<Float>(0);
        }
      case NULL_INDEX:
        {
          std::stringstream s{};
          s << cell;
          Igor::Panic("NULL index for cell {}, need to setup boundary condtions.", s.str());
        }
      default:
        {
          std::stringstream s{};
          s << cell;
          Igor::Panic("Unknown index type with value {} for cell {}.", idx, s.str());
        }
    }
    std::unreachable();
  }

  // -----------------------------------------------------------------------------------------------
  template <typename Float, size_t DIM>
  requires(DIM == 1)
  [[nodiscard]] constexpr auto
  apply_flux(const Grid<Float, DIM>& grid, const auto& cell, Side side) const noexcept -> Float {
    size_t idx = NULL_INDEX;
    switch (side) {
      case BOTTOM: idx = cell.bottom_idx; break;
      case RIGHT:  idx = cell.right_idx; break;
      case TOP:    idx = cell.top_idx; break;
      case LEFT:   idx = cell.left_idx; break;
      default:     Igor::Panic("Must be BOTTOM, RIGHT, TOP, or LEFT."); std::unreachable();
    }

    if (!grid.is_cell(idx)) {
      return apply_boundary_flux(cell, side, idx);
    }

    if (cell.is_cartesian()) {
      const auto& cell_value = cell.get_cartesian().value;

      if (grid.m_cells[idx].is_cartesian()) {
        const auto& other_cell_value = grid.m_cells[idx].get_cartesian().value;

        switch (side) {
          case LEFT:   return m_numerical_flux_x(other_cell_value(0), cell_value(0));
          case RIGHT:  return m_numerical_flux_x(cell_value(0), other_cell_value(0));

          case BOTTOM: return m_numerical_flux_y(other_cell_value(0), cell_value(0));
          case TOP:    return m_numerical_flux_y(cell_value(0), other_cell_value(0));

          default:     Igor::Panic("Must be BOTTOM, RIGHT, TOP, or LEFT."); std::unreachable();
        }
      } else if (grid.m_cells[idx].is_cut()) {
        Igor::Todo("Flux for cut cells is not implemented yet.");
        std::unreachable();
      } else {
        Igor::Panic("Unknown cell type with variant index {}", grid.m_cells[idx].value.index());
        std::unreachable();
      }
    } else if (cell.is_cut()) {
      Igor::Todo("Flux for cut cells is not implemented yet.");
      std::unreachable();
    } else {
      Igor::Panic("Unknown cell type with variant index {}", cell.value.index());
      std::unreachable();
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <typename Float, size_t DIM>
  requires(DIM > 1)
  [[nodiscard]] constexpr auto apply_flux(const Grid<Float, DIM>& grid,
                                          const auto& cell,
                                          Side side) const noexcept -> Eigen::Vector<Float, DIM> {
    (void)grid;
    (void)cell;
    (void)side;
    Igor::Todo("Not implemented for DIM={}", DIM);
    std::unreachable();
  }

 public:
  // -----------------------------------------------------------------------------------------------
  constexpr Solver(FluxX F, FluxY G)
      : m_numerical_flux_x(std::move(F)),
        m_numerical_flux_y(std::move(G)) {}

  // -----------------------------------------------------------------------------------------------
  template <typename Float, std::size_t DIM, typename GridWriter, typename TimeWriter>
  auto solve(Grid<Float, DIM> grid,
             Float tend,
             GridWriter& grid_writer,
             TimeWriter& time_writer) noexcept -> std::optional<Grid<Float, DIM>> {
    auto& grid_curr = grid;
    auto grid_next  = grid;

    if (!grid_writer.write_data(grid_curr) || !time_writer.write_data(Float{0})) {
      return std::nullopt;
    }

    for (Float t = 0.0; t < tend;) {
      const Float CFL_factor = grid_curr.abs_max_value();
      if (std::isnan(CFL_factor) || std::isinf(CFL_factor)) {
        Igor::Warn("CFL_factor is invalid at time t={}: CFL_factor = {}", t, CFL_factor);
        return std::nullopt;
      }
      const Float dt = std::min(0.5 * grid_curr.min_delta() / CFL_factor, tend - t);

#pragma omp parallel for
      for (std::size_t i = 0; i < grid_curr.m_cells.size(); ++i) {
        const auto& curr_cell = grid_curr.m_cells[i];
        auto& next_cell       = grid_next.m_cells[i];

        if constexpr (DIM == 1) {
          const auto F_minus = apply_flux(grid_curr, curr_cell, LEFT);
          const auto F_plus  = apply_flux(grid_curr, curr_cell, RIGHT);
          const auto G_minus = apply_flux(grid_curr, curr_cell, BOTTOM);
          const auto G_plus  = apply_flux(grid_curr, curr_cell, TOP);

          assert(curr_cell.is_cartesian() && "Non-cartesian cell is not implemented.");
          const auto& curr_cell_value = std::get<CartesianValue<Float, DIM>>(curr_cell.value).value;

          assert(next_cell.is_cartesian() && "Non-cartesian cell is not implemented.");
          auto& next_cell_value = std::get<CartesianValue<Float, DIM>>(next_cell.value).value;

          // clang-format off
          next_cell_value(0) = curr_cell_value(0) -
                               (dt / curr_cell.dx) * (F_plus - F_minus) -
                               (dt / curr_cell.dy) * (G_plus - G_minus);
          // clang-format on
        } else {
          Igor::Todo("Not implemented for DIM={}", DIM);
        }
      }

      // Update time
      t += dt;

      std::swap(grid_curr.m_cells, grid_next.m_cells);

      if (!grid_writer.write_data(grid_curr) || !time_writer.write_data(t)) {
        return std::nullopt;
      }
    }

    return grid_curr;
  }
};  // namespace Zap::CellBased

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_SOLVER_HPP_
