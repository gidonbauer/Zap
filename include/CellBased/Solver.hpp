#ifndef ZAP_CELL_BASED_SOLVER_HPP_
#define ZAP_CELL_BASED_SOLVER_HPP_

#include "CellBased/Grid.hpp"

namespace Zap::CellBased {

template <typename FluxX, typename FluxY>
class Solver {
  FluxX m_numerical_flux_x;
  FluxY m_numerical_flux_y;

  template <typename Float, size_t DIM>
  [[nodiscard]] constexpr auto apply_flux(const Grid<Float, DIM>& grid,
                                          auto numerical_flux,
                                          size_t idx,
                                          const auto& cell,
                                          bool cell_first) const noexcept {
    if constexpr (DIM == 1) {
      if (grid.is_cell(idx)) {
        assert(cell.is_cartesian() && "Non-cartesian cell not implemented.");
        assert(grid.m_cells[idx].is_cartesian() && "Non-cartesian cell not implemented.");

        const auto& cell_value = std::get<CartesianValue<Float, DIM>>(cell.value).value;
        const auto& other_cell_value =
            std::get<CartesianValue<Float, DIM>>(grid.m_cells[idx].value).value;

        if (cell_first) {
          return numerical_flux(cell_value(0), other_cell_value(0));
        } else {
          return numerical_flux(other_cell_value(0), cell_value(0));
        }
      } else {
        switch (idx) {
          case SAME_VALUE_INDEX:
            {
              assert(cell.is_cartesian() && "Non-cartesian cell not implemented.");
              const auto& cell_value = std::get<CartesianValue<Float, DIM>>(cell.value).value;
              return numerical_flux(cell_value(0), cell_value(0));
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
      }
      std::unreachable();
    } else {
      Igor::Todo("Not implemented for DIM={}", DIM);
    }
  }

 public:
  constexpr Solver(FluxX F, FluxY G)
      : m_numerical_flux_x(std::move(F)),
        m_numerical_flux_y(std::move(G)) {}

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
          const auto F_minus =
              apply_flux(grid_curr, m_numerical_flux_x, curr_cell.left_idx, curr_cell, false);
          const auto F_plus =
              apply_flux(grid_curr, m_numerical_flux_x, curr_cell.right_idx, curr_cell, true);
          const auto G_minus =
              apply_flux(grid_curr, m_numerical_flux_y, curr_cell.bottom_idx, curr_cell, false);
          const auto G_plus =
              apply_flux(grid_curr, m_numerical_flux_y, curr_cell.top_idx, curr_cell, true);

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
};

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_SOLVER_HPP_
