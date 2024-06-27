#ifndef ZAP_CELL_BASED_SOLVER_HPP_
#define ZAP_CELL_BASED_SOLVER_HPP_

#include "CellBased/Cell.hpp"

namespace Zap::CellBased {

template <typename FluxX, typename FluxY>
class Solver {
  FluxX m_numerical_flux_x;
  FluxY m_numerical_flux_y;

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

    if (!grid_writer.write_data(grid) || !time_writer.write_data(Float{0})) {
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
        auto& curr_cell = grid_curr.m_cells[i];
        auto& next_cell = grid_next.m_cells[i];
        assert((curr_cell.left_idx != CartesianCell<Float, DIM>::NULL_INDEX));
        assert((curr_cell.right_idx != CartesianCell<Float, DIM>::NULL_INDEX));
        assert((curr_cell.bottom_idx != CartesianCell<Float, DIM>::NULL_INDEX));
        assert((curr_cell.top_idx != CartesianCell<Float, DIM>::NULL_INDEX));

        if constexpr (DIM == 1) {
          // clang-format off
          const auto F_minus = m_numerical_flux_x(grid_curr.m_cells[curr_cell.left_idx].value[0],
                                                  curr_cell.value[0]);
          const auto F_plus  = m_numerical_flux_x(curr_cell.value[0],
                                                  grid_curr.m_cells[curr_cell.right_idx].value[0]);
          const auto G_minus = m_numerical_flux_x(grid_curr.m_cells[curr_cell.bottom_idx].value[0],
                                                  curr_cell.value[0]);
          const auto G_plus = m_numerical_flux_x(curr_cell.value[0],
                                                 grid_curr.m_cells[curr_cell.top_idx].value[0]);
          next_cell.value[0] = curr_cell.value[0] -
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

      if (!grid_writer.write_data(grid) || !time_writer.write_data(Float{0})) {
        return std::nullopt;
      }
    }

    return grid_curr;
  }
};

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_SOLVER_HPP_
