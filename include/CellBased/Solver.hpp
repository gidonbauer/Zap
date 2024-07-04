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

      auto get_flux = [&grid_curr](auto numerical_flux,
                                   size_t idx,
                                   const auto& cell,
                                   bool cell_first) noexcept {
        if constexpr (DIM == 1) {
          if (grid_curr.is_cell(idx)) {
            if (cell_first) {
              return numerical_flux(cell.value[0], grid_curr.m_cells[idx].value[0]);
            } else {
              return numerical_flux(grid_curr.m_cells[idx].value[0], cell.value[0]);
            }
          } else {
            switch (idx) {
              case SAME_VALUE_INDEX:
                return numerical_flux(cell.value[0], cell.value[0]);
              case ZERO_FLUX_INDEX:
                return static_cast<Float>(0);
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
      };

#pragma omp parallel for
      for (std::size_t i = 0; i < grid_curr.m_cells.size(); ++i) {
        const auto& curr_cell = grid_curr.m_cells[i];
        auto& next_cell       = grid_next.m_cells[i];

        if constexpr (DIM == 1) {
          const auto F_minus = get_flux(m_numerical_flux_x, curr_cell.left_idx, curr_cell, false);
          const auto F_plus  = get_flux(m_numerical_flux_x, curr_cell.right_idx, curr_cell, true);
          const auto G_minus = get_flux(m_numerical_flux_y, curr_cell.bottom_idx, curr_cell, false);
          const auto G_plus  = get_flux(m_numerical_flux_y, curr_cell.top_idx, curr_cell, true);
          // clang-format off
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

      if (!grid_writer.write_data(grid_curr) || !time_writer.write_data(t)) {
        return std::nullopt;
      }
    }

    return grid_curr;
  }
};

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_SOLVER_HPP_
