#ifndef ZAP_CELL_BASED_SOLVER_HPP_
#define ZAP_CELL_BASED_SOLVER_HPP_

#include "CellBased/Grid.hpp"

#ifndef OLD_SOLVER

namespace Zap::CellBased {

template <typename EigValsX, typename EigVecsX, typename EigValsY, typename EigVecsY>
class Solver {
  EigValsX m_eig_vals_x;
  EigVecsX m_eig_vecs_x;
  EigValsY m_eig_vals_y;
  EigVecsY m_eig_vecs_y;

 public:
  // -----------------------------------------------------------------------------------------------
  constexpr Solver(EigValsX eig_vals_x,
                   EigVecsX eig_vecs_x,
                   EigValsY eig_vals_y,
                   EigVecsY eig_vecs_y)
      : m_eig_vals_x(std::move(eig_vals_x)),
        m_eig_vecs_x(std::move(eig_vecs_x)),
        m_eig_vals_y(std::move(eig_vals_y)),
        m_eig_vecs_y(std::move(eig_vecs_y)) {}

  // -----------------------------------------------------------------------------------------------
  template <typename Float, size_t DIM, typename GridWriter, typename TimeWriter>
  [[nodiscard]] constexpr auto
  solve(Grid<Float, DIM> grid,
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
      constexpr Float CFL_safety_factor = 0.25;  // Must be in (0, 1]
      const Float dt = std::min(CFL_safety_factor * grid_curr.min_delta() / CFL_factor, tend - t);

      // TODO: Move shock in next grid
      // TODO: Remove old shock in next grid
      // TODO: Update values

      for (size_t i = 0; i < grid_curr.m_cells.size(); ++i) {
        const auto& curr_cell = grid_curr.m_cells[i];
        auto& next_cell       = grid_next.m_cells[i];
        assert(curr_cell.is_cartesian());
        assert(next_cell.is_cartesian());

        using FluxType       = std::conditional_t<DIM == 1, Float, Eigen::Vector<Float, DIM>>;
        using EigType        = std::conditional_t<DIM == 1, Float, Eigen::Matrix<Float, DIM, DIM>>;
        auto linearized_flux = [](const Cell<Float, DIM>& curr_cell,
                                  const Cell<Float, DIM>& other_cell,
                                  auto get_eig_vals,
                                  auto get_eig_vecs,
                                  bool from_left_or_bottom) -> FluxType {
          assert(curr_cell.is_cartesian());
          assert(other_cell.is_cartesian());

          const auto u_mid = [&] -> FluxType {
            if constexpr (DIM == 1) {
              return (curr_cell.get_cartesian().value(0) + other_cell.get_cartesian().value(0)) / 2;
            } else {
              return (curr_cell.get_cartesian().value + other_cell.get_cartesian().value) / 2;
            }
          }();
          EigType eig_vals = get_eig_vals(u_mid);
          EigType eig_vecs = get_eig_vecs(u_mid);

          if constexpr (DIM == 1) {
            if (std::abs(eig_vals) < 1e-8) { return Float{0}; }

            assert(std::abs(eig_vecs) > 1e-8);
            if (from_left_or_bottom) {
              const auto A_plus = eig_vals > 0 ? eig_vecs * eig_vals / eig_vecs : 0;
              return -A_plus *
                     (curr_cell.get_cartesian().value(0) - other_cell.get_cartesian().value(0));
            } else {
              const auto A_minus = eig_vals < 0 ? eig_vecs * eig_vals / eig_vecs : 0;
              return A_minus *
                     (other_cell.get_cartesian().value(0) - curr_cell.get_cartesian().value(0));
            }
          } else {
            bool all_zero = true;
            for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
              if (std::abs(eig_vals(i, i)) >= 1e-6) { all_zero = false; }
            }
            if (all_zero) { return FluxType::Zero(); }

            if (!(std::abs(eig_vecs.determinant()) >= 1e-8)) {
              std::cerr << "u_mid = " << u_mid.transpose() << '\n';
              std::cerr << "eig_vals =\n" << eig_vals << '\n';
              std::cerr << "eig_vecs =\n" << eig_vecs << '\n';
              std::cerr << "det(eig_vecs) = " << eig_vecs.determinant() << '\n';
            }
            assert(std::abs(eig_vecs.determinant()) >= 1e-8);
            if (from_left_or_bottom) {
              for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
                eig_vals(i, i) = eig_vals(i, i) > 0 ? eig_vals(i, i) : 0;
              }

              // const auto A_plus = (eigvecs * eigvals * eigvecs.transpose()).eval();
              const auto A_plus = (eig_vecs * eig_vals * eig_vecs.inverse()).eval();
              return (-A_plus *
                      (curr_cell.get_cartesian().value - other_cell.get_cartesian().value))
                  .eval();
            } else {
              for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
                eig_vals(i, i) = eig_vals(i, i) < 0 ? eig_vals(i, i) : 0;
              }

              // const auto A_minus = (eigvecs * eigvals * eigvecs.transpose()).eval();
              const auto A_minus = (eig_vecs * eig_vals * eig_vecs.inverse()).eval();
              return (A_minus *
                      (other_cell.get_cartesian().value - curr_cell.get_cartesian().value))
                  .eval();
            }
          }
        };

        // Left side: U_(i, j-1) -> U_(i, j)
        assert(grid_curr.is_cell(curr_cell.left_idx));
        const auto F_minus = linearized_flux(
            curr_cell, grid_curr.m_cells[curr_cell.left_idx], m_eig_vals_x, m_eig_vecs_x, true);

        // Right side: U_(i, j) -> U_(i, j+1)
        assert(grid_curr.is_cell(curr_cell.right_idx));
        const auto F_plus = linearized_flux(
            curr_cell, grid_curr.m_cells[curr_cell.right_idx], m_eig_vals_x, m_eig_vecs_x, false);

        // Bottom side: U_(i-1, j) -> U_(i, j)
        assert(grid_curr.is_cell(curr_cell.bottom_idx));
        const auto G_minus = linearized_flux(
            curr_cell, grid_curr.m_cells[curr_cell.bottom_idx], m_eig_vals_y, m_eig_vecs_y, true);

        // Top side: U_(i, j) -> U_(i+1, j)
        assert(grid_curr.is_cell(curr_cell.top_idx));
        const auto G_plus = linearized_flux(
            curr_cell, grid_curr.m_cells[curr_cell.top_idx], m_eig_vals_y, m_eig_vecs_y, false);

        // Left side: U_(i, j-1) -> U_(i, j)
        // const auto F_minus = [&] -> FluxType {
        //   assert(grid_curr.is_cell(curr_cell.left_idx));
        //   const auto& left_cell = grid_curr.m_cells[curr_cell.left_idx];
        //   assert(left_cell.is_cartesian());
        //   const auto u_mid = [&] -> FluxType {
        //     if constexpr (DIM == 1) {
        //       return (curr_cell.get_cartesian().value(0) + left_cell.get_cartesian().value(0)) /
        //       2;
        //     } else {
        //       return (curr_cell.get_cartesian().value + left_cell.get_cartesian().value) / 2;
        //     }
        //   }();

        //   auto eigvals = m_eig_vals_x(u_mid);
        //   auto eigvecs = m_eig_vecs_x(u_mid);

        //   if constexpr (DIM == 1) {
        //     const auto A_plus = eigvals > 0 ? eigvecs * eigvals / eigvecs : 0;
        //     return -A_plus *
        //            (curr_cell.get_cartesian().value(0) - left_cell.get_cartesian().value(0));
        //   } else {
        //     for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        //       eigvals(i, i) = eigvals(i, i) > 0 ? eigvals(i, i) : 0;
        //     }
        //     // const auto A_plus = (eigvecs * eigvals * eigvecs.transpose()).eval();
        //     const auto A_plus = (eigvecs * eigvals * eigvecs.inverse()).eval();
        //     return (-A_plus * (curr_cell.get_cartesian().value -
        //     left_cell.get_cartesian().value))
        //         .eval();
        //   }
        // }();

        // Left side: U_(i, j) -> U_(i, j+1)
        // const auto F_plus = [&] -> FluxType {
        //   assert(grid_curr.is_cell(curr_cell.right_idx));
        //   const auto& right_cell = grid_curr.m_cells[curr_cell.right_idx];
        //   assert(right_cell.is_cartesian());
        //   const auto u_mid = [&] -> FluxType {
        //     if constexpr (DIM == 1) {
        //       return (curr_cell.get_cartesian().value(0) + right_cell.get_cartesian().value(0)) /
        //       2;
        //     } else {
        //       return (curr_cell.get_cartesian().value + right_cell.get_cartesian().value) / 2;
        //     }
        //   }();

        //   auto eigvals = m_eig_vals_x(u_mid);
        //   auto eigvecs = m_eig_vecs_x(u_mid);

        //   if constexpr (DIM == 1) {
        //     const auto A_minus = eigvals < 0 ? eigvecs * eigvals / eigvecs : 0;
        //     return A_minus *
        //            (right_cell.get_cartesian().value(0) - curr_cell.get_cartesian().value(0));
        //   } else {
        //     for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        //       eigvals(i, i) = eigvals(i, i) < 0 ? eigvals(i, i) : 0;
        //     }
        //     // const auto A_minus = (eigvecs * eigvals * eigvecs.transpose()).eval();
        //     const auto A_minus = (eigvecs * eigvals * eigvecs.inverse()).eval();
        //     return (A_minus * (right_cell.get_cartesian().value -
        //     curr_cell.get_cartesian().value))
        //         .eval();
        //   }
        // }();

        // Left side: U_(i-1, j) -> U_(i, j)
        // const auto G_minus = [&] -> FluxType {
        //   assert(grid_curr.is_cell(curr_cell.bottom_idx));
        //   const auto& bottom_cell = grid_curr.m_cells[curr_cell.bottom_idx];
        //   assert(bottom_cell.is_cartesian());
        //   const auto u_mid = [&] -> FluxType {
        //     if constexpr (DIM == 1) {
        //       return (curr_cell.get_cartesian().value(0) + bottom_cell.get_cartesian().value(0))
        //       /
        //               2;
        //     } else {
        //       return (curr_cell.get_cartesian().value + bottom_cell.get_cartesian().value) / 2;
        //     }
        //   }();

        //   auto eigvals = m_eig_vals_y(u_mid);
        //   auto eigvecs = m_eig_vecs_y(u_mid);

        //   if constexpr (DIM == 1) {
        //     const auto A_plus = eigvals > 0 ? eigvecs * eigvals / eigvecs : 0;
        //     return -A_plus *
        //             (curr_cell.get_cartesian().value(0) - bottom_cell.get_cartesian().value(0));
        //   } else {
        //     for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        //       eigvals(i, i) = eigvals(i, i) > 0 ? eigvals(i, i) : 0;
        //     }
        //     // const auto A_plus = (eigvecs * eigvals * eigvecs.transpose()).eval();
        //     const auto A_plus = (eigvecs * eigvals * eigvecs.inverse()).eval();
        //     return (-A_plus * (curr_cell.get_cartesian().value -
        //     bottom_cell.get_cartesian().value))
        //         .eval();
        //   }
        // }();

        // Left side: U_(i, j) -> U_(i+1, j)
        // const auto G_plus = [&] -> FluxType {
        //   assert(grid_curr.is_cell(curr_cell.top_idx));
        //   const auto& top_cell = grid_curr.m_cells[curr_cell.top_idx];
        //   assert(top_cell.is_cartesian());
        //   const auto u_mid = [&] -> FluxType {
        //     if constexpr (DIM == 1) {
        //       return (curr_cell.get_cartesian().value(0) + top_cell.get_cartesian().value(0)) /
        //       2;
        //     } else {
        //       return (curr_cell.get_cartesian().value + top_cell.get_cartesian().value) / 2;
        //     }
        //   }();

        //   auto eigvals = m_eig_vals_y(u_mid);
        //   auto eigvecs = m_eig_vecs_y(u_mid);

        //   if constexpr (DIM == 1) {
        //     const auto A_minus = eigvals < 0 ? eigvecs * eigvals / eigvecs : 0;
        //     return A_minus *
        //             (top_cell.get_cartesian().value(0) - curr_cell.get_cartesian().value(0));
        //   } else {
        //     for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        //       eigvals(i, i) = eigvals(i, i) < 0 ? eigvals(i, i) : 0;
        //     }
        //     // const auto A_minus = (eigvecs * eigvals * eigvecs.transpose()).eval();
        //     const auto A_minus = (eigvecs * eigvals * eigvecs.inverse()).eval();
        //     return (A_minus * (top_cell.get_cartesian().value - curr_cell.get_cartesian().value))
        //         .eval();
        //   }
        // }();

        assert(curr_cell.is_cartesian() && next_cell.is_cartesian());
        if constexpr (DIM == 1) {
          next_cell.get_cartesian().value(0) = curr_cell.get_cartesian().value(0) -
                                               (dt / curr_cell.dx) * (F_plus - F_minus) -
                                               (dt / curr_cell.dy) * (G_plus - G_minus);
        } else {
          next_cell.get_cartesian().value = curr_cell.get_cartesian().value -
                                            (dt / curr_cell.dx) * (F_plus - F_minus) -
                                            (dt / curr_cell.dy) * (G_plus - G_minus);
        }
      }

      // Update time
      t += dt;

      grid_curr = grid_next;

      if (!grid_writer.write_data(grid_curr) || !time_writer.write_data(t)) { return std::nullopt; }
    }

    return grid_curr;
  }
};

}  // namespace Zap::CellBased

#else

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

    if (!grid.is_cell(idx)) { return apply_boundary_flux(cell, side, idx); }
    const auto& other_cell = grid.m_cells[idx];
    assert(cell.is_cartesian() || cell.is_cut() && "Invalid cell type.");
    assert(other_cell.is_cartesian() || other_cell.is_cut() && "Invalid cell type.");

    // This cell cartesian and other cell cartesian
    if (cell.is_cartesian() && other_cell.is_cartesian()) {
      const auto& cell_value       = cell.get_cartesian().value;
      const auto& other_cell_value = other_cell.get_cartesian().value;

      switch (side) {
        case LEFT:   return m_numerical_flux_x(other_cell_value(0), cell_value(0));
        case RIGHT:  return m_numerical_flux_x(cell_value(0), other_cell_value(0));

        case BOTTOM: return m_numerical_flux_y(other_cell_value(0), cell_value(0));
        case TOP:    return m_numerical_flux_y(cell_value(0), other_cell_value(0));

        default:     Igor::Panic("Must be BOTTOM, RIGHT, TOP, or LEFT."); std::unreachable();
      }
    }
    // This cell cartesian and other cell cut
    else if (cell.is_cartesian() && other_cell.is_cut()) {
      std::stringstream s{};
      s << cell;
      Igor::Debug("cell = {}", s.str());
      s.str("");
      s << grid.m_cells[idx];
      Igor::Debug("other_cell = {}", s.str());

      Igor::Todo("Flux for cartesian to cut cells is not implemented yet.");
      std::unreachable();
    }
    // This cell cut and other cell cartesian
    else if (cell.is_cut() && other_cell.is_cartesian()) {
      std::stringstream s{};
      s << cell;
      Igor::Debug("cell = {}", s.str());
      s.str("");
      s << grid.m_cells[idx];
      Igor::Debug("other_cell = {}", s.str());

      Igor::Todo("Flux for cut to cartesian cells is not implemented yet.");
      std::unreachable();
    }
    // This cell cut and other cell cut
    else if (cell.is_cut() && other_cell.is_cut()) {
      std::stringstream s{};
      s << cell;
      Igor::Debug("cell = {}", s.str());
      s.str("");
      s << grid.m_cells[idx];
      Igor::Debug("other_cell = {}", s.str());

      Igor::Todo("Flux for cut to cut cells is not implemented yet.");
      std::unreachable();
    }

    Igor::Panic("Unreachable.");
    std::unreachable();
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
          if (curr_cell.is_cartesian()) {
            const auto F_minus = apply_flux(grid_curr, curr_cell, LEFT);
            const auto F_plus  = apply_flux(grid_curr, curr_cell, RIGHT);
            const auto G_minus = apply_flux(grid_curr, curr_cell, BOTTOM);
            const auto G_plus  = apply_flux(grid_curr, curr_cell, TOP);

            const auto& curr_cell_value = curr_cell.get_cartesian().value;

            assert(next_cell.is_cartesian() && "Non-cartesian cell is not implemented.");
            auto& next_cell_value = next_cell.get_cartesian().value;

            // clang-format off
            next_cell_value(0) = curr_cell_value(0) -
                                 (dt / curr_cell.dx) * (F_plus - F_minus) -
                                 (dt / curr_cell.dy) * (G_plus - G_minus);
            // clang-format on
          } else if (curr_cell.is_cut()) {
            switch (curr_cell.get_cut().type) {
              case CutType::BOTTOM_LEFT:  Igor::Todo("BOTTOM_LEFT"); break;
              case CutType::BOTTOM_RIGHT: Igor::Todo("BOTTOM_RIGHT"); break;
              case CutType::TOP_RIGHT:    Igor::Todo("TOP_RIGHT"); break;
              case CutType::TOP_LEFT:     Igor::Todo("TOP_LEFT"); break;
              case CutType::MIDDLE_HORI:  Igor::Todo("MIDDLE_HORI"); break;
              case CutType::MIDDLE_VERT:  Igor::Todo("MIDDLE_VERT"); break;
            }
          } else {
            Igor::Panic("Invalid cell type with variant index {}", curr_cell.value.index());
          }
        } else {
          Igor::Todo("Not implemented for DIM={}", DIM);
        }
      }

      // Update time
      t += dt;

      std::swap(grid_curr.m_cells, grid_next.m_cells);

      if (!grid_writer.write_data(grid_curr) || !time_writer.write_data(t)) { return std::nullopt; }
    }

    return grid_curr;
  }
};

}  // namespace Zap::CellBased

#endif  // OLD_SOLVER

#endif  // ZAP_CELL_BASED_SOLVER_HPP_
