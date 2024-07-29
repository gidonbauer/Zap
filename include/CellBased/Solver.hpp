#ifndef ZAP_CELL_BASED_SOLVER_HPP_
#define ZAP_CELL_BASED_SOLVER_HPP_

#include <numbers>

#include "CellBased/Grid.hpp"

#ifndef OLD_SOLVER

namespace Zap::CellBased {

template <typename A, typename B>
class Solver {
  A m_A;
  B m_B;

 public:
  // -----------------------------------------------------------------------------------------------
  constexpr Solver(A a, B b)
      : m_A(std::move(a)),
        m_B(std::move(b)) {}

  // -----------------------------------------------------------------------------------------------
  template <typename Float, size_t DIM, typename GridWriter, typename TimeWriter>
  [[nodiscard]] constexpr auto
  solve(Grid<Float, DIM> grid,
        Float tend,
        GridWriter& grid_writer,
        TimeWriter& time_writer,
        Float CFL_safety_factor = 0.5) noexcept -> std::optional<Grid<Float, DIM>> {
    using VecOrScalar = std::conditional_t<DIM == 1, Float, Eigen::Vector<Float, DIM>>;
    using MatOrScalar = std::conditional_t<DIM == 1, Float, Eigen::Matrix<Float, DIM, DIM>>;

    assert(CFL_safety_factor > 0 && CFL_safety_factor <= 1);  // Must be in (0, 1]
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
      const Float dt = std::min(CFL_safety_factor * grid_curr.min_delta() / CFL_factor, tend - t);

      // TODO: Move shock in next grid
      // TODO: Remove old shock in next grid
      // TODO: Update values

      for (size_t i = 0; i < grid_curr.m_cells.size(); ++i) {
        const auto& curr_cell = grid_curr.m_cells[i];
        auto& next_cell       = grid_next.m_cells[i];

        auto get_u_mid = [](const Eigen::Vector<Float, DIM>& lhs_value,
                            const Eigen::Vector<Float, DIM>& rhs_value) -> VecOrScalar {
          if constexpr (DIM == 1) {
            return (lhs_value(0) + rhs_value(0)) / 2;
          } else {
            return (lhs_value + rhs_value) / 2;
          }
        };

        auto linearized_xy_flux_impl = [this,
                                        get_u_mid](const Eigen::Vector<Float, DIM>& curr_value,
                                                   const Eigen::Vector<Float, DIM>& other_value,
                                                   Float factor,
                                                   Side from) -> VecOrScalar {
          const auto u_mid     = get_u_mid(curr_value, other_value);
          MatOrScalar eig_vals = (from & (LEFT | RIGHT)) > 0 ? m_A.eig_vals(u_mid)  //
                                                             : m_B.eig_vals(u_mid);
          MatOrScalar eig_vecs = (from & (LEFT | RIGHT)) > 0 ? m_A.eig_vecs(u_mid)  //
                                                             : m_B.eig_vecs(u_mid);

          // Scalar case -------------------------------------------------------------------------
          if constexpr (DIM == 1) {
            if (std::abs(eig_vals) < 1e-8) { return Float{0}; }

            assert(std::abs(eig_vecs) > 1e-8);
            if ((from & (LEFT | BOTTOM)) > 0) {
              const auto A_plus = eig_vals > 0 ? eig_vecs * eig_vals / eig_vecs : 0;
              return -A_plus * factor * (curr_value(0) - other_value(0));
            } else {
              const auto A_minus = eig_vals < 0 ? eig_vecs * eig_vals / eig_vecs : 0;
              return A_minus * factor * (other_value(0) - curr_value(0));
            }
          }
          // Vector case -------------------------------------------------------------------------
          else {
            bool all_zero = true;
            for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
              if (std::abs(eig_vals(i, i)) >= 1e-6) { all_zero = false; }
            }
            if (all_zero) { return VecOrScalar::Zero(); }

            assert(std::abs(eig_vecs.determinant()) >= 1e-8);
            if ((from & (LEFT | BOTTOM)) > 0) {
              for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
                eig_vals(i, i) = eig_vals(i, i) > 0 ? eig_vals(i, i) : 0;
              }
              const auto A_plus = (eig_vecs * eig_vals * eig_vecs.inverse()).eval();
              return -A_plus * factor * (curr_value - other_value);
            } else {
              for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
                eig_vals(i, i) = eig_vals(i, i) < 0 ? eig_vals(i, i) : 0;
              }
              const auto A_minus = (eig_vecs * eig_vals * eig_vecs.inverse()).eval();
              return A_minus * factor * (other_value - curr_value);
            }
          }
        };

        // curr_cell must be cartesian, cut cells are handled differently
        auto linearized_xy_flux = [dt, linearized_xy_flux_impl](const Cell<Float, DIM>& curr_cell,
                                                                const Cell<Float, DIM>& other_cell,
                                                                Side from) -> VecOrScalar {
          assert(from == LEFT || from == BOTTOM || from == RIGHT || from == TOP);
          assert(curr_cell.is_cartesian());
          assert(other_cell.is_cartesian() || other_cell.is_cut());

          if (other_cell.is_cartesian()) {
            const auto factor = (from & (LEFT | RIGHT)) > 0 ? (dt / curr_cell.dx)  //
                                                            : (dt / curr_cell.dy);
            return linearized_xy_flux_impl(
                curr_cell.get_cartesian().value, other_cell.get_cartesian().value, factor, from);
          } else {
            switch (other_cell.get_cut().type) {
              // - BOTTOM_LEFT ---------------------------------------------------------------------
              case CutType::BOTTOM_LEFT:
                {
                  switch (from) {
                    case LEFT:
                      {
                        const auto factor = dt / curr_cell.dx;
                        return linearized_xy_flux_impl(curr_cell.get_cartesian().value,
                                                       other_cell.get_cut().right_value,
                                                       factor,
                                                       from);
                      }
                      break;
                    case RIGHT:
                      {
                        Igor::Todo("Cut type BOTTOM_LEFT not implemented for RIGHT yet.");
                      }
                      break;
                    case BOTTOM:
                      {
                        const auto factor = dt / curr_cell.dy;
                        return linearized_xy_flux_impl(curr_cell.get_cartesian().value,
                                                       other_cell.get_cut().right_value,
                                                       factor,
                                                       from);
                      }
                      break;
                    case TOP:
                      {
                        Igor::Todo("Cut type BOTTOM_LEFT not implemented for TOP yet.");
                      }
                      break;
                    case ALL: std::unreachable();  // To silence warning
                  }
                }
                break;

              // - BOTTOM_RIGHT --------------------------------------------------------------------
              case CutType::BOTTOM_RIGHT:
                {
                  switch (from) {
                    case LEFT:
                      {
                        Igor::Todo("Cut type BOTTOM_RIGHT not implemented for LEFT yet.");
                      }
                      break;
                    case RIGHT:
                      {
                        Igor::Todo("Cut type BOTTOM_RIGHT not implemented for RIGHT yet.");
                      }
                      break;
                    case BOTTOM:
                      {
                        Igor::Todo("Cut type BOTTOM_RIGHT not implemented for BOTTOM yet.");
                      }
                      break;
                    case TOP:
                      {
                        Igor::Todo("Cut type BOTTOM_RIGHT not implemented for TOP yet.");
                      }
                      break;
                    case ALL: std::unreachable();  // To silence warning
                  }
                }
                break;

              // - TOP_RIGHT -----------------------------------------------------------------------
              case CutType::TOP_RIGHT:
                {
                  switch (from) {
                    case LEFT:
                      {
                        Igor::Todo("Cut type TOP_RIGHT not implemented for LEFT yet.");
                      }
                      break;
                    case RIGHT:
                      {
                        const auto factor = dt / curr_cell.dx;
                        return linearized_xy_flux_impl(curr_cell.get_cartesian().value,
                                                       other_cell.get_cut().left_value,
                                                       factor,
                                                       from);
                      }
                      break;
                    case BOTTOM:
                      {
                        Igor::Todo("Cut type TOP_RIGHT not implemented for BOTTOM yet.");
                      }
                      break;
                    case TOP:
                      {
                        const auto factor = dt / curr_cell.dy;
                        return linearized_xy_flux_impl(curr_cell.get_cartesian().value,
                                                       other_cell.get_cut().left_value,
                                                       factor,
                                                       from);
                      }
                      break;
                    case ALL: std::unreachable();  // To silence warning
                  }
                }
                break;

              // - TOP_LEFT ------------------------------------------------------------------------
              case CutType::TOP_LEFT:
                {
                  switch (from) {
                    case LEFT:
                      {
                        Igor::Todo("Cut type TOP_LEFT not implemented for LEFT yet.");
                      }
                      break;
                    case RIGHT:
                      {
                        Igor::Todo("Cut type TOP_LEFT not implemented for RIGHT yet.");
                      }
                      break;
                    case BOTTOM:
                      {
                        Igor::Todo("Cut type TOP_LEFT not implemented for BOTTOM yet.");
                      }
                      break;
                    case TOP:
                      {
                        Igor::Todo("Cut type TOP_LEFT not implemented for TOP yet.");
                      }
                      break;
                    case ALL: std::unreachable();  // To silence warning
                  }
                }
                break;

              // - MIDDLE_HORI ---------------------------------------------------------------------
              case CutType::MIDDLE_HORI:
                {
                  switch (from) {
                    case LEFT:
                      {
                        Igor::Todo("Cut type MIDDLE_HORI not implemented for LEFT yet.");
                      }
                      break;
                    case RIGHT:
                      {
                        const auto dy_lower_half = other_cell.get_cut().y2_cut - other_cell.y_min;
                        const auto factor_lower_half =
                            (dy_lower_half / other_cell.dy) * dt / other_cell.dx;
                        const auto factor_upper_half =
                            (1 - (dy_lower_half / other_cell.dy)) * dt / other_cell.dx;
                        return linearized_xy_flux_impl(curr_cell.get_cartesian().value,
                                                       other_cell.get_cut().left_value,
                                                       factor_lower_half,
                                                       from) +
                               linearized_xy_flux_impl(curr_cell.get_cartesian().value,
                                                       other_cell.get_cut().right_value,
                                                       factor_upper_half,
                                                       from);
                      }
                      break;
                    case BOTTOM:
                      {
                        const auto factor = dt / curr_cell.dy;
                        return linearized_xy_flux_impl(curr_cell.get_cartesian().value,
                                                       other_cell.get_cut().right_value,
                                                       factor,
                                                       from);
                      }
                      break;
                    case TOP:
                      {
                        const auto factor = dt / curr_cell.dy;
                        return linearized_xy_flux_impl(curr_cell.get_cartesian().value,
                                                       other_cell.get_cut().left_value,
                                                       factor,
                                                       from);
                      }
                      break;
                    case ALL: std::unreachable();  // To silence warning
                  }
                }
                break;

              // - MIDDLE_VERT ---------------------------------------------------------------------
              case CutType::MIDDLE_VERT:
                {
                  switch (from) {
                    case LEFT:
                      {
                        const auto factor = dt / curr_cell.dx;
                        return linearized_xy_flux_impl(curr_cell.get_cartesian().value,
                                                       other_cell.get_cut().right_value,
                                                       factor,
                                                       from);
                      }
                      break;
                    case RIGHT:
                      {
                        const auto factor = dt / curr_cell.dx;
                        return linearized_xy_flux_impl(curr_cell.get_cartesian().value,
                                                       other_cell.get_cut().left_value,
                                                       factor,
                                                       from);
                      }
                      break;
                    case BOTTOM:
                      {
                        Igor::Todo("Cut type MIDDLE_VERT not implemented for BOTTOM yet.");
                      }
                      break;
                    case TOP:
                      {
                        const auto dx_left_half = other_cell.get_cut().x1_cut - other_cell.x_min;
                        const auto factor_left_half =
                            (dx_left_half / other_cell.dx) * dt / other_cell.dy;
                        const auto factor_right_half =
                            (1 - (dx_left_half / other_cell.dx)) * dt / other_cell.dy;
                        return linearized_xy_flux_impl(curr_cell.get_cartesian().value,
                                                       other_cell.get_cut().left_value,
                                                       factor_left_half,
                                                       from) +
                               linearized_xy_flux_impl(curr_cell.get_cartesian().value,
                                                       other_cell.get_cut().right_value,
                                                       factor_right_half,
                                                       from);
                      }
                      break;
                    case ALL: std::unreachable();  // To silence warning
                  }
                }
                break;

              // -----------------------------------------------------------------------------------
              default:
                Igor::Panic(
                    "Unknown cut type with value {}",
                    static_cast<std::underlying_type_t<CutType>>(other_cell.get_cut().type));
                std::unreachable();
            }

            Igor::Panic("Unreachable.");
            std::unreachable();
          }
        };

        if (curr_cell.is_cartesian()) {
          // Left side: U_(i, j-1) -> U_(i, j)
          assert(grid_curr.is_cell(curr_cell.left_idx));
          const auto F_minus = linearized_xy_flux(curr_cell,  //
                                                  grid_curr.m_cells[curr_cell.left_idx],
                                                  LEFT);

          // Right side: U_(i, j) -> U_(i, j+1)
          assert(grid_curr.is_cell(curr_cell.right_idx));
          const auto F_plus = linearized_xy_flux(curr_cell,  //
                                                 grid_curr.m_cells[curr_cell.right_idx],
                                                 RIGHT);

          // Bottom side: U_(i-1, j) -> U_(i, j)
          assert(grid_curr.is_cell(curr_cell.bottom_idx));
          const auto G_minus = linearized_xy_flux(curr_cell,  //
                                                  grid_curr.m_cells[curr_cell.bottom_idx],
                                                  BOTTOM);

          // Top side: U_(i, j) -> U_(i+1, j)
          assert(grid_curr.is_cell(curr_cell.top_idx));
          const auto G_plus = linearized_xy_flux(curr_cell,  //
                                                 grid_curr.m_cells[curr_cell.top_idx],
                                                 TOP);

          if (next_cell.is_cartesian()) {
            if constexpr (DIM == 1) {
              next_cell.get_cartesian().value(0) = curr_cell.get_cartesian().value(0) -  //
                                                   (F_plus - F_minus) -                  //
                                                   (G_plus - G_minus);
            } else {
              next_cell.get_cartesian().value = curr_cell.get_cartesian().value -  //
                                                (F_plus - F_minus) -               //
                                                (G_plus - G_minus);
            }
          } else {
            Igor::Todo("next_cell non-cartesian is not implemented yet.");
          }
        } else if (curr_cell.is_cut()) {
          // - Handle internal interface -----------------------------------------------------------
          const Eigen::Vector<Float, 2> cut_vector{
              curr_cell.get_cut().x2_cut - curr_cell.get_cut().x1_cut,
              curr_cell.get_cut().y2_cut - curr_cell.get_cut().y1_cut};
          const auto cut_angle = std::acos(cut_vector(1) / cut_vector.norm());
          Igor::Debug("cut_angle = {:.6f} rad ({:.6f}°)",
                      cut_angle,
                      cut_angle * 180 / std::numbers::pi_v<Float>);
          Igor::Debug("cos(cut_angle) = {}", std::cos(cut_angle));
          Igor::Debug("-sin(cut_angle) = {}", -std::sin(cut_angle));

          const Eigen::Vector<Float, 2> norm_vector{std::cos(cut_angle), std::sin(cut_angle)};

          const auto u_mid =
              get_u_mid(curr_cell.get_cut().left_value, curr_cell.get_cut().right_value);

          const MatOrScalar eta_mat =
              std::cos(cut_angle) * m_A.mat(u_mid) - std::sin(cut_angle) * m_B.mat(u_mid);

          if constexpr (DIM == 1) {
            IGOR_DEBUG_PRINT(eta_mat);

            const auto eig_vals = eta_mat;
            const auto eig_vecs = MatOrScalar{1};
            const auto alpha =
                curr_cell.get_cut().left_value(0) - curr_cell.get_cut().right_value(0);
            const auto wave_length = eig_vals * dt;

            IGOR_DEBUG_PRINT(eig_vals);
            IGOR_DEBUG_PRINT(eig_vecs);
            IGOR_DEBUG_PRINT(alpha);
            Igor::Debug("Wave length = {}", wave_length);

            const Eigen::Vector<Float, 2> start_point_1 =
                Eigen::Vector<Float, 2>{curr_cell.get_cut().x1_cut, curr_cell.get_cut().y1_cut};
            const Eigen::Vector<Float, 2> end_point_1 = start_point_1 + norm_vector * wave_length;

            const Eigen::Vector<Float, 2> start_point_2 =
                Eigen::Vector<Float, 2>{curr_cell.get_cut().x2_cut, curr_cell.get_cut().y2_cut};
            const Eigen::Vector<Float, 2> end_point_2 = start_point_2 + norm_vector * wave_length;

            IGOR_DEBUG_PRINT(start_point_1);
            IGOR_DEBUG_PRINT(end_point_1);
            IGOR_DEBUG_PRINT(start_point_2);
            IGOR_DEBUG_PRINT(end_point_2);

            {
              std::stringstream s{};
              s << curr_cell;
              Igor::Debug("curr_cell = {}", s.str());
            }
            {
              std::stringstream s{};
              s << grid_curr.m_cells.at(curr_cell.right_idx);
              Igor::Debug("right_cell = {}", s.str());
            }
            {
              std::stringstream s{};
              s << grid_curr.m_cells.at(curr_cell.top_idx);
              Igor::Debug("top_cell = {}", s.str());
            }
            {
              std::stringstream s{};
              s << grid_curr.m_cells.at(grid_curr.m_cells.at(curr_cell.top_idx).right_idx);
              Igor::Debug("diag_cell = {}", s.str());
            }
          } else {
            Igor::Todo("DIM = {} is not implemented yet.", DIM);
          }
          return std::nullopt;
          // Igor::Todo("Cut cells are not implemented yet.");
        } else {
          Igor::Warn("Unknown cell type with variant index {}", curr_cell.value.index());
          return std::nullopt;
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
