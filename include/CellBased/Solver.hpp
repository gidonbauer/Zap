#ifndef ZAP_CELL_BASED_SOLVER_HPP_
#define ZAP_CELL_BASED_SOLVER_HPP_

#include <numbers>
#include <numeric>
#include <unordered_set>

#include "CellBased/Geometry.hpp"
#include "CellBased/Grid.hpp"
#include "CellBased/Interface.hpp"

namespace Zap::CellBased {

/*
 * Solve 2D Burgers equation of from d_t u + u * d_x u + u * d_y u = 0
 */
template <ExtendType extend_type>
class Solver {
  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, typename PassiveFloat>
  [[nodiscard]] constexpr auto
  cfl_factor(const UniformGrid<ActiveFloat, PassiveFloat>& grid) const noexcept -> ActiveFloat {
    return std::transform_reduce(
        std::cbegin(grid.cells()),
        std::cend(grid.cells()),
        ActiveFloat{0},
        [](const auto& a, const auto& b) { return std::max(a, b); },
        [](const Cell<ActiveFloat, PassiveFloat>& cell) -> ActiveFloat {
          assert(cell.is_cartesian() || cell.is_cut());
          if (cell.is_cartesian()) {
            return std::abs(cell.get_cartesian().value);
          } else {
            return std::max<ActiveFloat>(std::abs(cell.get_cut().left_value),
                                         std::abs(cell.get_cut().right_value));
          }
        });
  }

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, Point2D_c PointType>
  struct WaveProperties {
    ActiveFloat mass;
    Geometry::Polygon<PointType> polygon;
    ActiveFloat sign;
  };

  // -----------------------------------------------------------------------------------------------
  // TODO: Use this as reference
  // - normal direction:     https://github.com/clawpack/riemann/blob/master/src/rpn2_burgers.f90
  // - tangential direction: https://github.com/clawpack/riemann/blob/master/src/rpt2_burgers.f90
  template <typename ActiveFloat, typename PassiveFloat, Point2D_c PointType>
  [[nodiscard]] constexpr auto
  calculate_interface(const FullInterface<ActiveFloat, PointType>& interface,
                      ActiveFloat dt,
                      PassiveFloat scale_x,
                      PassiveFloat scale_y) const noexcept
      -> WaveProperties<ActiveFloat, PointType> {
    if ((interface.end - interface.begin).norm() < EPS<PassiveFloat>) { return {}; }

    // Vector tangential to interface
    const PointType tangent_vector = (interface.end - interface.begin).normalized();

    // Vector normal to cut
    const PointType normal_vector{tangent_vector.y, -tangent_vector.x};
    IGOR_ASSERT(std::abs(tangent_vector.dot(normal_vector)) <= EPS<PassiveFloat>,
                "tangent_vector {} is not orthogonal to normal_vector {}, tangent_vector^T * "
                "normal_vector = {}.",
                tangent_vector,
                normal_vector,
                tangent_vector.dot(normal_vector));

    // Scale tangential vector when operating in grid coordinates
    auto scale_if_grid_c = [=](const PointType& p) {
      if constexpr (is_GridCoord_v<PointType>) {
        return PointType{
            p.x * scale_x,
            p.y * scale_y,
        };
      } else {
        return p;
      }
    };
#ifdef ZAP_TANGENTIAL_CORRECTION
    const PointType scaled_tangent_vector = scale_if_grid_c(tangent_vector);
#endif  // ZAP_TANGENTIAL_CORRECTION
    const PointType scaled_normal_vector = scale_if_grid_c(normal_vector);

    const ActiveFloat u_mid = (interface.left_value + interface.right_value) / 2;

    // Matrix for rotated PDE in normal direction to cut
    // TODO: Why cos(interface_angle) * A + sin(interface_angle) * B and not -?
    // Equivalent to `cos(interface_angle) * A + sin(interface_angle) * B`, but do not actually
    // calculate the interface angle
    const ActiveFloat normal_mat =
        tangent_vector.y * u_mid +  //
        -sign(tangent_vector.x) * std::sqrt(1 - tangent_vector.y * tangent_vector.y) * u_mid;

    const ActiveFloat normal_wave_length = normal_mat * dt;

    // Eigen expansion of jump; wave strength
    const ActiveFloat normal_alpha = interface.right_value - interface.left_value;

#ifndef ZAP_TANGENTIAL_CORRECTION
    const WaveProperties<ActiveFloat, PointType> wave{
        .mass    = normal_alpha,
        .polygon = Geometry::Polygon<PointType>{{
            interface.begin,
            interface.begin + scaled_normal_vector * normal_wave_length,
            interface.end,
            interface.end + scaled_normal_vector * normal_wave_length,
        }},
        .sign    = sign(normal_mat),
    };
#else
    // Matrix for rotated PDE in tangential direction to cut
    // Equivalent to `-std::sin(interface_angle) * A + std::cos(interface_angle) * B`
    const ActiveFloat tangent_mat =
        sign(tangent_vector.x) * std::sqrt(1 - tangent_vector.y * tangent_vector.y) * u_mid +
        tangent_vector.y * u_mid;

    const ActiveFloat tangent_wave_length = tangent_mat * dt;
    const ActiveFloat tangent_beta        = normal_alpha;

    const WaveProperties<ActiveFloat, PointType> wave = {
        .mass    = tangent_beta,
        .polygon = Geometry::Polygon<PointType>{{
            interface.begin,
            interface.begin + scaled_normal_vector * normal_wave_length +
                scaled_tangent_vector * tangent_wave_length,
            interface.end,
            interface.end + scaled_normal_vector * normal_wave_length +
                scaled_tangent_vector * tangent_wave_length,
        }},
        .sign    = sign(normal_mat),
    };
#endif  // ZAP_NOT_TANGENTIAL_CORRECTION

    return wave;
  }

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, typename PassiveFloat, Point2D_c PointType>
  constexpr void update_cell_by_wave(Cell<ActiveFloat, PassiveFloat>& cell,
                                     const WaveProperties<ActiveFloat, PointType>& wave) {
    // TODO: Do something about periodic boundary conditions
    assert(cell.is_cartesian() || cell.is_cut());
    auto update_value = [&](ActiveFloat& value, const Geometry::Polygon<PointType>& cell_polygon) {
      const auto cell_area = cell_polygon.area();
      assert(cell_area > 0 || std::abs(cell_area) <= EPS<PassiveFloat>);
      if (std::abs(cell_area) <= EPS<PassiveFloat>) { return; }

      const auto intersect      = cell_polygon & wave.polygon;
      const auto intersect_area = intersect.area();
      IGOR_ASSERT(intersect_area >= 0 || std::abs(intersect_area) < EPS<PassiveFloat>,
                  "Expect the intersection area to be greater or equal to zero but is {}",
                  intersect_area);
      IGOR_ASSERT(intersect_area - cell_area <= 50 * EPS<PassiveFloat>,
                  "Expected area intersection of intersection to be smaller or equal to the area "
                  "of the cell, but intersection area is {} and cell area is {}, intersection area "
                  "is {} larger than cell area",
                  intersect_area,
                  cell_area,
                  static_cast<ActiveFloat>(intersect_area - cell_area));

      value += -wave.sign * (intersect_area / cell_area) * wave.mass;
    };

    constexpr CoordType coord_type = PointType2CoordType<PointType>;
    if (cell.is_cartesian()) {
      update_value(cell.get_cartesian().value, cell.template get_cartesian_polygon<coord_type>());
    } else {
      // Left subcell
      update_value(cell.get_cut().left_value, cell.template get_cut_left_polygon<coord_type>());
      // Right subcell
      update_value(cell.get_cut().right_value, cell.template get_cut_right_polygon<coord_type>());
    }
  }

  // -----------------------------------------------------------------------------------------------
#ifndef ZAP_TANGENTIAL_CORRECTION
  template <typename ActiveFloat, typename PassiveFloat>
  constexpr void
  apply_side_interfaces_uncorrected(const UniformGrid<ActiveFloat, PassiveFloat>& curr_grid,
                                    Cell<ActiveFloat, PassiveFloat>& next_cell,
                                    const Cell<ActiveFloat, PassiveFloat>& curr_cell,
                                    size_t idx,
                                    Side side,
                                    ActiveFloat dt) {
    if (curr_grid.is_cell(idx)) {
      const auto& other_cell = curr_grid[idx];
      const auto interfaces =
          get_shared_interfaces<ActiveFloat, PassiveFloat, GridCoord<ActiveFloat>>(
              curr_cell, other_cell, side);
      for (const auto& interface : interfaces) {
        const auto wave =
            calculate_interface(interface, dt, curr_grid.scale_x(), curr_grid.scale_y());
        update_cell_by_wave(next_cell, wave);
      }
    } else if (idx == SAME_VALUE_INDEX || idx == ZERO_FLUX_INDEX) {
      // TODO: Noop? Same?
    } else if (idx == NULL_INDEX) {
      Igor::Debug("cell = {}", curr_cell);
      Igor::Panic("cell has NULL index.");
    } else {
      Igor::Debug("cell = {}", curr_cell);
      Igor::Panic("cell has unknown index with value {}.", idx);
    }
  };
#endif  // ZAP_TANGENTIAL_CORRECTION

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, typename PassiveFloat, Point2D_c PointType>
  [[nodiscard]] constexpr auto
  get_internal_interface(const Cell<ActiveFloat, PassiveFloat>& cell) const noexcept
      -> FullInterface<ActiveFloat, PointType> {
    assert(cell.is_cut());

    constexpr CoordType coord_type = PointType2CoordType<PointType>;

    return FullInterface<ActiveFloat, PointType>{
        .left_value  = cell.get_cut().left_value,
        .right_value = cell.get_cut().right_value,
        .begin       = cell.template cut_entry<coord_type>(),
        .end         = cell.template cut_exit<coord_type>(),
    };
  }

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, typename PassiveFloat>
  [[nodiscard]] constexpr auto
  move_wave_front(const UniformGrid<ActiveFloat, PassiveFloat>& curr_grid,
                  [[maybe_unused]] const UniformGrid<ActiveFloat, PassiveFloat>& next_grid,
                  ActiveFloat dt) const noexcept
      -> std::optional<std::vector<GridCoord<ActiveFloat>>> {
    using PointType = GridCoord<ActiveFloat>;

    // Move old cuts according to strongest wave
    std::vector<std::pair<PointType, PointType>> new_shock_points;
    for (size_t cell_idx : curr_grid.cut_cell_idxs()) {
      const auto& curr_cell = curr_grid[cell_idx];
      assert(curr_cell.is_cut());

      const auto interface =
          get_internal_interface<ActiveFloat, PassiveFloat, PointType>(curr_cell);

      const auto tangent_vector = (interface.end - interface.begin).normalized();
      assert(std::abs(tangent_vector.norm() - 1) <= 1e-8);

      PointType normal_vector{tangent_vector.y, -tangent_vector.x};
      IGOR_ASSERT(std::abs(tangent_vector.dot(normal_vector)) <= EPS<PassiveFloat>,
                  "tangent_vector {} is not orthogonal to normal_vector {}, tangent_vector^T * "
                  "normal_vector = {}",
                  tangent_vector,
                  normal_vector,
                  tangent_vector.dot(normal_vector));
      // Scale normal vector when operating in grid coordinates
      if constexpr (is_GridCoord_v<PointType>) {
        normal_vector.x *= curr_grid.scale_x();
        normal_vector.y *= curr_grid.scale_y();
      }

      const ActiveFloat u_mid = (interface.left_value + interface.right_value) / 2;

      // Equivalent to `cos(interface_angle) * A + sin(interface_angle) * B`, but do not actually
      // calculate the interface angle
      const ActiveFloat eta_mat =
          tangent_vector.y * u_mid +
          -sign(tangent_vector.x) * std::sqrt(1 - tangent_vector.y * tangent_vector.y) * u_mid;
      const ActiveFloat wave_length = eta_mat * dt;

      new_shock_points.emplace_back(interface.begin + normal_vector * wave_length,
                                    interface.end + normal_vector * wave_length);
    }

    std::vector<PointType> avg_new_shock_points(new_shock_points.size() + 1);
    avg_new_shock_points[0] = new_shock_points[0].first;
    for (size_t i = 1; i < avg_new_shock_points.size() - 1; ++i) {
      avg_new_shock_points[i] = (new_shock_points[i - 1].second + new_shock_points[i].first) / 2;
    }
    avg_new_shock_points[new_shock_points.size()] = new_shock_points.back().second;

    if ((avg_new_shock_points[0] - avg_new_shock_points[1]).norm() < 1e-2) {
      avg_new_shock_points.erase(std::next(std::begin(avg_new_shock_points)));
    }
    if ((avg_new_shock_points[avg_new_shock_points.size() - 1] -
         avg_new_shock_points[avg_new_shock_points.size() - 2])
            .norm() < 1e-2) {
      avg_new_shock_points.erase(std::prev(std::end(avg_new_shock_points), 2));
    }

    return avg_new_shock_points;
  }

 public:
  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, typename PassiveFloat, typename GridWriter, typename TimeWriter>
  [[nodiscard]] auto solve(UniformGrid<ActiveFloat, PassiveFloat> grid,
                           PassiveFloat tend,
                           GridWriter& grid_writer,
                           TimeWriter& time_writer,
                           PassiveFloat CFL_safety_factor = 0.5) noexcept
      -> std::optional<UniformGrid<ActiveFloat, PassiveFloat>> {
    if (!(CFL_safety_factor > 0 && CFL_safety_factor <= 1)) {
      Igor::Warn("CFL_safety_factor must be in (0, 1], is {}", CFL_safety_factor);
      return std::nullopt;
    }
    auto& curr_grid = grid;
    auto next_grid  = grid;

    if (!grid_writer.write_data(curr_grid) || !time_writer.write_data(PassiveFloat{0})) {
      return std::nullopt;
    }

    for (ActiveFloat t = 0.0; t < tend;) {
      const ActiveFloat CFL_factor = cfl_factor(curr_grid);

      if (std::isnan(CFL_factor) || std::isinf(CFL_factor)) {
        Igor::Warn("CFL_factor is invalid at time t={}: CFL_factor = {}", t, CFL_factor);
        return std::nullopt;
      }
      const ActiveFloat dt =
          std::min(CFL_safety_factor * curr_grid.min_delta() / CFL_factor, tend - t);

      next_grid = curr_grid;

#ifndef ZAP_STATIC_CUT
      if (!curr_grid.cut_cell_idxs().empty()) {
        next_grid.merge_cut_cells();

        auto avg_new_shock_points = move_wave_front(curr_grid, next_grid, dt);
        if (!avg_new_shock_points.has_value()) { return std::nullopt; }
        if (!next_grid.template cut_piecewise_linear<extend_type>(
                std::move(*avg_new_shock_points))) {
          Igor::Warn("Could not cut on new shock curve.");
          return std::nullopt;
        }
      }

      // Re-calculate value for newly cut cells
      for (size_t new_cut_idx : next_grid.cut_cell_idxs()) {
        auto& next_cell       = next_grid[new_cut_idx];
        const auto& curr_cell = curr_grid[new_cut_idx];
        assert(next_cell.is_cut());
        assert(curr_cell.is_cartesian() || curr_cell.is_cut());

        if (curr_cell.is_cut()) {
          const auto curr_cell_left_polygon  = curr_cell.template get_cut_left_polygon<GRID_C>();
          const auto curr_cell_right_polygon = curr_cell.template get_cut_right_polygon<GRID_C>();

          // Left subcell
          {
            const auto next_subcell_polygon = next_cell.template get_cut_left_polygon<GRID_C>();
            assert(next_subcell_polygon.area() > 0 ||
                   std::abs(next_subcell_polygon.area()) <= EPS<PassiveFloat>);
            if (std::abs(next_subcell_polygon.area()) > EPS<PassiveFloat>) {
              const auto left_intersect_area =
                  (next_subcell_polygon & curr_cell_left_polygon).area();
              const auto right_intersect_area =
                  (next_subcell_polygon & curr_cell_right_polygon).area();

              IGOR_ASSERT(std::abs(left_intersect_area + right_intersect_area -
                                   next_subcell_polygon.area()) < 1e-6,
                          "Expect the sum of the intersected areas of the subcells to be equal to "
                          "the area of the entire subcell, but {} + {} = {} != {}",
                          left_intersect_area,
                          right_intersect_area,
                          static_cast<ActiveFloat>(left_intersect_area + right_intersect_area),
                          next_subcell_polygon.area());

              next_cell.get_cut().left_value =
                  (curr_cell.get_cut().left_value * left_intersect_area +
                   curr_cell.get_cut().right_value * right_intersect_area) /
                  next_subcell_polygon.area();
            }
          }

          // Right subcell
          {
            const auto next_subcell_polygon = next_cell.template get_cut_right_polygon<GRID_C>();
            assert(next_subcell_polygon.area() > 0 ||
                   std::abs(next_subcell_polygon.area()) <= EPS<PassiveFloat>);
            if (std::abs(next_subcell_polygon.area()) > EPS<PassiveFloat>) {
              const auto left_intersect_area =
                  (next_subcell_polygon & curr_cell_left_polygon).area();
              const auto right_intersect_area =
                  (next_subcell_polygon & curr_cell_right_polygon).area();
              assert(std::abs(left_intersect_area + right_intersect_area -
                              next_subcell_polygon.area()) < 1e-6);

              next_cell.get_cut().right_value =
                  (curr_cell.get_cut().left_value * left_intersect_area +
                   curr_cell.get_cut().right_value * right_intersect_area) /
                  next_subcell_polygon.area();
            }
          }
        }
      }

      // TODO: Consider purging cuts that are no longer at the shock front using the
      //       RK-Jump-Condition
#endif  // ZAP_STATIC_CUT

      // TODO: What happens when a subcell has area 0?
      //   -> cannot "uncut" the cell because that would loose shock position information

#pragma omp parallel for
      for (size_t cell_idx = 0; cell_idx < curr_grid.size(); ++cell_idx) {
        const auto& curr_cell = curr_grid[cell_idx];
        auto& next_cell       = next_grid[cell_idx];

        // = Cut cell ==============================================================================
        if (next_cell.is_cut()) {
#ifndef ZAP_TANGENTIAL_CORRECTION
          // Left interface
          apply_side_interfaces_uncorrected(
              curr_grid, next_cell, curr_cell, curr_cell.left_idx, LEFT, dt);
          // Right interface
          apply_side_interfaces_uncorrected(
              curr_grid, next_cell, curr_cell, curr_cell.right_idx, RIGHT, dt);
          // Bottom interface
          apply_side_interfaces_uncorrected(
              curr_grid, next_cell, curr_cell, curr_cell.bottom_idx, BOTTOM, dt);
          // Top interface
          apply_side_interfaces_uncorrected(
              curr_grid, next_cell, curr_cell, curr_cell.top_idx, TOP, dt);
#else
          continue;
#endif  // ZAP_TANGENTIAL_CORRECTION
        }
        // = Cut cell ==============================================================================
        // = Cartesian cell ========================================================================
        else {
          // - Handle left interface ---------------------------------------------------------------
          if (curr_grid.is_cell(curr_cell.left_idx)) {
            const auto& other_cell = curr_grid[curr_cell.left_idx];
            const auto interfaces =
                get_shared_interfaces<ActiveFloat, PassiveFloat, SimCoord<ActiveFloat>>(
                    curr_cell, other_cell, LEFT);

            for (const auto& interface : interfaces) {
              const ActiveFloat u_mid  = (interface.left_value + interface.right_value) / 2;
              const ActiveFloat A_plus = std::max(u_mid, ActiveFloat{0});

              next_cell.get_cartesian().value -=
                  A_plus *
                  ((interface.end - interface.begin).norm() / curr_cell.template dy<SIM_C>()) *
                  (dt / curr_cell.template dx<SIM_C>()) *
                  (interface.right_value - interface.left_value);
            }
          } else if (curr_cell.left_idx == SAME_VALUE_INDEX ||
                     curr_cell.left_idx == ZERO_FLUX_INDEX) {
            // No-op
          } else {
            Igor::Debug("cell = {}", curr_cell);
            Igor::Panic("Invalid left-index.");
          }
          // - Handle left interface ---------------------------------------------------------------

          // - Handle right interface --------------------------------------------------------------
          if (curr_grid.is_cell(curr_cell.right_idx)) {
            const auto& other_cell = curr_grid[curr_cell.right_idx];
            const auto interfaces =
                get_shared_interfaces<ActiveFloat, PassiveFloat, SimCoord<ActiveFloat>>(
                    curr_cell, other_cell, RIGHT);

            for (const auto& interface : interfaces) {
              const ActiveFloat u_mid   = (interface.left_value + interface.right_value) / 2;
              const ActiveFloat A_minus = std::min(u_mid, ActiveFloat{0.0});

              // TODO: This might need to be `-=`, but the flux is always zero so it does not show
              // up
              next_cell.get_cartesian().value +=
                  A_minus *
                  ((interface.end - interface.begin).norm() / curr_cell.template dy<SIM_C>()) *
                  (dt / curr_cell.template dx<SIM_C>()) *
                  (interface.right_value - interface.left_value);
            }
          } else if (curr_cell.right_idx == SAME_VALUE_INDEX ||
                     curr_cell.right_idx == ZERO_FLUX_INDEX) {
            // No-op
          } else {
            Igor::Debug("cell = {}", curr_cell);
            Igor::Panic("Invalid right-index.");
          }
          // - Handle right interface --------------------------------------------------------------

          // - Handle bottom interface -------------------------------------------------------------
          if (curr_grid.is_cell(curr_cell.bottom_idx)) {
            const auto& other_cell = curr_grid[curr_cell.bottom_idx];
            const auto interfaces =
                get_shared_interfaces<ActiveFloat, PassiveFloat, SimCoord<ActiveFloat>>(
                    curr_cell, other_cell, BOTTOM);

            for (const auto& interface : interfaces) {
              const ActiveFloat u_mid  = (interface.left_value + interface.right_value) / 2;
              const ActiveFloat A_plus = std::max(u_mid, ActiveFloat{0.0});

              next_cell.get_cartesian().value -=
                  A_plus *
                  ((interface.end - interface.begin).norm() / curr_cell.template dx<SIM_C>()) *
                  (dt / curr_cell.template dy<SIM_C>()) *
                  (interface.right_value - interface.left_value);
            }
          } else if (curr_cell.bottom_idx == SAME_VALUE_INDEX ||
                     curr_cell.bottom_idx == ZERO_FLUX_INDEX) {
            // No-op
          } else {
            Igor::Debug("cell = {}", curr_cell);
            Igor::Panic("Invalid bottom-index.");
          }
          // - Handle bottom interface -------------------------------------------------------------

          // - Handle top interface ----------------------------------------------------------------
          if (curr_grid.is_cell(curr_cell.top_idx)) {
            const auto& other_cell = curr_grid[curr_cell.top_idx];
            const auto interfaces =
                get_shared_interfaces<ActiveFloat, PassiveFloat, SimCoord<ActiveFloat>>(
                    curr_cell, other_cell, TOP);

            for (const auto& interface : interfaces) {
              const ActiveFloat u_mid   = (interface.left_value + interface.right_value) / 2;
              const ActiveFloat A_minus = std::min(u_mid, 0.0);

              // TODO: This might need to be `-=`, but the flux is always zero so it does not show
              next_cell.get_cartesian().value +=
                  A_minus *
                  ((interface.end - interface.begin).norm() / curr_cell.template dx<SIM_C>()) *
                  (dt / curr_cell.template dy<SIM_C>()) *
                  (interface.right_value - interface.left_value);
            }
          } else if (curr_cell.top_idx == SAME_VALUE_INDEX ||
                     curr_cell.top_idx == ZERO_FLUX_INDEX) {
            // No-op
          } else {
            Igor::Debug("cell = {}", curr_cell);
            Igor::Panic("Invalid top-index.");
          }
          // - Handle top interface ----------------------------------------------------------------
        }
        // = Cartesian cell ========================================================================
      }

#ifdef ZAP_TANGENTIAL_CORRECTION
      auto pair_hasher = [](const std::pair<size_t, size_t>& p) { return p.first ^ p.second; };
      auto pair_equal  = [](const std::pair<size_t, size_t>& p1,
                           const std::pair<size_t, size_t>& p2) {
        return (p1.first == p2.first && p1.second == p2.second) ||
               (p1.first == p2.second && p1.second == p2.first);
      };
      std::unordered_set<std::pair<size_t, size_t>, decltype(pair_hasher), decltype(pair_equal)>
          calculated_interfaces;
      for (size_t cell_idx : next_grid.cut_cell_idxs()) {
        const auto& curr_cell      = curr_grid[cell_idx];
        auto& next_cell            = next_grid[cell_idx];
        const auto cell_neighbours = curr_grid.get_existing_neighbours(cell_idx);

        // - Left outer interface ------------------------------------------------------------------
        if (curr_grid.is_cell(curr_cell.left_idx) &&
            !calculated_interfaces.contains({cell_idx, curr_cell.left_idx})) {
          const auto& other_cell = curr_grid[curr_cell.left_idx];
          const auto interfaces =
              get_shared_interfaces<ActiveFloat, PassiveFloat, GridCoord<ActiveFloat>>(
                  curr_cell, other_cell, LEFT);

          for (const auto& interface : interfaces) {
            const auto wave =
                calculate_interface(interface, dt, curr_grid.scale_x(), curr_grid.scale_y());
            for (size_t neighbour_idx : cell_neighbours) {
              update_cell_by_wave(next_grid[neighbour_idx], wave);
            }
            update_cell_by_wave(next_cell, wave);
          }
          calculated_interfaces.emplace(cell_idx, curr_cell.left_idx);
        }

        // - Right outer interface -----------------------------------------------------------------
        if (curr_grid.is_cell(curr_cell.right_idx) &&
            !calculated_interfaces.contains({cell_idx, curr_cell.right_idx})) {
          const auto& other_cell = curr_grid[curr_cell.right_idx];
          const auto interfaces =
              get_shared_interfaces<ActiveFloat, PassiveFloat, GridCoord<ActiveFloat>>(
                  curr_cell, other_cell, RIGHT);

          for (const auto& interface : interfaces) {
            const auto wave =
                calculate_interface(interface, dt, curr_grid.scale_x(), curr_grid.scale_y());
            for (size_t neighbour_idx : cell_neighbours) {
              update_cell_by_wave(next_grid[neighbour_idx], wave);
            }
            update_cell_by_wave(next_cell, wave);
          }
          calculated_interfaces.emplace(cell_idx, curr_cell.right_idx);
        }

        // - Bottom outer interface ----------------------------------------------------------------
        if (curr_grid.is_cell(curr_cell.bottom_idx) &&
            !calculated_interfaces.contains({cell_idx, curr_cell.bottom_idx})) {
          const auto& other_cell = curr_grid[curr_cell.bottom_idx];
          const auto interfaces =
              get_shared_interfaces<ActiveFloat, PassiveFloat, GridCoord<ActiveFloat>>(
                  curr_cell, other_cell, BOTTOM);

          for (const auto& interface : interfaces) {
            const auto wave =
                calculate_interface(interface, dt, curr_grid.scale_x(), curr_grid.scale_y());
            for (size_t neighbour_idx : cell_neighbours) {
              update_cell_by_wave(next_grid[neighbour_idx], wave);
            }
            update_cell_by_wave(next_cell, wave);
          }
          calculated_interfaces.emplace(cell_idx, curr_cell.bottom_idx);
        }

        // - Top outer interface -------------------------------------------------------------------
        if (curr_grid.is_cell(curr_cell.top_idx) &&
            !calculated_interfaces.contains({cell_idx, curr_cell.top_idx})) {
          const auto& other_cell = curr_grid[curr_cell.top_idx];
          const auto interfaces =
              get_shared_interfaces<ActiveFloat, PassiveFloat, GridCoord<ActiveFloat>>(
                  curr_cell, other_cell, TOP);

          for (const auto& interface : interfaces) {
            const auto wave =
                calculate_interface(interface, dt, curr_grid.scale_x(), curr_grid.scale_y());
            for (size_t neighbour_idx : cell_neighbours) {
              update_cell_by_wave(next_grid[neighbour_idx], wave);
            }
            update_cell_by_wave(next_cell, wave);
          }
          calculated_interfaces.emplace(cell_idx, curr_cell.top_idx);
        }
      }
#endif  // ZAP_TANGENTIAL_CORRECTION

      // - Handle internal interface ---------------------------------------------------------------
      for (size_t cell_idx : curr_grid.cut_cell_idxs()) {
        const auto& curr_cell = curr_grid[cell_idx];
        assert(curr_cell.is_cut());
        auto& next_cell = next_grid[cell_idx];
        assert(next_cell.is_cartesian() || next_cell.is_cut());

        const auto internal_interface =
            get_internal_interface<ActiveFloat, PassiveFloat, GridCoord<ActiveFloat>>(curr_cell);
        const auto internal_wave =
            calculate_interface<ActiveFloat, PassiveFloat, GridCoord<ActiveFloat>>(
                internal_interface, dt, curr_grid.scale_x(), curr_grid.scale_y());

        const auto cell_neighbours = curr_grid.get_existing_neighbours(cell_idx);

        // Neighbouring cells
        for (size_t neighbour_idx : cell_neighbours) {
          update_cell_by_wave(next_grid[neighbour_idx], internal_wave);
        }
        // This cell
        update_cell_by_wave(next_cell, internal_wave);
      }
      // - Handle internal interface ---------------------------------------------------------------

      // Update time
      t += dt;

      curr_grid = next_grid;

      if (!grid_writer.write_data(curr_grid) || !time_writer.write_data(ad::value(t))) {
        return std::nullopt;
      }
    }

    return curr_grid;
  }
};

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_SOLVER_HPP_
