#ifndef ZAP_CELL_BASED_SOLVER_HPP_
#define ZAP_CELL_BASED_SOLVER_HPP_

#include <numeric>
#include <unordered_set>

#include "CellBased/Geometry.hpp"
#include "CellBased/Grid.hpp"
#include "CellBased/Interface.hpp"
#include "CellBased/Wave.hpp"

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
  template <Side side, Point2D_c PointType, typename ActiveFloat, typename Wave>
  constexpr void update_value_by_overlap(ActiveFloat& value,
                                         const Geometry::Polygon<PointType>& cell_polygon,
                                         const Geometry::Polygon<PointType>& wave_polygon,
                                         const Wave& wave) {
    // clang-format off
    static_assert(
      // x-aligned wave must come from left or right
      (std::is_same_v<Wave, AxisAlignedWave<ActiveFloat, PointType, X>> && (side == LEFT || side == RIGHT)) ||
      // y-aligned wave must come from bottom or top
      (std::is_same_v<Wave, AxisAlignedWave<ActiveFloat, PointType, Y>> && (side == BOTTOM || side == TOP)) ||
      // free wave cannot come from any side
      (std::is_same_v<Wave, FreeWave<ActiveFloat, PointType>> && count_sides(side) == 0),
       "Invalid combination of wave-type and side parameter.");
    // clang-format on

    const auto cell_area = cell_polygon.area();
    assert(cell_area > 0 || std::abs(cell_area) <= EPS<ActiveFloat>);
    if (std::abs(cell_area) <= EPS<ActiveFloat>) { return; }

    const auto intersect_area = (cell_polygon & wave_polygon).area();
    IGOR_ASSERT(intersect_area >= 0 || std::abs(intersect_area) < EPS<ActiveFloat>,
                "Expect the intersection area to be greater or equal to zero but is {}",
                intersect_area);
    IGOR_ASSERT(intersect_area - cell_area <= 50 * EPS<ActiveFloat>,
                "Expected area of intersection to be smaller or equal to the area of the cell, "
                "but intersection area is {} and cell area is {}, intersection area is {} "
                "larger than cell area.",
                intersect_area,
                cell_area,
                static_cast<ActiveFloat>(intersect_area - cell_area));

    value -= (intersect_area / cell_area) * wave.first_order_update;
    // TODO: Maybe dont do that.
    // if constexpr ((side == LEFT || side == BOTTOM)) {
    //   value -= (intersect_area / cell_area) * wave.second_order_update;
    // } else if constexpr ((side == RIGHT) || (side == TOP)) {
    //   value -= (intersect_area / cell_area) * (-wave.second_order_update);
    // }
  }

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, typename PassiveFloat, Point2D_c PointType>
  constexpr void update_cell_by_free_wave(Cell<ActiveFloat, PassiveFloat>& cell,
                                          const FreeWave<ActiveFloat, PointType>& wave,
                                          const ActiveFloat& dt) {
    constexpr CoordType coord_type                  = PointType2CoordType<PointType>;
    const Geometry::Polygon<PointType> wave_polygon = calc_wave_polygon(wave, dt);

    if (cell.is_cartesian()) {
      update_value_by_overlap<static_cast<Side>(0)>(
          cell.get_cartesian().value,
          cell.template get_cartesian_polygon<coord_type>(),
          wave_polygon,
          wave);
    } else {
      // Left subcell
      update_value_by_overlap<static_cast<Side>(0)>(
          cell.get_cut().left_value,
          cell.template get_cut_left_polygon<coord_type>(),
          wave_polygon,
          wave);
      // Right subcell
      update_value_by_overlap<static_cast<Side>(0)>(
          cell.get_cut().right_value,
          cell.template get_cut_right_polygon<coord_type>(),
          wave_polygon,
          wave);
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat,
            typename PassiveFloat,
            Point2D_c PointType,
            Orientation orientation,
            Side side>
  constexpr void
  update_cell_by_axis_aligned_wave(Cell<ActiveFloat, PassiveFloat>& cell,
                                   const AxisAlignedWave<ActiveFloat, PointType, orientation>& wave,
                                   const ActiveFloat& dt) {
    // clang-format off
    static_assert(
      // x-aligned wave must come from left or right
      (orientation == X && (side == LEFT || side == RIGHT)) ||
      // y-aligned wave must come from bottom or top
      (orientation == Y && (side == BOTTOM || side == TOP)),
       "Invalid combination of orientation and side parameter.");
    // clang-format on

    const auto ds = orientation == X ? cell.template dx<SIM_C>() : cell.template dy<SIM_C>();
    if (cell.is_cartesian()) {
      if constexpr ((side == LEFT || side == BOTTOM)) {
        if (wave.is_right_going) {
          cell.get_cartesian().value -= (dt / ds) * wave.first_order_update;
        }
        cell.get_cartesian().value -= (dt / ds) * wave.second_order_update;
      } else {
        if (!wave.is_right_going) {
          cell.get_cartesian().value -= (dt / ds) * wave.first_order_update;
        }
        cell.get_cartesian().value -= (dt / ds) * (-wave.second_order_update);
      }
    } else {
      const Geometry::Polygon<PointType> wave_polygon =
          calc_wave_polygon(wave, cell.template dx<SIM_C>(), cell.template dy<SIM_C>(), dt);

      constexpr CoordType coord_type = PointType2CoordType<PointType>;
      // Left subcell
      update_value_by_overlap<side>(cell.get_cut().left_value,
                                    cell.template get_cut_left_polygon<coord_type>(),
                                    wave_polygon,
                                    wave);
      // Right subcell
      update_value_by_overlap<side>(cell.get_cut().right_value,
                                    cell.template get_cut_right_polygon<coord_type>(),
                                    wave_polygon,
                                    wave);
    }
  }

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
    using PointType = GridCoord<ActiveFloat>;

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

      // = Handle outer interfaces: LEFT, RIGHT, BOTTOM, and TOP ===================================
#pragma omp parallel for
      for (size_t cell_idx = 0; cell_idx < curr_grid.size(); ++cell_idx) {
        const auto& curr_cell = curr_grid[cell_idx];
        auto& next_cell       = next_grid[cell_idx];

        // - Handle left interface ---------------------------------------------------------------
        if (curr_grid.is_cell(curr_cell.left_idx)) {
          const auto& other_cell = curr_grid[curr_cell.left_idx];
          const auto interfaces  = get_shared_interfaces<ActiveFloat, PassiveFloat, PointType>(
              curr_cell, other_cell, LEFT);

          for (const auto& interface : interfaces) {
            const auto wave = normal_wave<X>(
                interface, curr_cell.template dx<SIM_C>(), curr_cell.template dy<SIM_C>(), dt);
            update_cell_by_axis_aligned_wave<ActiveFloat, PassiveFloat, PointType, X, LEFT>(
                next_cell, wave, dt);
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
          const auto interfaces  = get_shared_interfaces<ActiveFloat, PassiveFloat, PointType>(
              curr_cell, other_cell, RIGHT);

          for (const auto& interface : interfaces) {
            const auto wave = normal_wave<X>(
                interface, curr_cell.template dx<SIM_C>(), curr_cell.template dy<SIM_C>(), dt);
            update_cell_by_axis_aligned_wave<ActiveFloat, PassiveFloat, PointType, X, RIGHT>(
                next_cell, wave, dt);
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
          const auto interfaces  = get_shared_interfaces<ActiveFloat, PassiveFloat, PointType>(
              curr_cell, other_cell, BOTTOM);

          for (const auto& interface : interfaces) {
            const auto wave = normal_wave<Y>(
                interface, curr_cell.template dx<SIM_C>(), curr_cell.template dy<SIM_C>(), dt);
            update_cell_by_axis_aligned_wave<ActiveFloat, PassiveFloat, PointType, Y, BOTTOM>(
                next_cell, wave, dt);
          }
        } else if (curr_cell.bottom_idx == SAME_VALUE_INDEX ||
                   curr_cell.bottom_idx == ZERO_FLUX_INDEX) {
          // No-op
        } else {
          Igor::Debug("cell = {}", curr_cell);
          Igor::Panic("Invalid bottom-index.");
        }
        // - Handle bottom interface -------------------------------------------------------------

        // - Handle top interface ---------------------------------------------------------------
        if (curr_grid.is_cell(curr_cell.top_idx)) {
          const auto& other_cell = curr_grid[curr_cell.top_idx];
          const auto interfaces  = get_shared_interfaces<ActiveFloat, PassiveFloat, PointType>(
              curr_cell, other_cell, TOP);

          for (const auto& interface : interfaces) {
            const auto wave = normal_wave<Y>(
                interface, curr_cell.template dx<SIM_C>(), curr_cell.template dy<SIM_C>(), dt);
            update_cell_by_axis_aligned_wave<ActiveFloat, PassiveFloat, PointType, Y, TOP>(
                next_cell, wave, dt);
          }
        } else if (curr_cell.top_idx == SAME_VALUE_INDEX || curr_cell.top_idx == ZERO_FLUX_INDEX) {
          // No-op
        } else {
          Igor::Debug("cell = {}", curr_cell);
          Igor::Panic("Invalid top-index.");
        }
        // - Handle top interface ---------------------------------------------------------------
      }
      // = Handle outer interfaces: LEFT, RIGHT, BOTTOM, and TOP ===================================

      // = Handle internal interface ===============================================================
      for (size_t cell_idx : curr_grid.cut_cell_idxs()) {
        const auto& curr_cell = curr_grid[cell_idx];
        assert(curr_cell.is_cut());
        auto& next_cell = next_grid[cell_idx];
        assert(next_cell.is_cartesian() || next_cell.is_cut());

        const auto internal_interface =
            get_internal_interface<ActiveFloat, PassiveFloat, PointType>(curr_cell);
        const auto internal_wave = normal_wave<FREE>(
            internal_interface, curr_cell.template dx<SIM_C>(), curr_cell.template dy<SIM_C>(), dt);

        const auto cell_neighbours = curr_grid.get_existing_neighbours(cell_idx);
        // Neighbouring cells
        for (size_t neighbour_idx : cell_neighbours) {
          update_cell_by_free_wave(next_grid[neighbour_idx], internal_wave, dt);
        }
        // This cell
        update_cell_by_free_wave(next_cell, internal_wave, dt);
      }
      // = Handle internal interface ===============================================================

      // Update time
      t += dt;

      curr_grid = next_grid;

      if (!grid_writer.write_data(curr_grid) || !time_writer.write_data(ad::value(t))) {
        return std::nullopt;
      }
      break;
    }

    return curr_grid;
  }
};

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_SOLVER_HPP_
