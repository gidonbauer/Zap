#ifndef ZAP_CELL_BASED_SOLVER_HPP_
#define ZAP_CELL_BASED_SOLVER_HPP_

#include <numeric>
#include <unordered_set>

#ifndef ZAP_NO_PARALLEL
#include <omp.h>
#endif  // ZAP_NO_PARALLEL

#include "CellBased/Geometry.hpp"
#include "CellBased/Grid.hpp"
#include "CellBased/Interface.hpp"
#include "CellBased/Wave.hpp"

namespace Zap::CellBased {

// TODO: Fix periodic boundary conditions
// TODO: Fix bug in parallelization -> order of cell that are updated matters

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +              Solve 2D Burgers equation of from d_t u + u * d_x u + u * d_y u = 0              +
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
                                         const Wave& wave) const noexcept {
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
                                          const ActiveFloat& dt) const noexcept {
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
                                   const ActiveFloat& dt) const noexcept {
    // clang-format off
    static_assert(
      // x-aligned wave must come from left or right
      (orientation == X && (side == LEFT || side == RIGHT)) ||
      // y-aligned wave must come from bottom or top
      (orientation == Y && (side == BOTTOM || side == TOP)),
       "Invalid combination of orientation and side parameter.");
    // clang-format on

    IGOR_ASSERT(approx_eq(wave.tangent_speed, ActiveFloat{0}),
                "This function may only be called if the tangential speed is zero, but the "
                "tangential speed is {}",
                wave.tangent_speed);

    if (cell.is_cartesian()) {
      const auto ds = orientation == X ? cell.template dx<SIM_C>() : cell.template dy<SIM_C>();
      if constexpr ((side == LEFT || side == BOTTOM)) {
        if (wave.normal_speed >= 0) {
          cell.get_cartesian().value -=
              (dt / ds) * std::abs(wave.normal_speed) * wave.first_order_update;
        }
        cell.get_cartesian().value -= (dt / ds) * wave.second_order_update;
      } else {
        if (!(wave.normal_speed >= 0)) {
          cell.get_cartesian().value -=
              (dt / ds) * std::abs(wave.normal_speed) * wave.first_order_update;
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
  template <Point2D_c PointType, Side side, typename ActiveFloat, typename PassiveFloat>
  requires(side == BOTTOM || side == RIGHT || side == TOP || side == LEFT)
  void handle_side_iterface(Cell<ActiveFloat, PassiveFloat>& next_cell,
                            const Cell<ActiveFloat, PassiveFloat>& curr_cell,
                            UniformGrid<ActiveFloat, PassiveFloat>& next_grid,
                            const UniformGrid<ActiveFloat, PassiveFloat>& curr_grid,
                            const ActiveFloat& dt) const noexcept {
    size_t idx = NULL_INDEX;
    switch (side) {
      case BOTTOM: idx = curr_cell.bottom_idx; break;
      case RIGHT:  idx = curr_cell.right_idx; break;
      case TOP:    idx = curr_cell.top_idx; break;
      case LEFT:   idx = curr_cell.left_idx; break;
    }

    constexpr auto orientation = [] -> Orientation {
      switch (side) {
        case RIGHT:
        case LEFT:   return X;
        case BOTTOM:
        case TOP:    return Y;
      }
      std::unreachable();
    }();

    // - Handle side interface ---------------------------------------------------------------
    if (curr_grid.is_cell(idx)) {
      const auto& other_cell = curr_grid[idx];
      const auto interfaces =
          get_shared_interfaces<ActiveFloat, PassiveFloat, PointType>(curr_cell, other_cell, side);

      for (const auto& interface : interfaces) {
        const auto wave = normal_wave<orientation>(
            interface, curr_cell.template dx<SIM_C>(), curr_cell.template dy<SIM_C>(), dt);

        if (!wave.has_value() || ((side == LEFT || side == BOTTOM) && wave->normal_speed <= 0) ||
            ((side == RIGHT || side == TOP) && wave->normal_speed >= 0)) {
          continue;
        }

        // = Update this cell =========================================
        if (next_cell.is_cartesian()) {
          const auto ds =
              orientation == X ? curr_cell.template dy<SIM_C>() : curr_cell.template dx<SIM_C>();
          const ActiveFloat overlap_area =
              std::abs(wave->normal_speed) * dt * (ds - 0.5 * std::abs(wave->tangent_speed) * dt);
          const ActiveFloat cell_area =
              curr_cell.template dx<SIM_C>() * curr_cell.template dy<SIM_C>();

          next_cell.get_cartesian().value -= (overlap_area / cell_area) * wave->first_order_update;
        } else {
          const Geometry::Polygon<PointType> wave_polygon = calc_wave_polygon(
              *wave, curr_cell.template dx<SIM_C>(), curr_cell.template dy<SIM_C>(), dt);

          constexpr CoordType coord_type = PointType2CoordType<PointType>;
          // Left subcell
          update_value_by_overlap<side>(next_cell.get_cut().left_value,
                                        next_cell.template get_cut_left_polygon<coord_type>(),
                                        wave_polygon,
                                        *wave);
          // Right subcell
          update_value_by_overlap<side>(next_cell.get_cut().right_value,
                                        next_cell.template get_cut_right_polygon<coord_type>(),
                                        wave_polygon,
                                        *wave);
        }
        // = Update this cell ==============================

        // = Update tangentally affected cell ==============
        const size_t tangent_update_idx = [&curr_cell, &wave] {
          if constexpr (orientation == X) {
            if (wave->tangent_speed >= 0.0) {
              return curr_cell.top_idx;
            } else {
              return curr_cell.bottom_idx;
            }
          } else {
            if (wave->tangent_speed >= 0.0) {
              return curr_cell.right_idx;
            } else {
              return curr_cell.left_idx;
            }
          }
        }();

        if (curr_grid.is_cell(tangent_update_idx)) {
          auto& tangent_cell = next_grid[tangent_update_idx];
          if (tangent_cell.is_cartesian()) {
            const ActiveFloat overlap_area =
                0.5 * std::abs(wave->tangent_speed) * dt * std::abs(wave->normal_speed) * dt;
            const ActiveFloat cell_area =
                curr_cell.template dx<SIM_C>() * curr_cell.template dy<SIM_C>();

            tangent_cell.get_cartesian().value -=
                (overlap_area / cell_area) * wave->first_order_update;
          } else {
            const Geometry::Polygon<PointType> wave_polygon = calc_wave_polygon(
                *wave, curr_cell.template dx<SIM_C>(), curr_cell.template dy<SIM_C>(), dt);

            constexpr CoordType coord_type = PointType2CoordType<PointType>;
            // Left subcell
            update_value_by_overlap<side>(tangent_cell.get_cut().left_value,
                                          tangent_cell.template get_cut_left_polygon<coord_type>(),
                                          wave_polygon,
                                          *wave);
            // Right subcell
            update_value_by_overlap<side>(tangent_cell.get_cut().right_value,
                                          tangent_cell.template get_cut_right_polygon<coord_type>(),
                                          wave_polygon,
                                          *wave);
          }
        }
        // = Update tangentally affected cell ==============
      }
    } else if (idx == SAME_VALUE_INDEX || idx == ZERO_FLUX_INDEX) {
      // No-op
    } else {
      Igor::Debug("cell = {}", curr_cell);
      Igor::Panic("Invalid index on {} side.", side);
    }
    // - Handle side interface ---------------------------------------------------------------
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
  template <typename ActiveFloat, Point2D_c PointType>
  constexpr void
  move_wave_front(std::vector<PointType>& avg_new_shock_points,
                  const std::vector<FreeWave<ActiveFloat, PointType>>& internal_waves,
                  ActiveFloat dt) const noexcept {
    avg_new_shock_points.clear();
    if (internal_waves.empty()) { return; }

    avg_new_shock_points.push_back(
        internal_waves.front().begin +
        dt * internal_waves.front().normal_speed * internal_waves.front().normal +
        dt * internal_waves.front().tangent_speed * internal_waves.front().tangent);

    for (size_t i = 0; i < internal_waves.size() - 1; ++i) {
      const auto p1 = internal_waves[i].end +
                      dt * internal_waves[i].normal_speed * internal_waves[i].normal +
                      dt * internal_waves[i].tangent_speed * internal_waves[i].tangent;
      const auto p2 = internal_waves[i + 1].begin +
                      dt * internal_waves[i + 1].normal_speed * internal_waves[i + 1].normal +
                      dt * internal_waves[i + 1].tangent_speed * internal_waves[i + 1].tangent;
      avg_new_shock_points.push_back((p1 + p2) / 2);
    }
    avg_new_shock_points.push_back(
        internal_waves.back().end +
        dt * internal_waves.back().normal_speed * internal_waves.back().normal +
        dt * internal_waves.back().tangent_speed * internal_waves.back().tangent);

    // TODO: Is this necessary?
    // Remove first and last point if they are to close to the next point
    if ((avg_new_shock_points[0] - avg_new_shock_points[1]).norm() < 1e-2) {
      avg_new_shock_points.erase(std::next(std::begin(avg_new_shock_points)));
    }
    if (avg_new_shock_points.size() >= 2 && (avg_new_shock_points[avg_new_shock_points.size() - 1] -
                                             avg_new_shock_points[avg_new_shock_points.size() - 2])
                                                    .norm() < 1e-2) {
      avg_new_shock_points.erase(std::prev(std::end(avg_new_shock_points), 2));
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, typename PassiveFloat>
  void recalculate_cut_cell_values(
      UniformGrid<ActiveFloat, PassiveFloat>& next_grid,
      const UniformGrid<ActiveFloat, PassiveFloat>& curr_grid) const noexcept {

#pragma omp parallel for
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
            const auto left_intersect_area = (next_subcell_polygon & curr_cell_left_polygon).area();
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
            const auto left_intersect_area = (next_subcell_polygon & curr_cell_left_polygon).area();
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

    // Buffer for internal waves, defined here to avoid some allocations
    std::vector<FreeWave<ActiveFloat, PointType>> internal_waves{};
    internal_waves.reserve(grid.cut_cell_idxs().size());

    // Buffer for the points on the new buffer
    std::vector<PointType> avg_new_shock_points{};
    internal_waves.reserve(internal_waves.size() + 1);

    if (!(CFL_safety_factor > 0 && CFL_safety_factor <= 1)) {
      Igor::Warn("CFL_safety_factor must be in (0, 1], is {}", CFL_safety_factor);
      return std::nullopt;
    }
    auto& curr_grid = grid;
    auto next_grid  = grid;

#ifndef ZAP_NO_PARALLEL
    // ~ Chunks for parallelization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    struct Chunk {
      size_t xi_min, xi_max, yi_min, yi_max;
    };

    constexpr size_t NUM_CHUNKS   = 8;
    constexpr size_t NUM_CHUNKS_X = 4;
    constexpr size_t NUM_CHUNKS_Y = 2;
    static_assert(NUM_CHUNKS == NUM_CHUNKS_X * NUM_CHUNKS_Y);
    const auto chunk_width     = curr_grid.nx() / NUM_CHUNKS_X;
    const auto chunk_height    = curr_grid.ny() / NUM_CHUNKS_Y;
    const auto chunk_missing_x = curr_grid.nx() - chunk_width * NUM_CHUNKS_X;
    const auto chunk_missing_y = curr_grid.ny() - chunk_height * NUM_CHUNKS_Y;

    std::array<Chunk, NUM_CHUNKS> chunks;
    for (size_t chunk_yi = 0; chunk_yi < NUM_CHUNKS_Y; ++chunk_yi) {
      for (size_t chunk_xi = 0; chunk_xi < NUM_CHUNKS_X; ++chunk_xi) {
        const auto i = chunk_xi + chunk_yi * NUM_CHUNKS_X;
        chunks[i]    = Chunk{
               .xi_min = chunk_xi * chunk_width,
               .xi_max =
                (chunk_xi + 1) * chunk_width + chunk_missing_x * (chunk_xi + 1 == NUM_CHUNKS_X),
               .yi_min = chunk_yi * chunk_height,
               .yi_max =
                (chunk_yi + 1) * chunk_height + chunk_missing_y * (chunk_yi + 1 == NUM_CHUNKS_Y),
        };
      }
    }

#pragma omp parallel
    {
#pragma omp single
      if (omp_get_num_threads() != static_cast<int>(NUM_CHUNKS)) {
        Igor::Warn("The solver is optimized for {} threads, found {} threads.",
                   NUM_CHUNKS,
                   omp_get_num_threads());
      }
    }
    // ~ Chunks for parallelization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#endif  // ZAP_NO_PARALLEL

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

      // = Pre-calculate internal interface waves ==================================================
      internal_waves.resize(curr_grid.cut_cell_idxs().size());
#pragma omp parallel for
      for (size_t i = 0; i < curr_grid.cut_cell_idxs().size(); ++i) {
        const size_t cell_idx = curr_grid.cut_cell_idxs()[i];
        const auto& curr_cell = curr_grid[cell_idx];
        const auto internal_interface =
            get_internal_interface<ActiveFloat, PassiveFloat, PointType>(curr_cell);
        auto wave = normal_wave<FREE>(
            internal_interface, curr_cell.template dx<SIM_C>(), curr_cell.template dy<SIM_C>(), dt);
        IGOR_ASSERT(wave.has_value(),
                    "Internal wave must exist, i.e. internal interface must have length > EPS but "
                    "internal interface has begin={} and end={}",
                    internal_interface.begin,
                    internal_interface.end);
        internal_waves[i] = std::move(*wave);
      }
      // = Pre-calculate internal interface waves ==================================================

#ifndef ZAP_STATIC_CUT
      // = Move the wave front =====================================================================
      if (!curr_grid.cut_cell_idxs().empty()) {
        next_grid.merge_cut_cells();

        move_wave_front(avg_new_shock_points, internal_waves, dt);
        if (!next_grid.template cut_piecewise_linear<extend_type>(avg_new_shock_points)) {
          Igor::Warn("Could not cut on new shock curve.");
          return std::nullopt;
        }
      }

      recalculate_cut_cell_values(next_grid, curr_grid);
      // = Move the wave front =====================================================================

      // TODO: Consider purging cuts that are no longer at the shock front using the
      //       RK-Jump-Condition
#endif  // ZAP_STATIC_CUT
      assert(curr_grid.cut_cell_idxs().size() == internal_waves.size());

      // TODO: What happens when a subcell has area 0?
      //   -> cannot "uncut" the cell because that would loose shock position information

      // TODO: Do parallelization that handles tangential correction

#ifdef ZAP_NO_PARALLEL
      for (size_t cell_idx = 0; cell_idx < curr_grid.size(); ++cell_idx) {
        const auto& curr_cell = curr_grid[cell_idx];
        auto& next_cell       = next_grid[cell_idx];

        handle_side_iterface<PointType, LEFT>(next_cell, curr_cell, next_grid, curr_grid, dt);
        handle_side_iterface<PointType, RIGHT>(next_cell, curr_cell, next_grid, curr_grid, dt);
        handle_side_iterface<PointType, BOTTOM>(next_cell, curr_cell, next_grid, curr_grid, dt);
        handle_side_iterface<PointType, TOP>(next_cell, curr_cell, next_grid, curr_grid, dt);
      }
      // = Handle outer interfaces: LEFT, RIGHT, BOTTOM, and TOP ===================================

      // = Handle internal interface ===============================================================
      assert(curr_grid.cut_cell_idxs().size() == internal_waves.size());
      for (size_t i = 0; i < internal_waves.size(); ++i) {
        const auto cell_idx   = curr_grid.cut_cell_idxs()[i];
        const auto& curr_cell = curr_grid[cell_idx];
        assert(curr_cell.is_cut());
        auto& next_cell = next_grid[cell_idx];
        assert(next_cell.is_cartesian() || next_cell.is_cut());

        const auto& internal_wave  = internal_waves[i];
        const auto cell_neighbours = curr_grid.get_existing_neighbours(cell_idx);
        // Neighbouring cells
        for (size_t neighbour_idx : cell_neighbours) {
          update_cell_by_free_wave(next_grid[neighbour_idx], internal_wave, dt);
        }
        // This cell
        update_cell_by_free_wave(next_cell, internal_wave, dt);
      }
      // = Handle internal interface ===============================================================
#else
      // ~ Parallelization of grid update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      auto do_update = [&](size_t xi, size_t yi) {
        const size_t cell_idx = curr_grid.to_vec_idx(xi, yi);
        const auto& curr_cell = curr_grid[cell_idx];
        auto& next_cell       = next_grid[cell_idx];

        // = Handle outer interfaces: LEFT, RIGHT, BOTTOM, and TOP =================================
        handle_side_iterface<PointType, LEFT>(next_cell, curr_cell, next_grid, curr_grid, dt);
        handle_side_iterface<PointType, RIGHT>(next_cell, curr_cell, next_grid, curr_grid, dt);
        handle_side_iterface<PointType, BOTTOM>(next_cell, curr_cell, next_grid, curr_grid, dt);
        handle_side_iterface<PointType, TOP>(next_cell, curr_cell, next_grid, curr_grid, dt);
        // = Handle outer interfaces: LEFT, RIGHT, BOTTOM, and TOP =================================

        // = Handle internal interface =============================================================
        // TODO: Why does the result change if I merge the loops?
        if (curr_cell.is_cut()) {
          // TODO: Can we do better than linear search?
          const auto it = std::find(
              curr_grid.cut_cell_idxs().cbegin(), curr_grid.cut_cell_idxs().cend(), cell_idx);
          IGOR_ASSERT(it != curr_grid.cut_cell_idxs().cend(),
                      "Did not find cell_idx {} in array of cut cell indices {}",
                      cell_idx,
                      curr_grid.cut_cell_idxs());
          const auto iw_idx =
              static_cast<size_t>(std::distance(curr_grid.cut_cell_idxs().cbegin(), it));

          const auto& internal_wave  = internal_waves[iw_idx];
          const auto cell_neighbours = curr_grid.get_existing_neighbours(cell_idx);
          // Neighbouring cells
          for (size_t neighbour_idx : cell_neighbours) {
            update_cell_by_free_wave(next_grid[neighbour_idx], internal_wave, dt);
          }
          // This cell
          update_cell_by_free_wave(next_cell, internal_wave, dt);
        }
        // = Handle internal interface =============================================================
      };

      auto update_chunk = [&](const Chunk& chunk) {
        for (size_t yi = chunk.yi_min + 1; yi < chunk.yi_max - 1; ++yi) {
          for (size_t xi = chunk.xi_min + 1; xi < chunk.xi_max - 1; ++xi) {
            do_update(xi, yi);
          }
        }
      };

      auto stich_chunk = [&](const Chunk& chunk) {
        for (size_t xi = chunk.xi_min; xi < chunk.xi_max; ++xi) {
          do_update(xi, chunk.yi_min);
          do_update(xi, chunk.yi_max - 1);
        }
        for (size_t yi = chunk.yi_min; yi < chunk.yi_max; ++yi) {
          do_update(chunk.xi_min, yi);
          do_update(chunk.xi_max - 1, yi);
        }
      };

#pragma omp parallel
#pragma omp single
      for (const auto& chunk : chunks) {
#pragma omp task
        update_chunk(chunk);
      }

      for (const auto& chunk : chunks) {
        stich_chunk(chunk);
      }
      // ~ Parallelization of grid update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#endif  // ZAP_NO_PARALLEL

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
