#ifndef ZAP_CELL_BASED_SOLVER_HPP_
#define ZAP_CELL_BASED_SOLVER_HPP_

#include <numbers>
#include <numeric>

#include "CellBased/Geometry.hpp"
#include "CellBased/Grid.hpp"
#include "CellBased/Interface.hpp"

namespace Zap::CellBased {

// TODO: Use grid coordinates for all geometry operations

template <typename Float, int DIM>
requires(DIM > 0)
constexpr void get_eigen_decomp(const Eigen::Matrix<Float, DIM, DIM>& mat,
                                Eigen::Matrix<Float, DIM, DIM>& eig_vals,
                                Eigen::Matrix<Float, DIM, DIM>& eig_vecs) noexcept {
  if (DIM == 1) {
    eig_vals = mat;
    eig_vecs = Eigen::Matrix<Float, DIM, DIM>::Identity();
  } else if (DIM == 2) {
    assert(std::pow((mat(0, 0) + mat(1, 1)) / 2, 2) >= mat.determinant() &&
           "mat has complex eigenvalues");
    Igor::Todo("DIM = 2 is not implemented yet.");
    std::unreachable();
  } else {
    Igor::Todo("DIM = {} is not implemented yet.", DIM);
    std::unreachable();
  }
}

template <ExtendType extend_type, typename A, typename B>
class Solver {
  A m_A;
  B m_B;

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, typename PassiveFloat, size_t DIM>
  [[nodiscard]] constexpr auto
  cfl_factor(const UniformGrid<ActiveFloat, PassiveFloat, DIM>& grid) const noexcept
      -> ActiveFloat {
    return std::transform_reduce(
        std::cbegin(grid.cells()),
        std::cend(grid.cells()),
        ActiveFloat{0},
        [](const auto& a, const auto& b) { return std::max(a, b); },
        [this](const Cell<ActiveFloat, PassiveFloat, DIM>& cell) {
          assert(cell.is_cartesian() || cell.is_cut());
          if (cell.is_cartesian()) {
            return std::max(m_A.max_abs_eig_val(cell.get_cartesian().value),
                            m_B.max_abs_eig_val(cell.get_cartesian().value));
          } else {
            return std::max({m_A.max_abs_eig_val(cell.get_cut().left_value),
                             m_A.max_abs_eig_val(cell.get_cut().right_value),
                             m_B.max_abs_eig_val(cell.get_cut().left_value),
                             m_B.max_abs_eig_val(cell.get_cut().right_value)});
          }
        });
  }

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, size_t DIM, Point2D_c PointType>
  requires(DIM > 0)
  struct WaveProperties {
    Eigen::Vector<ActiveFloat, DIM> value;
    Geometry::Polygon<PointType> polygon;
  };

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, typename PassiveFloat, size_t DIM, Point2D_c PointType>
  requires(DIM > 0)
  [[nodiscard]] constexpr auto
  calculate_interface(const FullInterface<ActiveFloat, DIM, PointType>& interface,
                      ActiveFloat dt,
                      PassiveFloat scale_x,
                      PassiveFloat scale_y) const noexcept
      -> SmallVector<WaveProperties<ActiveFloat, DIM, PointType>> {
    if ((interface.end - interface.begin).norm() < EPS<PassiveFloat>) { return {}; }

    // Vector tangential to interface
    const PointType tangent_vector = (interface.end - interface.begin).normalized();

    // Angle between cut_vector and y-axis (0, 1)
    //     cut_angle = arccos(cut_vector^T * (0, 1) / ||cut_vector|| * ||(0, 1)||)
    // <=> cut_angle = arccos(cut_vector_1 / ||cut_vector||)
    ActiveFloat interface_angle = std::acos(tangent_vector.y);
    if (tangent_vector.x > 0) {
      interface_angle = 2 * std::numbers::pi_v<PassiveFloat> - interface_angle;
    }

    // Vector normal to cut
    PointType normal_vector{std::cos(interface_angle), std::sin(interface_angle)};
    IGOR_ASSERT(std::abs(tangent_vector.dot(normal_vector)) <= EPS<PassiveFloat>,
                "tangent_vector {} is not orthogonal to normal_vector {}, tangent_vector^T * "
                "normal_vector = {}",
                tangent_vector,
                normal_vector,
                tangent_vector.dot(normal_vector));
    // Scale normal vector when operating in grid coordinates
    if constexpr (is_GridCoord_v<PointType>) {
      normal_vector.x *= scale_x;
      normal_vector.y *= scale_y;
    }

    const Eigen::Vector<ActiveFloat, DIM> u_mid =
        (interface.left_value + interface.right_value) / 2;

    // Matrix for rotated PDE in normal direction to cut
    // TODO: Why cos(interface_angle) * A + sin(interface_angle) * B and not -?
    const Eigen::Matrix<ActiveFloat, DIM, DIM> eta_mat =
        std::cos(interface_angle) * m_A.mat(u_mid) + std::sin(interface_angle) * m_B.mat(u_mid);

    Eigen::Matrix<ActiveFloat, DIM, DIM> eig_vals;
    Eigen::Matrix<ActiveFloat, DIM, DIM> eig_vecs;
    get_eigen_decomp(eta_mat, eig_vals, eig_vecs);
    const Eigen::Matrix<ActiveFloat, DIM, DIM> wave_lengths = eig_vals * dt;

    // Eigen expansion of jump; wave strength
    const Eigen::Vector<ActiveFloat, DIM> alpha =
        eig_vecs.inverse() * (interface.right_value - interface.left_value);

    SmallVector<WaveProperties<ActiveFloat, DIM, PointType>> waves(DIM);

    for (Eigen::Index p = 0; p < static_cast<Eigen::Index>(DIM); ++p) {
      auto& wave = waves[static_cast<size_t>(p)];

      wave.value = alpha(p, p) * eig_vecs.col(p);

      wave.polygon = Geometry::Polygon<PointType>{{
          interface.begin,
          interface.begin + normal_vector * wave_lengths(p, p),
          interface.end,
          interface.end + normal_vector * wave_lengths(p, p),
      }};
    }

    return waves;
  }

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, typename PassiveFloat, size_t DIM, Point2D_c PointType>
  constexpr void update_cell(Cell<ActiveFloat, PassiveFloat, DIM>& cell,
                             const WaveProperties<ActiveFloat, DIM, PointType>& wave) {
    // TODO: Do something about periodic boundary conditions
    assert(cell.is_cartesian() || cell.is_cut());
    auto update_value = [&](Eigen::Vector<ActiveFloat, DIM>& value,
                            const Geometry::Polygon<PointType>& cell_polygon) {
      const auto cell_area = cell_polygon.area();
      assert(cell_area > 0 || std::abs(cell_area) <= EPS<PassiveFloat>);
      if (std::abs(cell_area) <= EPS<PassiveFloat>) { return; }

      const auto intersect      = Geometry::intersection(cell_polygon, wave.polygon);
      const auto intersect_area = intersect.area();
      assert(intersect_area >= 0 || std::abs(intersect_area) < EPS<PassiveFloat>);
      IGOR_ASSERT(intersect_area - cell_area <= 50 * EPS<PassiveFloat>,
                  "Expected area intersection of intersection to be smaller or equal to the area "
                  "of the cell, but intersection area is {} and cell area is {}, intersection area "
                  "is {} larger than cell area",
                  intersect_area,
                  cell_area,
                  static_cast<ActiveFloat>(intersect_area - cell_area));

      value -= (intersect_area / cell_area) * wave.value;
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
  template <typename ActiveFloat, typename PassiveFloat, size_t DIM>
  constexpr void apply_side_interfaces(const UniformGrid<ActiveFloat, PassiveFloat, DIM>& curr_grid,
                                       Cell<ActiveFloat, PassiveFloat, DIM>& next_cell,
                                       const Cell<ActiveFloat, PassiveFloat, DIM>& curr_cell,
                                       size_t idx,
                                       Side side,
                                       ActiveFloat dt) {
    if (curr_grid.is_cell(idx)) {
      const auto& other_cell = curr_grid[idx];
      const auto interfaces =
          get_shared_interfaces<ActiveFloat, PassiveFloat, DIM, GridCoord<ActiveFloat>>(
              curr_cell, other_cell, side);
      for (const auto& interface : interfaces) {
        const auto waves =
            calculate_interface(interface, dt, curr_grid.scale_x(), curr_grid.scale_y());
        for (const auto& wave : waves) {
          update_cell(next_cell, wave);
        }
      }
    } else if (idx == SAME_VALUE_INDEX) {
      // TODO: Noop?
    } else if (idx == ZERO_FLUX_INDEX) {
      // TODO: Noop?
    } else if (idx == NULL_INDEX) {
      std::stringstream s;
      s << curr_cell;
      Igor::Debug("cell = {}", s.str());
      Igor::Panic("cell has NULL index.");
    } else {
      std::stringstream s;
      s << curr_cell;
      Igor::Debug("cell = {}", s.str());
      Igor::Panic("cell has unknown index with value {}.", idx);
    }
  };

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, typename PassiveFloat, size_t DIM, Point2D_c PointType>
  [[nodiscard]] constexpr auto
  get_internal_interface(const Cell<ActiveFloat, PassiveFloat, DIM>& cell) const noexcept
      -> FullInterface<ActiveFloat, DIM, PointType> {
    assert(cell.is_cut());

    constexpr CoordType coord_type = PointType2CoordType<PointType>;

    return FullInterface<ActiveFloat, DIM, PointType>{
        .left_value  = cell.get_cut().left_value,
        .right_value = cell.get_cut().right_value,
        .begin       = cell.template cut1<coord_type>(),
        .end         = cell.template cut2<coord_type>(),
    };
  }

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat, typename PassiveFloat, size_t DIM>
  [[nodiscard]] constexpr auto
  move_wave_front(const UniformGrid<ActiveFloat, PassiveFloat, DIM>& curr_grid,
                  [[maybe_unused]] const UniformGrid<ActiveFloat, PassiveFloat, DIM>& next_grid,
                  ActiveFloat dt) const noexcept
      -> std::optional<std::vector<GridCoord<ActiveFloat>>> {
    using PointType = GridCoord<ActiveFloat>;

    // Move old cuts according to strongest wave
    std::vector<std::pair<PointType, PointType>> new_shock_points;
    for (size_t cell_idx : curr_grid.cut_cell_idxs()) {
      const auto& curr_cell = curr_grid[cell_idx];
      assert(curr_cell.is_cut());

      const auto interface =
          get_internal_interface<ActiveFloat, PassiveFloat, DIM, PointType>(curr_cell);

      const auto tangent_vector = (interface.end - interface.begin).normalized();
      assert(std::abs(tangent_vector.norm() - 1) <= 1e-8);

      ActiveFloat interface_angle = std::acos(tangent_vector.y);
      if (tangent_vector.x > 0) {
        interface_angle = 2 * std::numbers::pi_v<PassiveFloat> - interface_angle;
      }
      PointType normal_vector{std::cos(interface_angle), std::sin(interface_angle)};
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

      const Eigen::Vector<ActiveFloat, DIM> u_mid =
          (interface.left_value + interface.right_value) / 2;

      const Eigen::Matrix<ActiveFloat, DIM, DIM> eta_mat =
          std::cos(interface_angle) * m_A.mat(u_mid) + std::sin(interface_angle) * m_B.mat(u_mid);

      Eigen::Matrix<ActiveFloat, DIM, DIM> eig_vals;
      Eigen::Matrix<ActiveFloat, DIM, DIM> eig_vecs;
      get_eigen_decomp(eta_mat, eig_vals, eig_vecs);
      // Wave strength
      const Eigen::Vector<ActiveFloat, DIM> alpha =
          eig_vecs.inverse() * (interface.right_value - interface.left_value);

      Eigen::Index strong_wave_idx = 0;
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        // TODO: Abs or not?
        if (std::abs(alpha(i)) > std::abs(alpha(strong_wave_idx))) { strong_wave_idx = i; }
      }

      const ActiveFloat wave_length = eig_vals(strong_wave_idx) * dt;

      new_shock_points.emplace_back(interface.begin + normal_vector * wave_length,
                                    interface.end + normal_vector * wave_length);
    }

    std::vector<PointType> avg_new_shock_points(new_shock_points.size() + 1);
    assert(avg_new_shock_points.size() > 2);
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
  constexpr Solver(A a, B b)
      : m_A(std::move(a)),
        m_B(std::move(b)) {}

  // -----------------------------------------------------------------------------------------------
  template <typename ActiveFloat,
            typename PassiveFloat,
            size_t DIM,
            typename GridWriter,
            typename TimeWriter>
  [[nodiscard]] auto solve(UniformGrid<ActiveFloat, PassiveFloat, DIM> grid,
                           PassiveFloat tend,
                           GridWriter& grid_writer,
                           TimeWriter& time_writer,
                           PassiveFloat CFL_safety_factor = 0.5) noexcept
      -> std::optional<UniformGrid<ActiveFloat, PassiveFloat, DIM>> {
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

#ifdef ZAP_ASSERT_CONSERVATIVE
      const auto mass_before_moving = next_grid.mass();
#endif  // ZAP_ASSERT_CONSERVATIVE

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
          // TODO: Change the type for one to SIM to see compilation error
          const auto curr_cell_left_polygon  = curr_cell.template get_cut_left_polygon<GRID_C>();
          const auto curr_cell_right_polygon = curr_cell.template get_cut_right_polygon<GRID_C>();

          // Left subcell
          {
            const auto next_subcell_polygon = next_cell.template get_cut_left_polygon<GRID_C>();
            assert(next_subcell_polygon.area() > 0 ||
                   std::abs(next_subcell_polygon.area()) <= EPS<PassiveFloat>);
            if (std::abs(next_subcell_polygon.area()) > EPS<PassiveFloat>) {
              const auto left_intersect_area =
                  Geometry::intersection(next_subcell_polygon, curr_cell_left_polygon).area();
              const auto right_intersect_area =
                  Geometry::intersection(next_subcell_polygon, curr_cell_right_polygon).area();
              assert(std::abs(left_intersect_area + right_intersect_area -
                              next_subcell_polygon.area()) < 1e-6);

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
                  Geometry::intersection(next_subcell_polygon, curr_cell_left_polygon).area();
              const auto right_intersect_area =
                  Geometry::intersection(next_subcell_polygon, curr_cell_right_polygon).area();
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

#ifdef ZAP_ASSERT_CONSERVATIVE
      const auto mass_after_moving = next_grid.mass();
      IGOR_ASSERT(
          (mass_before_moving - mass_after_moving).norm() <= EPS<PassiveFloat>,
          "Mass has to be conserved while moving the cuts, but mass before is {} and after {}.",
          mass_before_moving,
          mass_after_moving);
#endif  // ZAP_ASSERT_CONSERVATIVE
#endif  // ZAP_STATIC_CUT

      // TODO: What happens when a subcell has area 0?
      //   -> cannot "uncut" the cell because that would loose shock position information

#pragma omp parallel for
      for (size_t cell_idx = 0; cell_idx < curr_grid.size(); ++cell_idx) {
        const auto& curr_cell = curr_grid[cell_idx];
        auto& next_cell       = next_grid[cell_idx];

        // = Cut cell ==============================================================================
        if (curr_cell.is_cut() || next_cell.is_cut()) {
          // Left interface
          apply_side_interfaces(curr_grid, next_cell, curr_cell, curr_cell.left_idx, LEFT, dt);
          // Right interface
          apply_side_interfaces(curr_grid, next_cell, curr_cell, curr_cell.right_idx, RIGHT, dt);
          // Bottom interface
          apply_side_interfaces(curr_grid, next_cell, curr_cell, curr_cell.bottom_idx, BOTTOM, dt);
          // Top interface
          apply_side_interfaces(curr_grid, next_cell, curr_cell, curr_cell.top_idx, TOP, dt);
        }
        // = Cut cell ==============================================================================
        // = Cartesian cell ========================================================================
        else {
          // - Handle left interface ---------------------------------------------------------------
          if (curr_grid.is_cell(curr_cell.left_idx)) {
            const auto& other_cell = curr_grid[curr_cell.left_idx];
            const auto interfaces =
                get_shared_interfaces<ActiveFloat, PassiveFloat, DIM, SimCoord<ActiveFloat>>(
                    curr_cell, other_cell, LEFT);

            for (const auto& interface : interfaces) {
              const Eigen::Vector<ActiveFloat, DIM> u_mid =
                  (interface.left_value + interface.right_value) / 2;
              Eigen::Matrix<ActiveFloat, DIM, DIM> eig_vals = m_A.eig_vals(u_mid);
              Eigen::Matrix<ActiveFloat, DIM, DIM> eig_vecs = m_A.eig_vecs(u_mid);

              bool all_zero = true;
              for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
                if (std::abs(eig_vals(i, i)) >= 1e-6) { all_zero = false; }
              }
              if (all_zero) { continue; }

              assert(std::abs(eig_vecs.determinant()) >= 1e-8);
              for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
                eig_vals(i, i) = eig_vals(i, i) > 0 ? eig_vals(i, i) : 0;
              }
              const Eigen::Matrix<ActiveFloat, DIM, DIM> A_plus =
                  eig_vecs * eig_vals * eig_vecs.inverse();

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
                get_shared_interfaces<ActiveFloat, PassiveFloat, DIM, SimCoord<ActiveFloat>>(
                    curr_cell, other_cell, RIGHT);

            for (const auto& interface : interfaces) {
              const Eigen::Vector<ActiveFloat, DIM> u_mid =
                  (interface.left_value + interface.right_value) / 2;
              Eigen::Matrix<ActiveFloat, DIM, DIM> eig_vals = m_A.eig_vals(u_mid);
              Eigen::Matrix<ActiveFloat, DIM, DIM> eig_vecs = m_A.eig_vecs(u_mid);

              bool all_zero = true;
              for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
                if (std::abs(eig_vals(i, i)) >= 1e-6) { all_zero = false; }
              }
              if (all_zero) { continue; }

              assert(std::abs(eig_vecs.determinant()) >= 1e-8);
              for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
                eig_vals(i, i) = eig_vals(i, i) < 0 ? eig_vals(i, i) : 0;
              }
              const Eigen::Matrix<ActiveFloat, DIM, DIM> A_minus =
                  eig_vecs * eig_vals * eig_vecs.inverse();
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
                get_shared_interfaces<ActiveFloat, PassiveFloat, DIM, SimCoord<ActiveFloat>>(
                    curr_cell, other_cell, BOTTOM);

            for (const auto& interface : interfaces) {
              const Eigen::Vector<ActiveFloat, DIM> u_mid =
                  (interface.left_value + interface.right_value) / 2;
              Eigen::Matrix<ActiveFloat, DIM, DIM> eig_vals = m_B.eig_vals(u_mid);
              Eigen::Matrix<ActiveFloat, DIM, DIM> eig_vecs = m_B.eig_vecs(u_mid);

              bool all_zero = true;
              for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
                if (std::abs(eig_vals(i, i)) >= 1e-6) { all_zero = false; }
              }
              if (all_zero) { continue; }

              assert(std::abs(eig_vecs.determinant()) >= 1e-8);
              for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
                eig_vals(i, i) = eig_vals(i, i) > 0 ? eig_vals(i, i) : 0;
              }
              const Eigen::Matrix<ActiveFloat, DIM, DIM> A_plus =
                  eig_vecs * eig_vals * eig_vecs.inverse();

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
                get_shared_interfaces<ActiveFloat, PassiveFloat, DIM, SimCoord<ActiveFloat>>(
                    curr_cell, other_cell, TOP);

            for (const auto& interface : interfaces) {
              const Eigen::Vector<ActiveFloat, DIM> u_mid =
                  (interface.left_value + interface.right_value) / 2;
              Eigen::Matrix<ActiveFloat, DIM, DIM> eig_vals = m_B.eig_vals(u_mid);
              Eigen::Matrix<ActiveFloat, DIM, DIM> eig_vecs = m_B.eig_vecs(u_mid);

              bool all_zero = true;
              for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
                if (std::abs(eig_vals(i, i)) >= 1e-6) { all_zero = false; }
              }
              if (all_zero) { continue; }

              assert(std::abs(eig_vecs.determinant()) >= 1e-8);
              for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
                eig_vals(i, i) = eig_vals(i, i) < 0 ? eig_vals(i, i) : 0;
              }
              const Eigen::Matrix<ActiveFloat, DIM, DIM> A_minus =
                  eig_vecs * eig_vals * eig_vecs.inverse();
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

#ifdef ZAP_ASSERT_CONSERVATIVE
      const auto mass_after_inter_cell_update = next_grid.mass();
#endif  // ZAP_ASSERT_CONSERVATIVE

      for (size_t cell_idx : curr_grid.cut_cell_idxs()) {
        const auto& curr_cell = curr_grid[cell_idx];
        assert(curr_cell.is_cut());
        auto& next_cell = next_grid[cell_idx];
        assert(next_cell.is_cartesian() || next_cell.is_cut());

        // - Handle internal interface -----------------------------------------------------------
        const auto internal_interface =
            get_internal_interface<ActiveFloat, PassiveFloat, DIM, GridCoord<ActiveFloat>>(
                curr_cell);
        const auto internal_waves =
            calculate_interface<ActiveFloat, PassiveFloat, DIM, GridCoord<ActiveFloat>>(
                internal_interface, dt, curr_grid.scale_x(), curr_grid.scale_y());

        const auto cell_neighbours = curr_grid.get_existing_neighbours(cell_idx);
        for (const auto& wave : internal_waves) {
          // Neighbouring cells
          for (size_t neighbour_idx : cell_neighbours) {
            update_cell(next_grid[neighbour_idx], wave);
          }
          // This cell
          update_cell(next_cell, wave);
        }
        // - Handle internal interface -----------------------------------------------------------
      }

#ifdef ZAP_ASSERT_CONSERVATIVE
      const auto mass_after_inner_cell_update = next_grid.mass();
      IGOR_ASSERT((mass_after_inner_cell_update - mass_after_moving).norm() <= EPS<PassiveFloat>,
                  "Mass has to be conserved while doing the update, but mass before is {} and "
                  "after {}. Difference is {}. Mass after just doing inter-cell updates {}. Number "
                  "of small subcells is {}.",
                  mass_after_moving,
                  mass_after_inner_cell_update,
                  (mass_after_inner_cell_update - mass_after_moving).eval(),
                  mass_after_inter_cell_update,
                  std::count_if(std::cbegin(next_grid.cells()),
                                std::cend(next_grid.cells()),
                                [](const Cell<ActiveFloat, PassiveFloat, DIM>& cell) {
                                  return cell.is_cut() &&
                                         (cell.template get_cut_left_polygon<SIM_C>().area() <
                                              EPS<PassiveFloat> ||
                                          cell.template get_cut_right_polygon<SIM_C>().area() <
                                              EPS<PassiveFloat>);
                                }));
#endif  // ZAP_ASSERT_CONSERVATIVE

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

// -------------------------------------------------------------------------------------------------
template <ExtendType extend_type, typename A, typename B>
[[nodiscard]] constexpr auto make_solver(A&& a, B&& b) noexcept -> Solver<extend_type, A, B> {
  return Solver<extend_type, A, B>(std::forward<A>(a), std::forward<B>(b));
}

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_SOLVER_HPP_
