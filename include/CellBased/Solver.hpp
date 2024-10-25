#ifndef ZAP_CELL_BASED_SOLVER_HPP_
#define ZAP_CELL_BASED_SOLVER_HPP_

#include <numbers>
#include <numeric>
#include <unordered_set>

#include "CellBased/Geometry.hpp"
#include "CellBased/Grid.hpp"
#include "CellBased/Interface.hpp"

namespace Zap::CellBased {

// TODO: Use grid coordinates for all geometry operations

template <typename Float, int DIM>
requires(DIM > 0)
struct EigenDecomp {
  Eigen::Matrix<Float, DIM, DIM> vals;
  Eigen::Matrix<Float, DIM, DIM> vecs;
};

template <typename Float, int DIM>
requires(DIM > 0)
[[nodiscard]] constexpr auto get_eigen_decomp(const Eigen::Matrix<Float, DIM, DIM>& mat) noexcept
    -> EigenDecomp<Float, DIM> {
  if (DIM == 1) {
    return {.vals = mat, .vecs = Eigen::Matrix<Float, DIM, DIM>::Identity()};
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
    Eigen::Vector<ActiveFloat, DIM> mass;
    Geometry::Polygon<PointType> polygon;
    ActiveFloat sign;
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
    // TODO: Do tangential splitting correction.

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

    const Eigen::Vector<ActiveFloat, DIM> u_mid =
        (interface.left_value + interface.right_value) / 2;

    // Matrix for rotated PDE in normal direction to cut
    // TODO: Why cos(interface_angle) * A + sin(interface_angle) * B and not -?
    // Equivalent to `cos(interface_angle) * A + sin(interface_angle) * B`, but do not actually
    // calculate the interface angle
    const Eigen::Matrix<ActiveFloat, DIM, DIM> normal_mat =
        tangent_vector.y * m_A.mat(u_mid) +  //
        -sign(tangent_vector.x) * std::sqrt(1 - tangent_vector.y * tangent_vector.y) *
            m_B.mat(u_mid);

    const auto normal_eig                                          = get_eigen_decomp(normal_mat);
    const Eigen::Matrix<ActiveFloat, DIM, DIM> normal_wave_lengths = normal_eig.vals * dt;

    // Eigen expansion of jump; wave strength
    const Eigen::Vector<ActiveFloat, DIM> normal_alpha =
        normal_eig.vecs.inverse() * (interface.right_value - interface.left_value);

#ifndef ZAP_TANGENTIAL_CORRECTION
    SmallVector<WaveProperties<ActiveFloat, DIM, PointType>> waves(DIM);
    for (Eigen::Index p = 0; p < static_cast<Eigen::Index>(DIM); ++p) {
      auto& wave = waves[static_cast<size_t>(p)];

      wave.mass = normal_alpha(p) * normal_eig.vecs.col(p);

      wave.polygon = Geometry::Polygon<PointType>{{
          interface.begin,
          interface.begin + scaled_normal_vector * normal_wave_lengths(p, p),
          interface.end,
          interface.end + scaled_normal_vector * normal_wave_lengths(p, p),
      }};

      wave.sign = sign(normal_eig.vals(p, p));
    }
#else
    // Matrix for rotated PDE in tangential direction to cut
    // Equivalent to `-std::sin(interface_angle) * A + std::cos(interface_angle) * B`
    const Eigen::Matrix<ActiveFloat, DIM, DIM> tangent_mat =
        sign(tangent_vector.x) * std::sqrt(1 - tangent_vector.y * tangent_vector.y) *
            m_A.mat(u_mid) +
        tangent_vector.y * m_B.mat(u_mid);

    const auto tangent_eig                                          = get_eigen_decomp(tangent_mat);
    const Eigen::Matrix<ActiveFloat, DIM, DIM> tangent_wave_lengths = tangent_eig.vals * dt;

    SmallVector<WaveProperties<ActiveFloat, DIM, PointType>> waves(DIM * DIM);
    size_t idx = 0;
    for (Eigen::Index p = 0; p < static_cast<Eigen::Index>(DIM); ++p) {
      for (Eigen::Index q = 0; q < static_cast<Eigen::Index>(DIM); ++q) {
        auto& wave = waves[idx];

        // Eigen expansion of jump; wave strength
        const Eigen::Vector<ActiveFloat, DIM> tangent_beta =
            tangent_eig.vecs.inverse() * (normal_alpha(p) * normal_eig.vecs.col(p));

        wave.mass = tangent_beta(q) * tangent_eig.vecs.col(q);

        wave.polygon = Geometry::Polygon<PointType>{{
            interface.begin,
            interface.begin + scaled_normal_vector * normal_wave_lengths(p, p) +
                scaled_tangent_vector * tangent_wave_lengths(q, q),
            interface.end,
            interface.end + scaled_normal_vector * normal_wave_lengths(p, p) +
                scaled_tangent_vector * tangent_wave_lengths(q, q),
        }};

        wave.sign = static_cast<PassiveFloat>(normal_eig.vals(p, p) > EPS<PassiveFloat>) -
                    static_cast<PassiveFloat>(normal_eig.vals(p, p) < -EPS<PassiveFloat>);

        idx += 1;
      }
    }
#endif  // ZAP_NOT_TANGENTIAL_CORRECTION

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
  template <typename ActiveFloat, typename PassiveFloat, size_t DIM>
  constexpr void
  apply_side_interfaces_uncorrected(const UniformGrid<ActiveFloat, PassiveFloat, DIM>& curr_grid,
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

      const Eigen::Vector<ActiveFloat, DIM> u_mid =
          (interface.left_value + interface.right_value) / 2;

      // Equivalent to `cos(interface_angle) * A + sin(interface_angle) * B`, but do not actually
      // calculate the interface angle
      const Eigen::Matrix<ActiveFloat, DIM, DIM> eta_mat =
          tangent_vector.y * m_A.mat(u_mid) +
          -sign(tangent_vector.x) * std::sqrt(1 - tangent_vector.y * tangent_vector.y) *
              m_B.mat(u_mid);

      const auto eig = get_eigen_decomp(eta_mat);
      // Wave strength
      const Eigen::Vector<ActiveFloat, DIM> alpha =
          eig.vecs.inverse() * (interface.right_value - interface.left_value);

      Eigen::Index strong_wave_idx = 0;
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        // TODO: Abs or not?
        if (std::abs(alpha(i)) > std::abs(alpha(strong_wave_idx))) { strong_wave_idx = i; }
      }

      const ActiveFloat wave_length = eig.vals(strong_wave_idx) * dt;

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
        // TODO: Calculate outer interfaces

        // - Left outer interface ------------------------------------------------------------------
        if (curr_grid.is_cell(curr_cell.left_idx) &&
            !calculated_interfaces.contains({cell_idx, curr_cell.left_idx})) {
          const auto& other_cell = curr_grid[curr_cell.left_idx];
          const auto interfaces =
              get_shared_interfaces<ActiveFloat, PassiveFloat, DIM, GridCoord<ActiveFloat>>(
                  curr_cell, other_cell, LEFT);

          for (const auto& interface : interfaces) {
            const auto waves =
                calculate_interface(interface, dt, curr_grid.scale_x(), curr_grid.scale_y());
            for (const auto& wave : waves) {
              for (size_t neighbour_idx : cell_neighbours) {
                update_cell(next_grid[neighbour_idx], wave);
              }
              update_cell(next_cell, wave);
            }
          }
          calculated_interfaces.emplace(cell_idx, curr_cell.left_idx);
        }

        // - Right outer interface -----------------------------------------------------------------
        if (curr_grid.is_cell(curr_cell.right_idx) &&
            !calculated_interfaces.contains({cell_idx, curr_cell.right_idx})) {
          const auto& other_cell = curr_grid[curr_cell.right_idx];
          const auto interfaces =
              get_shared_interfaces<ActiveFloat, PassiveFloat, DIM, GridCoord<ActiveFloat>>(
                  curr_cell, other_cell, RIGHT);

          for (const auto& interface : interfaces) {
            const auto waves =
                calculate_interface(interface, dt, curr_grid.scale_x(), curr_grid.scale_y());
            for (const auto& wave : waves) {
              for (size_t neighbour_idx : cell_neighbours) {
                update_cell(next_grid[neighbour_idx], wave);
              }
              update_cell(next_cell, wave);
            }
          }
          calculated_interfaces.emplace(cell_idx, curr_cell.right_idx);
        }

        // - Bottom outer interface ----------------------------------------------------------------
        if (curr_grid.is_cell(curr_cell.bottom_idx) &&
            !calculated_interfaces.contains({cell_idx, curr_cell.bottom_idx})) {
          const auto& other_cell = curr_grid[curr_cell.bottom_idx];
          const auto interfaces =
              get_shared_interfaces<ActiveFloat, PassiveFloat, DIM, GridCoord<ActiveFloat>>(
                  curr_cell, other_cell, BOTTOM);

          for (const auto& interface : interfaces) {
            const auto waves =
                calculate_interface(interface, dt, curr_grid.scale_x(), curr_grid.scale_y());
            for (const auto& wave : waves) {
              for (size_t neighbour_idx : cell_neighbours) {
                update_cell(next_grid[neighbour_idx], wave);
              }
              update_cell(next_cell, wave);
            }
          }
          calculated_interfaces.emplace(cell_idx, curr_cell.bottom_idx);
        }

        // - Top outer interface -------------------------------------------------------------------
        if (curr_grid.is_cell(curr_cell.top_idx) &&
            !calculated_interfaces.contains({cell_idx, curr_cell.top_idx})) {
          const auto& other_cell = curr_grid[curr_cell.top_idx];
          const auto interfaces =
              get_shared_interfaces<ActiveFloat, PassiveFloat, DIM, GridCoord<ActiveFloat>>(
                  curr_cell, other_cell, TOP);

          for (const auto& interface : interfaces) {
            const auto waves =
                calculate_interface(interface, dt, curr_grid.scale_x(), curr_grid.scale_y());
            for (const auto& wave : waves) {
              for (size_t neighbour_idx : cell_neighbours) {
                update_cell(next_grid[neighbour_idx], wave);
              }
              update_cell(next_cell, wave);
            }
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

// -------------------------------------------------------------------------------------------------
template <ExtendType extend_type, typename A, typename B>
[[nodiscard]] constexpr auto make_solver(A&& a, B&& b) noexcept -> Solver<extend_type, A, B> {
  return Solver<extend_type, A, B>(std::forward<A>(a), std::forward<B>(b));
}

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_SOLVER_HPP_
