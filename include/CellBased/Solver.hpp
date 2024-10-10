#ifndef ZAP_CELL_BASED_SOLVER_HPP_
#define ZAP_CELL_BASED_SOLVER_HPP_

#include <numbers>

#include "CellBased/Geometry.hpp"
#include "CellBased/Grid.hpp"
#include "CellBased/Interface.hpp"

#ifndef OLD_SOLVER

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

template <typename A, typename B>
class Solver {
  A m_A;
  B m_B;

  // -----------------------------------------------------------------------------------------------
  template <typename Float, size_t DIM>
  [[nodiscard]] constexpr auto cfl_factor(const UniformGrid<Float, DIM>& grid) const noexcept
      -> Float {
    return std::transform_reduce(
        std::cbegin(grid.cells()),
        std::cend(grid.cells()),
        Float{0},
        [](const Float a, const Float b) { return std::max(a, b); },
        [this](const Cell<Float, DIM>& cell) {
          assert(cell.is_cartesian() || cell.is_cut());
          if (cell.is_cartesian()) {
            return std::max(m_A.max_abs_eig_val(cell.get_cartesian().value),
                            m_B.max_abs_eig_val(cell.get_cartesian().value));
          } else {
            return std::max(std::max(m_A.max_abs_eig_val(cell.get_cut().left_value),
                                     m_A.max_abs_eig_val(cell.get_cut().right_value)),
                            std::max(m_B.max_abs_eig_val(cell.get_cut().left_value),
                                     m_B.max_abs_eig_val(cell.get_cut().right_value)));
          }
        });
  }

  // -----------------------------------------------------------------------------------------------
  template <typename Float, size_t DIM, Point2D_c PointType>
  requires(DIM > 0)
  struct WaveProperties {
    Eigen::Vector<Float, DIM> value;
    Geometry::Polygon<PointType> polygon;
  };

  // -----------------------------------------------------------------------------------------------
  template <typename Float, size_t DIM, Point2D_c PointType>
  requires(DIM > 0)
  [[nodiscard]] constexpr auto
  calculate_interface(const FullInterface<Float, DIM, PointType>& interface,
                      Float dt,
                      Float scale_x,
                      Float scale_y) const noexcept
      -> SmallVector<WaveProperties<Float, DIM, PointType>> {
    // Vector tangential to interface
    const PointType tangent_vector = (interface.end - interface.begin).normalized();

    // Angle between cut_vector and y-axis (0, 1)
    //     cut_angle = arccos(cut_vector^T * (0, 1) / ||cut_vector|| * ||(0, 1)||)
    // <=> cut_angle = arccos(cut_vector_1 / ||cut_vector||)
    auto interface_angle = std::acos(tangent_vector.y);
    if (tangent_vector.x > 0) { interface_angle = 2 * std::numbers::pi_v<Float> - interface_angle; }

    // Vector normal to cut
    PointType normal_vector{std::cos(interface_angle), std::sin(interface_angle)};
    assert(std::abs(tangent_vector.dot(normal_vector)) <= 1e-8);
    // Scale normal vector when operating in grid coordinates
    if constexpr (is_GridCoord_v<PointType>) {
      normal_vector.x *= scale_x;
      normal_vector.y *= scale_y;
    }

    const Eigen::Vector<Float, DIM> u_mid = (interface.left_value + interface.right_value) / 2;

    // Matrix for rotated PDE in normal direction to cut
    // TODO: Why cos(interface_angle) * A + sin(interface_angle) * B and not -?
    const Eigen::Matrix<Float, DIM, DIM> eta_mat =
        std::cos(interface_angle) * m_A.mat(u_mid) + std::sin(interface_angle) * m_B.mat(u_mid);

    Eigen::Matrix<Float, DIM, DIM> eig_vals;
    Eigen::Matrix<Float, DIM, DIM> eig_vecs;
    get_eigen_decomp(eta_mat, eig_vals, eig_vecs);
    const Eigen::Matrix<Float, DIM, DIM> wave_lengths = eig_vals * dt;

    // Eigen expansion of jump; wave strength
    const Eigen::Vector<Float, DIM> alpha =
        eig_vecs.inverse() * (interface.right_value - interface.left_value);

    SmallVector<WaveProperties<Float, DIM, PointType>> waves(DIM);

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
  template <typename Float, size_t DIM, Point2D_c PointType>
  constexpr void update_cell(Cell<Float, DIM>& cell,
                             const WaveProperties<Float, DIM, PointType>& wave) {
    // TODO: Do something about periodic boundary conditions
    assert(cell.is_cartesian() || cell.is_cut());
    auto update_value = [&](Eigen::Vector<Float, DIM>& value,
                            const Geometry::Polygon<PointType>& cell_polygon) {
      const auto cell_area = cell_polygon.area();
      assert(cell_area > 0 || std::abs(cell_area) <= EPS<Float>);
      if (std::abs(cell_area) <= EPS<Float>) { return; }

      const auto intersect      = Geometry::intersection(cell_polygon, wave.polygon);
      const auto intersect_area = intersect.area();
      assert(intersect_area >= 0 || std::abs(intersect_area) < 1e-8);
      assert(intersect_area <= cell_area &&
             "Something went wrong in the calculation of the intersection.");

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
  template <typename Float, size_t DIM>
  constexpr void apply_side_interfaces(const UniformGrid<Float, DIM>& curr_grid,
                                       Cell<Float, DIM>& next_cell,
                                       const Cell<Float, DIM>& curr_cell,
                                       size_t idx,
                                       Side side,
                                       Float dt) {
    if (curr_grid.is_cell(idx)) {
      const auto& other_cell = curr_grid[idx];
      const auto interfaces =
          get_shared_interfaces<Float, DIM, GridCoord<Float>>(curr_cell, other_cell, side);
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
  template <typename Float, size_t DIM, Point2D_c PointType>
  [[nodiscard]] constexpr auto get_internal_interface(const Cell<Float, DIM>& cell) const noexcept
      -> FullInterface<Float, DIM, PointType> {
    assert(cell.is_cut());

    constexpr CoordType coord_type = PointType2CoordType<PointType>;

    return FullInterface<Float, DIM, PointType>{
        .left_value  = cell.get_cut().left_value,
        .right_value = cell.get_cut().right_value,
        .begin       = cell.template cut1<coord_type>(),
        .end         = cell.template cut2<coord_type>(),
    };
  }

  // -----------------------------------------------------------------------------------------------
  template <typename Float, size_t DIM>
  [[nodiscard]] constexpr auto
  move_wave_front(const UniformGrid<Float, DIM>& curr_grid,
                  [[maybe_unused]] const UniformGrid<Float, DIM>& next_grid,
                  Float dt) const noexcept -> std::optional<std::vector<GridCoord<Float>>> {
    using PointType = GridCoord<Float>;

    // Move old cuts according to strongest wave
    std::vector<std::pair<PointType, PointType>> new_shock_points;
    for (size_t cell_idx : curr_grid.cut_cell_idxs()) {
      const auto& curr_cell = curr_grid[cell_idx];
      assert(curr_cell.is_cut());

      const auto interface = get_internal_interface<Float, DIM, PointType>(curr_cell);

      const auto tangent_vector = (interface.end - interface.begin).normalized();
      assert(std::abs(tangent_vector.norm() - 1) <= 1e-8);

      auto interface_angle = std::acos(tangent_vector.y);
      if (tangent_vector.x > 0) {
        interface_angle = 2 * std::numbers::pi_v<Float> - interface_angle;
      }
      PointType normal_vector{.x = std::cos(interface_angle), .y = std::sin(interface_angle)};
      assert(std::abs(tangent_vector.dot(normal_vector)) <= 1e-8);
      // Scale normal vector when operating in grid coordinates
      if constexpr (is_GridCoord_v<PointType>) {
        normal_vector.x *= curr_grid.scale_x();
        normal_vector.y *= curr_grid.scale_y();
      }

      const Eigen::Vector<Float, DIM> u_mid = (interface.left_value + interface.right_value) / 2;

      const Eigen::Matrix<Float, DIM, DIM> eta_mat =
          std::cos(interface_angle) * m_A.mat(u_mid) + std::sin(interface_angle) * m_B.mat(u_mid);

      Eigen::Matrix<Float, DIM, DIM> eig_vals;
      Eigen::Matrix<Float, DIM, DIM> eig_vecs;
      get_eigen_decomp(eta_mat, eig_vals, eig_vecs);

      // Wave strength
      const Eigen::Vector<Float, DIM> alpha =
          eig_vecs.inverse() * (interface.right_value - interface.left_value);

      Eigen::Index strong_wave_idx = 0;
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        // TODO: Abs or not?
        if (std::abs(alpha(i)) > std::abs(alpha(strong_wave_idx))) { strong_wave_idx = i; }
      }

      const auto wave_length = eig_vals(strong_wave_idx) * dt;

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
  template <typename Float, size_t DIM, typename GridWriter, typename TimeWriter>
  [[nodiscard]] auto solve(UniformGrid<Float, DIM> grid,
                           Float tend,
                           GridWriter& grid_writer,
                           TimeWriter& time_writer,
                           Float CFL_safety_factor = 0.5) noexcept
      -> std::optional<UniformGrid<Float, DIM>> {
    if (!(CFL_safety_factor > 0 && CFL_safety_factor <= 1)) {
      Igor::Warn("CFL_safety_factor must be in (0, 1], is {}", CFL_safety_factor);
      return std::nullopt;
    }
    auto& curr_grid = grid;
    auto next_grid  = grid;

    if (!grid_writer.write_data(curr_grid) || !time_writer.write_data(Float{0})) {
      return std::nullopt;
    }

    for (Float t = 0.0; t < tend;) {
      const Float CFL_factor = cfl_factor(curr_grid);

      if (std::isnan(CFL_factor) || std::isinf(CFL_factor)) {
        Igor::Warn("CFL_factor is invalid at time t={}: CFL_factor = {}", t, CFL_factor);
        return std::nullopt;
      }
      const Float dt = std::min(CFL_safety_factor * curr_grid.min_delta() / CFL_factor, tend - t);

      next_grid = curr_grid;
#ifndef ZAP_STATIC_CUT
      if (!curr_grid.cut_cell_idxs().empty()) {
        next_grid.merge_cut_cells();

        auto avg_new_shock_points = move_wave_front(curr_grid, next_grid, dt);
        if (!avg_new_shock_points.has_value()) { return std::nullopt; }

        if (!next_grid.cut_piecewise_linear(std::move(*avg_new_shock_points))) {
          Igor::Warn("Could not cut on new shock curve.");
          return std::nullopt;
        }
        // Igor::Debug("t = {}", t);
        // Igor::Debug("#cut cells curr_grid = {}", curr_grid.cut_cell_idxs().size());
        // Igor::Debug("#cut cells next_grid = {}", next_grid.cut_cell_idxs().size());
        // std::cout << "----------------------------------------\n";
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
                   std::abs(next_subcell_polygon.area()) <= EPS<Float>);
            if (std::abs(next_subcell_polygon.area()) > EPS<Float>) {
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
                   std::abs(next_subcell_polygon.area()) <= EPS<Float>);
            if (std::abs(next_subcell_polygon.area()) > EPS<Float>) {
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
#endif  // ZAP_STATIC_CUT

      // TODO: What happens when a subcell has area 0?
      //   -> cannot "uncut" the cell because that would loose shock position information

#ifndef USE_FLUX_FOR_CARTESIAN
#pragma omp parallel for
      for (size_t cell_idx = 0; cell_idx < curr_grid.size(); ++cell_idx) {
        const auto& curr_cell = curr_grid[cell_idx];
        auto& next_cell       = next_grid[cell_idx];

        // Left interface
        apply_side_interfaces(curr_grid, next_cell, curr_cell, curr_cell.left_idx, LEFT, dt);
        // Right interface
        apply_side_interfaces(curr_grid, next_cell, curr_cell, curr_cell.right_idx, RIGHT, dt);
        // Bottom interface
        apply_side_interfaces(curr_grid, next_cell, curr_cell, curr_cell.bottom_idx, BOTTOM, dt);
        // Top interface
        apply_side_interfaces(curr_grid, next_cell, curr_cell, curr_cell.top_idx, TOP, dt);
      }

      for (size_t cell_idx : curr_grid.cut_cell_idxs()) {
        const auto& curr_cell = curr_grid[cell_idx];
        assert(curr_cell.is_cut());
        auto& next_cell = next_grid[cell_idx];
        assert(next_cell.is_cartesian() || next_cell.is_cut());

        // - Handle internal interface -----------------------------------------------------------
        const auto internal_interface =
            get_internal_interface<Float, DIM, GridCoord<Float>>(curr_cell);
        const auto internal_waves = calculate_interface<Float, DIM, GridCoord<Float>>(
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
#else
      for (size_t cell_idx = 0; cell_idx < curr_grid.size(); ++cell_idx) {
        const auto& curr_cell = curr_grid[cell_idx];
        assert(curr_cell.is_cartesian() || curr_cell.is_cut());
        auto& next_cell = next_grid[cell_idx];
        assert(next_cell.is_cartesian() || next_cell.is_cut());

        // curr_cell must be cartesian, cut cells are handled differently
        if (curr_cell.is_cartesian() && next_cell.is_cartesian()) {
          auto calc_flux = [this, &curr_grid, dt](const Cell<Float, DIM>& cell,
                                                  size_t idx,
                                                  Side from) -> Eigen::Vector<Float, DIM> {
            if (curr_grid.is_cell(idx)) {
              return linearized_xy_flux(cell, curr_grid[idx], from, dt);
            } else if (idx == SAME_VALUE_INDEX) {
              const auto factor = (from & (LEFT | RIGHT)) > 0 ? (dt / cell.dx)  //
                                                              : (dt / cell.dy);
              return linearized_xy_flux_impl(
                  cell.get_cartesian().value, cell.get_cartesian().value, factor, from);
            } else if (idx == ZERO_FLUX_INDEX) {
              return Eigen::Vector<Float, DIM>::Zero();
            } else if (idx == NULL_INDEX) {
              std::stringstream s;
              s << cell;
              Igor::Debug("cell = {}", s.str());
              Igor::Panic("cell has NULL index.");
            } else {
              std::stringstream s;
              s << cell;
              Igor::Debug("cell = {}", s.str());
              Igor::Panic("cell has unknown index with value {}.", idx);
            }
            std::unreachable();
          };

          // Left side: U_(i, j-1) -> U_(i, j)
          const auto F_minus = calc_flux(curr_cell, curr_cell.left_idx, LEFT);

          // Right side: U_(i, j) -> U_(i, j+1)
          const auto F_plus = calc_flux(curr_cell, curr_cell.right_idx, RIGHT);

          // Bottom side: U_(i-1, j) -> U_(i, j)
          const auto G_minus = calc_flux(curr_cell, curr_cell.bottom_idx, BOTTOM);

          // Top side: U_(i, j) -> U_(i+1, j)
          const auto G_plus = calc_flux(curr_cell, curr_cell.top_idx, TOP);

          next_cell.get_cartesian().value = curr_cell.get_cartesian().value -  //
                                            (F_plus - F_minus) -               //
                                            (G_plus - G_minus);
        } else {
          // - Handle internal interface -----------------------------------------------------------
          {
            const auto internal_interface = get_internal_interface(curr_cell);
            const auto internal_waves     = calculate_interface(internal_interface, dt);

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
          // - Handle internal interface -----------------------------------------------------------

          // Left interface
          apply_side_interfaces(curr_grid, next_cell, curr_cell, curr_cell.left_idx, LEFT, dt);
          // Right interface
          apply_side_interfaces(curr_grid, next_cell, curr_cell, curr_cell.right_idx, RIGHT, dt);
          // Bottom interface
          apply_side_interfaces(curr_grid, next_cell, curr_cell, curr_cell.bottom_idx, BOTTOM, dt);
          // Top interface
          apply_side_interfaces(curr_grid, next_cell, curr_cell, curr_cell.top_idx, TOP, dt);
        }
      }
#endif

      // Update time
      t += dt;

      curr_grid = next_grid;

      if (!grid_writer.write_data(curr_grid) || !time_writer.write_data(t)) { return std::nullopt; }
      // break;
    }

    return curr_grid;
  }

#ifdef USE_FLUX_FOR_CARTESIAN
  // -----------------------------------------------------------------------------------------------
  template <typename Float, int DIM>
  requires(DIM > 0)
  [[nodiscard]] constexpr auto linearized_xy_flux_impl(const Eigen::Vector<Float, DIM>& curr_value,
                                                       const Eigen::Vector<Float, DIM>& other_value,
                                                       Float factor,
                                                       Side from) -> Eigen::Vector<Float, DIM> {
    const Eigen::Vector<Float, DIM> u_mid   = (curr_value + other_value) / 2;
    Eigen::Matrix<Float, DIM, DIM> eig_vals = (from & (LEFT | RIGHT)) > 0 ? m_A.eig_vals(u_mid)  //
                                                                          : m_B.eig_vals(u_mid);
    Eigen::Matrix<Float, DIM, DIM> eig_vecs = (from & (LEFT | RIGHT)) > 0 ? m_A.eig_vecs(u_mid)  //
                                                                          : m_B.eig_vecs(u_mid);

    bool all_zero = true;
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
      if (std::abs(eig_vals(i, i)) >= 1e-6) { all_zero = false; }
    }
    if (all_zero) { return Eigen::Vector<Float, DIM>::Zero(); }

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

  // -----------------------------------------------------------------------------------------------
  template <typename Float, size_t DIM>
  requires(DIM > 0)
  [[nodiscard]] constexpr auto linearized_xy_flux(const Cell<Float, DIM>& curr_cell,
                                                  const Cell<Float, DIM>& other_cell,
                                                  Side from,
                                                  Float dt) -> Eigen::Vector<Float, DIM> {
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
                  const auto dx_left_half     = other_cell.get_cut().x1_cut - other_cell.x_min;
                  const auto factor_left_half = (dx_left_half / other_cell.dx) * dt / other_cell.dy;
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
          Igor::Panic("Unknown cut type with value {}",
                      static_cast<std::underlying_type_t<CutType>>(other_cell.get_cut().type));
          std::unreachable();
      }

      Igor::Panic("Unreachable.");
      std::unreachable();
    }
  }
#endif  // USE_FLUX_FOR_CARTESIAN
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
  [[nodiscard]] constexpr auto apply_flux(const UniformGrid<Float, DIM>& grid,
                                          const auto& cell,
                                          Side side) const noexcept -> Float {
    size_t idx = NULL_INDEX;
    switch (side) {
      case BOTTOM: idx = cell.bottom_idx; break;
      case RIGHT:  idx = cell.right_idx; break;
      case TOP:    idx = cell.top_idx; break;
      case LEFT:   idx = cell.left_idx; break;
      default:     Igor::Panic("Must be BOTTOM, RIGHT, TOP, or LEFT."); std::unreachable();
    }

    if (!grid.is_cell(idx)) { return apply_boundary_flux(cell, side, idx); }
    const auto& other_cell = grid[idx];
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
      s << grid[idx];
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
      s << grid[idx];
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
      s << grid[idx];
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
  [[nodiscard]] constexpr auto apply_flux(const UniformGrid<Float, DIM>& grid,
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
  auto solve(UniformGrid<Float, DIM> grid,
             Float tend,
             GridWriter& grid_writer,
             TimeWriter& time_writer) noexcept -> std::optional<UniformGrid<Float, DIM>> {
    auto& curr_grid = grid;
    auto next_grid  = grid;

    if (!grid_writer.write_data(curr_grid) || !time_writer.write_data(Float{0})) {
      return std::nullopt;
    }

    for (Float t = 0.0; t < tend;) {
      const Float CFL_factor = curr_grid.abs_max_value();
      if (std::isnan(CFL_factor) || std::isinf(CFL_factor)) {
        Igor::Warn("CFL_factor is invalid at time t={}: CFL_factor = {}", t, CFL_factor);
        return std::nullopt;
      }
      const Float dt = std::min(0.5 * curr_grid.min_delta() / CFL_factor, tend - t);

#pragma omp parallel for
      for (std::size_t i = 0; i < curr_grid.size(); ++i) {
        const auto& curr_cell = curr_grid[i];
        auto& next_cell       = next_grid[i];

        if constexpr (DIM == 1) {
          if (curr_cell.is_cartesian()) {
            const auto F_minus = apply_flux(curr_grid, curr_cell, LEFT);
            const auto F_plus  = apply_flux(curr_grid, curr_cell, RIGHT);
            const auto G_minus = apply_flux(curr_grid, curr_cell, BOTTOM);
            const auto G_plus  = apply_flux(curr_grid, curr_cell, TOP);

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

      std::swap(curr_grid.cells(), next_grid.cells());

      if (!grid_writer.write_data(curr_grid) || !time_writer.write_data(t)) { return std::nullopt; }
    }

    return curr_grid;
  }
};

}  // namespace Zap::CellBased

#endif  // OLD_SOLVER

#endif  // ZAP_CELL_BASED_SOLVER_HPP_
