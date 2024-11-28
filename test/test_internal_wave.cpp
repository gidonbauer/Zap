#include <gtest/gtest.h>

#include <numbers>

#include "CellBased/Solver.hpp"

#include "Igor/Logging.hpp"

using namespace Zap::CellBased;
using Float           = double;
constexpr Float X_MIN = 0.0;
constexpr Float X_MAX = 5.0;
constexpr Float Y_MIN = 0.0;
constexpr Float Y_MAX = 5.0;
#if 1

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto sqr(auto x) noexcept { return x * x; }

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto u0(Float x, Float y) noexcept -> Float {
  return (sqr(x - X_MIN) + sqr(y - Y_MIN)) *
         static_cast<Float>((sqr(x - X_MIN) + sqr(y - Y_MIN)) <=
                            sqr((X_MIN + X_MAX + Y_MIN + Y_MAX) / 4));
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto init_shock(Float t) noexcept -> Zap::CellBased::SimCoord<Float> {
  const auto r = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;
  return {
      r * std::cos(std::numbers::pi_v<Float> / 2 * t),
      r * std::sin(std::numbers::pi_v<Float> / 2 * t),
  };
}

// -------------------------------------------------------------------------------------------------
TEST(InternalWave, QuarterCircle) {
  UniformGrid<Float, Float> grid(X_MIN, X_MAX, 5UZ, Y_MIN, Y_MAX, 5UZ);
  grid.periodic_boundary();
  const bool cut_success = grid.cut_curve(init_shock);
  ASSERT_TRUE(cut_success);
  grid.fill_four_point(u0);

  {
    Float avg_left_value = 0.0;
    for (size_t cell_idx : grid.cut_cell_idxs()) {
      avg_left_value += grid[cell_idx].get_cut().left_value;
    }
    avg_left_value /= static_cast<Float>(grid.cut_cell_idxs().size());

    for (size_t cell_idx : grid.cut_cell_idxs()) {
      grid[cell_idx].get_cut().left_value = avg_left_value;
    }
  }

  const auto dx = grid[0].template dx<SIM_C>();
  const auto dy = grid[0].template dy<SIM_C>();
  const auto dt = 0.1;

  Igor::Info("cut cells = {}", grid.cut_cell_idxs());

  std::vector<FullInterface<Float, GridCoord<Float>>> internal_interfaces(
      grid.cut_cell_idxs().size());
  for (size_t i = 0; i < grid.cut_cell_idxs().size(); ++i) {
    const auto interface =
        get_internal_interface<Float, Float, GridCoord<Float>>(grid[grid.cut_cell_idxs()[i]]);
    internal_interfaces[i] = interface;
  }

  std::cout << "------------------------------------------------------------\n";
  Igor::Info("internal_interfaces = {}", internal_interfaces);

  std::vector<FreeWave<Float, GridCoord<Float>>> internal_waves(internal_interfaces.size());
  for (size_t i = 0; i < internal_interfaces.size(); ++i) {
    const auto wave = calc_wave<FREE>(internal_interfaces[i], dx, dy, dt);
    ASSERT_TRUE(wave.has_value());
    internal_waves[i] = *wave;
  }

  std::cout << "------------------------------------------------------------\n";
  Igor::Info("internal_waves = {}", internal_waves);

  std::vector<GridCoord<Float>> avg_new_shock_points;
  move_wave_front(avg_new_shock_points, internal_waves, dt);

  std::cout << "------------------------------------------------------------\n";
  Igor::Info("avg_new_shock_points = {}", avg_new_shock_points);
}

#else

TEST(InternalWave, A) {
  UniformGrid<Float, Float> grid(0.0, 1.0, 1, 0.0, 1.0, 1);
  grid.periodic_boundary();
  grid[0].cell_type = CutValue<Float>{
      .left_value    = 2.0,
      .right_value   = 0.0,
      .rel_cut_entry = {1.0, 0.5},
      .entry_loc     = RIGHT,
      .rel_cut_exit  = {0.0, 0.5},
      .exit_loc      = LEFT,

  };

  const auto internal_interface = get_internal_interface<Float, Float, GridCoord<Float>>(grid[0]);
  Igor::Info("internal_interface = {{ .left_value = {}, .right_value = {}, .begin = {}, .end = {} "
             "}}",
             internal_interface.left_value,
             internal_interface.right_value,
             internal_interface.begin,
             internal_interface.end);

  const Float dt           = 0.1;
  const auto internal_wave = calc_wave<FREE>(
      internal_interface, grid[0].template dx<SIM_C>(), grid[0].template dy<SIM_C>(), dt);
  ASSERT_TRUE(internal_wave.has_value());
  Igor::Info("internal_wave = {}", *internal_wave);
}

TEST(InternalWave, B) {
  UniformGrid<Float, Float> grid(0.0, 1.0, 1, 0.0, 1.0, 1);
  grid.periodic_boundary();
  grid[0].cell_type = CutValue<Float>{
      .left_value    = 2.0,
      .right_value   = 0.0,
      .rel_cut_entry = {0.5, 1.0},
      .entry_loc     = TOP,
      .rel_cut_exit  = {0.0, 0.5},
      .exit_loc      = LEFT,

  };

  const auto internal_interface = get_internal_interface<Float, Float, GridCoord<Float>>(grid[0]);
  Igor::Info("internal_interface = {{ .left_value = {}, .right_value = {}, .begin = {}, .end = {} "
             "}}",
             internal_interface.left_value,
             internal_interface.right_value,
             internal_interface.begin,
             internal_interface.end);

  const Float dt           = 0.1;
  const auto internal_wave = calc_wave<FREE>(
      internal_interface, grid[0].template dx<SIM_C>(), grid[0].template dy<SIM_C>(), dt);
  ASSERT_TRUE(internal_wave.has_value());
  Igor::Info("internal_wave = {}", *internal_wave);
}

#endif
