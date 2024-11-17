#include <gtest/gtest.h>

#include "CellBased/Grid.hpp"
#include "CellBased/Interface.hpp"
#include "CellBased/Wave.hpp"

using namespace Zap::CellBased;

using Float = double;

// -------------------------------------------------------------------------------------------------
TEST(PeriodicBoundary, CartesianBottom) {
  UniformGrid<Float, Float> grid(0.0, 1.0, 2, 0.0, 1.0, 2);
  grid.periodic_boundary();

  for (size_t i = 0; i < grid.size(); ++i) {
    grid[i].get_cartesian().value = static_cast<Float>(i);
  }

  // - 0 to 2 => bottom periodic boundary ----------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[0], grid[2], BOTTOM);

    ASSERT_EQ(interfaces.size(), 1UZ);
    const auto interface = interfaces[0];

    EXPECT_DOUBLE_EQ(interface.left_value, 2.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 0.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{1.0, 0.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{0.0, 0.0}).norm(), EPS<Float>);
  }

  // - 1 to 3 => bottom periodic boundary ----------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[1], grid[3], BOTTOM);

    ASSERT_EQ(interfaces.size(), 1UZ);
    const auto interface = interfaces[0];

    EXPECT_DOUBLE_EQ(interface.left_value, 3.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 1.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{2.0, 0.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{1.0, 0.0}).norm(), EPS<Float>);
  }
}

// -------------------------------------------------------------------------------------------------
TEST(PeriodicBoundary, CartesianTop) {
  UniformGrid<Float, Float> grid(0.0, 1.0, 2, 0.0, 1.0, 2);
  grid.periodic_boundary();

  for (size_t i = 0; i < grid.size(); ++i) {
    grid[i].get_cartesian().value = static_cast<Float>(i);
  }

  // - 2 to 0 => top periodic boundary -------------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[2], grid[0], TOP);

    ASSERT_EQ(interfaces.size(), 1UZ);
    const auto interface = interfaces[0];

    EXPECT_DOUBLE_EQ(interface.left_value, 2.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 0.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{1.0, 2.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{0.0, 2.0}).norm(), EPS<Float>);
  }

  // - 3 to 1 => top periodic boundary -------------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[3], grid[1], TOP);

    ASSERT_EQ(interfaces.size(), 1UZ);
    const auto interface = interfaces[0];

    EXPECT_DOUBLE_EQ(interface.left_value, 3.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 1.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{2.0, 2.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{1.0, 2.0}).norm(), EPS<Float>);
  }
}

// -------------------------------------------------------------------------------------------------
TEST(PeriodicBoundary, CartesianLeft) {
  UniformGrid<Float, Float> grid(0.0, 1.0, 2, 0.0, 1.0, 2);
  grid.periodic_boundary();

  for (size_t i = 0; i < grid.size(); ++i) {
    grid[i].get_cartesian().value = static_cast<Float>(i);
  }

  // - 0 to 1 => left periodic boundary ------------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[0], grid[1], LEFT);

    ASSERT_EQ(interfaces.size(), 1UZ);
    const auto interface = interfaces[0];

    EXPECT_DOUBLE_EQ(interface.left_value, 1.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 0.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{0.0, 0.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{0.0, 1.0}).norm(), EPS<Float>);
  }

  // - 2 to 3 => left periodic boundary ------------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[2], grid[3], LEFT);

    ASSERT_EQ(interfaces.size(), 1UZ);
    const auto interface = interfaces[0];

    EXPECT_DOUBLE_EQ(interface.left_value, 3.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 2.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{0.0, 1.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{0.0, 2.0}).norm(), EPS<Float>);
  }
}

// -------------------------------------------------------------------------------------------------
TEST(PeriodicBoundary, CartesianRight) {
  UniformGrid<Float, Float> grid(0.0, 1.0, 2, 0.0, 1.0, 2);
  grid.periodic_boundary();

  for (size_t i = 0; i < grid.size(); ++i) {
    grid[i].get_cartesian().value = static_cast<Float>(i);
  }

  // - 1 to 0 => right periodic boundary -----------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[1], grid[0], RIGHT);

    ASSERT_EQ(interfaces.size(), 1UZ);
    const auto interface = interfaces[0];

    EXPECT_DOUBLE_EQ(interface.left_value, 1.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 0.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{2.0, 0.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{2.0, 1.0}).norm(), EPS<Float>);
  }

  // - 3 to 2 => right periodic boundary -----------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[3], grid[2], RIGHT);

    ASSERT_EQ(interfaces.size(), 1UZ);
    const auto interface = interfaces[0];

    EXPECT_DOUBLE_EQ(interface.left_value, 3.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 2.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{2.0, 1.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{2.0, 2.0}).norm(), EPS<Float>);
  }
}

// -------------------------------------------------------------------------------------------------
TEST(PeriodicBoundary, CutLeft) {
  UniformGrid<Float, Float> grid(0.0, 1.0, 2, 0.0, 1.0, 2);
  grid.periodic_boundary();

  for (size_t i = 0; i < grid.size(); ++i) {
    grid[i].get_cartesian().value = static_cast<Float>(i);
  }
  grid[0].cell_type = CutValue<Float>{
      .left_value    = 0.25,
      .right_value   = 0.75,
      .rel_cut_entry = {1.0, 0.5},
      .entry_loc     = RIGHT,
      .rel_cut_exit  = {0.0, 0.5},
      .exit_loc      = LEFT,
  };

  // - 0 to 1 => left periodic boundary ------------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[0], grid[1], LEFT);

    ASSERT_EQ(interfaces.size(), 2UZ);

    auto interface = interfaces[0];
    EXPECT_DOUBLE_EQ(interface.left_value, 1.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 0.25);
    EXPECT_LE((interface.begin - GridCoord<Float>{0.0, 0.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{0.0, 0.5}).norm(), EPS<Float>);

    interface = interfaces[1];
    EXPECT_DOUBLE_EQ(interface.left_value, 1.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 0.75);
    EXPECT_LE((interface.begin - GridCoord<Float>{0.0, 0.5}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{0.0, 1.0}).norm(), EPS<Float>);
  }

  // - 1 to 0 => right periodic boundary -----------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[1], grid[0], RIGHT);

    ASSERT_EQ(interfaces.size(), 2UZ);

    auto interface = interfaces[0];
    EXPECT_DOUBLE_EQ(interface.left_value, 1.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 0.25);
    EXPECT_LE((interface.begin - GridCoord<Float>{2.0, 0.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{2.0, 0.5}).norm(), EPS<Float>);

    interface = interfaces[1];
    EXPECT_DOUBLE_EQ(interface.left_value, 1.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 0.75);
    EXPECT_LE((interface.begin - GridCoord<Float>{2.0, 0.5}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{2.0, 1.0}).norm(), EPS<Float>);
  }

  // - 0 to 1 => right regular boundary ------------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[0], grid[1], RIGHT);

    ASSERT_EQ(interfaces.size(), 2UZ);

    auto interface = interfaces[0];
    EXPECT_DOUBLE_EQ(interface.left_value, 0.25);
    EXPECT_DOUBLE_EQ(interface.right_value, 1.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{1.0, 0.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{1.0, 0.5}).norm(), EPS<Float>);

    interface = interfaces[1];
    EXPECT_DOUBLE_EQ(interface.left_value, 0.75);
    EXPECT_DOUBLE_EQ(interface.right_value, 1.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{1.0, 0.5}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{1.0, 1.0}).norm(), EPS<Float>);
  }
}

// -------------------------------------------------------------------------------------------------
TEST(PeriodicBoundary, CutTop) {
  UniformGrid<Float, Float> grid(0.0, 1.0, 2, 0.0, 1.0, 2);
  grid.periodic_boundary();

  for (size_t i = 0; i < grid.size(); ++i) {
    grid[i].get_cartesian().value = static_cast<Float>(i);
  }
  grid[3].cell_type = CutValue<Float>{
      .left_value    = 3.25,
      .right_value   = 3.75,
      .rel_cut_entry = {0.5, 0.0},
      .entry_loc     = BOTTOM,
      .rel_cut_exit  = {0.5, 1.0},
      .exit_loc      = TOP,
  };

  // - 3 to 1 => top periodic boundary -------------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[3], grid[1], TOP);

    ASSERT_EQ(interfaces.size(), 2UZ);

    auto interface = interfaces[0];
    EXPECT_DOUBLE_EQ(interface.left_value, 3.25);
    EXPECT_DOUBLE_EQ(interface.right_value, 1.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{1.5, 2.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{1.0, 2.0}).norm(), EPS<Float>);

    interface = interfaces[1];
    EXPECT_DOUBLE_EQ(interface.left_value, 3.75);
    EXPECT_DOUBLE_EQ(interface.right_value, 1.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{2.0, 2.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{1.5, 2.0}).norm(), EPS<Float>);
  }

  // - 1 to 3 => bottom periodic boundary ----------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[1], grid[3], BOTTOM);

    ASSERT_EQ(interfaces.size(), 2UZ);

    auto interface = interfaces[0];
    EXPECT_DOUBLE_EQ(interface.left_value, 3.25);
    EXPECT_DOUBLE_EQ(interface.right_value, 1.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{1.5, 0.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{1.0, 0.0}).norm(), EPS<Float>);

    interface = interfaces[1];
    EXPECT_DOUBLE_EQ(interface.left_value, 3.75);
    EXPECT_DOUBLE_EQ(interface.right_value, 1.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{2.0, 0.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{1.5, 0.0}).norm(), EPS<Float>);
  }

  // - 1 to 3 => top regular boundary --------------------------------------------------------------
  {
    const auto interfaces =
        get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[1], grid[3], TOP);

    ASSERT_EQ(interfaces.size(), 2UZ);

    auto interface = interfaces[0];
    EXPECT_DOUBLE_EQ(interface.left_value, 1.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 3.25);
    EXPECT_LE((interface.begin - GridCoord<Float>{1.5, 1.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{1.0, 1.0}).norm(), EPS<Float>);

    interface = interfaces[1];
    EXPECT_DOUBLE_EQ(interface.left_value, 1.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 3.75);
    EXPECT_LE((interface.begin - GridCoord<Float>{2.0, 1.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{1.5, 1.0}).norm(), EPS<Float>);
  }
}

TEST(PeriodicBoundary, Wave) {
  UniformGrid<Float, Float> grid(0.0, 1.0, 2, 0.0, 1.0, 2);
  grid.periodic_boundary();
  const Float dt = 0.1;
  const Float dx = grid[0].template dx<SIM_C>();
  const Float dy = grid[0].template dy<SIM_C>();

  grid[0].get_cartesian().value = 1.0;
  grid[1].get_cartesian().value = 2.0;
  grid[2].get_cartesian().value = 2.0;
  grid[3].cell_type             = CutValue<Float>{
                  .left_value    = 2.5,
                  .right_value   = 0.0,
                  .rel_cut_entry = {1.0, 0.25},
                  .entry_loc     = RIGHT,
                  .rel_cut_exit  = {0.25, 1.0},
                  .exit_loc      = TOP,
  };
  ASSERT_EQ(grid.size(), 4UZ);

  EXPECT_EQ(grid.on_boundary(0), LEFT | BOTTOM);
  EXPECT_EQ(grid.on_boundary(1), RIGHT | BOTTOM);
  EXPECT_EQ(grid.on_boundary(2), LEFT | TOP);
  EXPECT_EQ(grid.on_boundary(3), RIGHT | TOP);

  const auto interfaces =
      get_shared_interfaces<Float, Float, GridCoord<Float>>(grid[3], grid[2], RIGHT);
  ASSERT_EQ(interfaces.size(), 2UZ);
  {
    auto interface = interfaces[0];
    EXPECT_DOUBLE_EQ(interface.left_value, 2.5);
    EXPECT_DOUBLE_EQ(interface.right_value, 2.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{2.0, 1.0}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{2.0, 1.25}).norm(), EPS<Float>);

    interface = interfaces[1];
    EXPECT_DOUBLE_EQ(interface.left_value, 0.0);
    EXPECT_DOUBLE_EQ(interface.right_value, 2.0);
    EXPECT_LE((interface.begin - GridCoord<Float>{2.0, 1.25}).norm(), EPS<Float>);
    EXPECT_LE((interface.end - GridCoord<Float>{2.0, 2.0}).norm(), EPS<Float>);
  }

  SmallVector<AxisAlignedWave<Float, GridCoord<Float>, X>> waves(interfaces.size());
  for (size_t i = 0; i < interfaces.size(); ++i) {
    const auto wave = calc_wave<X>(interfaces[i], dx, dy, dt);
    ASSERT_TRUE(wave.has_value());
    waves[i] = *wave;
  }

  SmallVector<Geometry::Polygon<GridCoord<Float>>> wave_polygons(waves.size());
  for (size_t i = 0; i < waves.size(); ++i) {
    wave_polygons[i] = calc_wave_polygon(waves[i], dx, dy, dt);
  }

  EXPECT_DOUBLE_EQ((wave_polygons[0] & grid[2].get_cartesian_polygon<GRID_C>()).area(), 0.0);
  EXPECT_DOUBLE_EQ((wave_polygons[0] & grid[0].get_cartesian_polygon<GRID_C>()).area(), 0.0);

  EXPECT_DOUBLE_EQ((wave_polygons[1] & grid[2].get_cartesian_polygon<GRID_C>()).area(), 0.0);
  EXPECT_DOUBLE_EQ((wave_polygons[1] & grid[0].get_cartesian_polygon<GRID_C>()).area(), 0.0);

  EXPECT_GT((wave_polygons[0] &
             Geometry::Polygon(grid.translate(grid[2].get_cartesian_points<GRID_C>(), RIGHT)))
                .area(),
            0.0);
  EXPECT_GE(
      (wave_polygons[0] & Geometry::Polygon(grid.translate(grid[0].get_cartesian_points<GRID_C>(),
                                                           static_cast<Side>(RIGHT | TOP))))
          .area(),
      0.0);

  EXPECT_GT((wave_polygons[1] &
             Geometry::Polygon(grid.translate(grid[2].get_cartesian_points<GRID_C>(), RIGHT)))
                .area(),
            0.0);
  EXPECT_GT(
      (wave_polygons[1] & Geometry::Polygon(grid.translate(grid[0].get_cartesian_points<GRID_C>(),
                                                           static_cast<Side>(RIGHT | TOP))))
          .area(),
      0.0);
}
