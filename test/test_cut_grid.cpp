#include <gtest/gtest.h>

#include <numbers>

#include <CellBased/Grid.hpp>

using namespace Zap::CellBased;

// =================================================================================================
TEST(CutGrid, CurveQuarterCircle) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, Float> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  auto quarter_circle = [=]<typename T>(T t) -> Zap::CellBased::SimCoord<T> {
    assert(t >= 0 && t <= 1);
    const auto r = (x_min + x_max + y_min + y_max) / 4;
    return {
        r * std::cos(std::numbers::pi_v<Float> / 2 * t),
        r * std::sin(std::numbers::pi_v<Float> / 2 * t),
    };
  };

  const auto success = grid.cut_curve(quarter_circle);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 5UZ) << "Expected exactly five cut cells.";

  EXPECT_EQ(cut_idxs[0], 2UZ);
  EXPECT_EQ(cut_idxs[1], 7UZ);
  EXPECT_EQ(cut_idxs[2], 6UZ);
  EXPECT_EQ(cut_idxs[3], 11UZ);
  EXPECT_EQ(cut_idxs[4], 10UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[0]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[1]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[2]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[3]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[4]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().exit_loc, LEFT);

  for (size_t i = 0; i < cut_idxs.size() - 1; ++i) {
    const auto p1 = grid[cut_idxs[i]].cut_exit<SIM_C>();
    const auto p2 = grid[cut_idxs[i + 1]].cut_entry<SIM_C>();
    EXPECT_DOUBLE_EQ((p1 - p2).norm(), 0.0);
  }
}

// =================================================================================================
TEST(CutGrid, PiecewiseLinearEmpty) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, Float> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points;
  const auto success = grid.cut_piecewise_linear<ExtendType::MAX>(points);
  ASSERT_FALSE(success) << "Should not be able to cut a grid with empty points.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 0UZ) << "Expected exactly zero cut cells.";
}

// =================================================================================================
TEST(CutGrid, PiecewiseLinearNoExtendInsideSingleCell) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, Float> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<GridCoord<Float>> points = {{2.5, 2.5}, {2.7, 2.7}};
  const auto success                   = grid.cut_piecewise_linear<ExtendType::NONE>(points);
  ASSERT_TRUE(success) << "Could not cut the grid.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 0UZ) << "Expected exactly zero cut cells.";
}

// =================================================================================================
TEST(CutGrid, PiecewiseLinearSingleLine) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, Float> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {0.7, 0.0},
      {0.0, 0.5},
  };
  const auto success = grid.cut_piecewise_linear<ExtendType::MAX>(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 6UZ) << "Expected exactly six cut cells.";

  EXPECT_EQ(cut_idxs[0], 3UZ);
  EXPECT_EQ(cut_idxs[1], 2UZ);
  EXPECT_EQ(cut_idxs[2], 7UZ);
  EXPECT_EQ(cut_idxs[3], 6UZ);
  EXPECT_EQ(cut_idxs[4], 5UZ);
  EXPECT_EQ(cut_idxs[5], 10UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[0]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[1]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[2]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[3]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[4]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[5]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[5]].get_cut().exit_loc, LEFT);

  for (size_t i = 0; i < cut_idxs.size() - 1; ++i) {
    const auto p1 = grid[cut_idxs[i]].cut_exit<SIM_C>();
    const auto p2 = grid[cut_idxs[i + 1]].cut_entry<SIM_C>();
    EXPECT_DOUBLE_EQ((p1 - p2).norm(), 0.0);
  }
}

// =================================================================================================
TEST(CutGrid, PiecewiseLinearVertLineOnGrid) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, Float> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {0.8, 0.0},
      {0.8, 1.0},
  };
  const auto success = grid.cut_piecewise_linear<ExtendType::MAX>(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 5UZ) << "Expected exactly five cut cells.";

  EXPECT_EQ(cut_idxs[0], 4UZ);
  EXPECT_EQ(cut_idxs[1], 9UZ);
  EXPECT_EQ(cut_idxs[2], 14UZ);
  EXPECT_EQ(cut_idxs[3], 19UZ);
  EXPECT_EQ(cut_idxs[4], 24UZ);

  for (size_t idx : cut_idxs) {
    EXPECT_EQ(grid[idx].get_cut().entry_loc, BOTTOM);
    EXPECT_EQ(grid[idx].get_cut().exit_loc, TOP);
  }

  for (size_t i = 0; i < cut_idxs.size() - 1; ++i) {
    const auto p1 = grid[cut_idxs[i]].cut_exit<SIM_C>();
    const auto p2 = grid[cut_idxs[i + 1]].cut_entry<SIM_C>();
    EXPECT_DOUBLE_EQ((p1 - p2).norm(), 0.0);
  }
}

// =================================================================================================
TEST(CutGrid, PiecewiseLinearHoriLineOnGrid) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, Float> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {0.0, 0.8},
      {1.0, 0.8},
  };
  const auto success = grid.cut_piecewise_linear<ExtendType::MAX>(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 5UZ) << "Expected exactly five cut cells.";

  EXPECT_EQ(cut_idxs[0], 20UZ);
  EXPECT_EQ(cut_idxs[1], 21UZ);
  EXPECT_EQ(cut_idxs[2], 22UZ);
  EXPECT_EQ(cut_idxs[3], 23UZ);
  EXPECT_EQ(cut_idxs[4], 24UZ);

  for (size_t idx : cut_idxs) {
    EXPECT_EQ(grid[idx].get_cut().entry_loc, LEFT);
    EXPECT_EQ(grid[idx].get_cut().exit_loc, RIGHT);
  }

  for (size_t i = 0; i < cut_idxs.size() - 1; ++i) {
    const auto p1 = grid[cut_idxs[i]].cut_exit<SIM_C>();
    const auto p2 = grid[cut_idxs[i + 1]].cut_entry<SIM_C>();
    EXPECT_DOUBLE_EQ((p1 - p2).norm(), 0.0);
  }
}

// =================================================================================================
TEST(CutGrid, PiecewiseLinearMultiplePointsNoExtend) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, Float> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {0.85, 0.35},
      {0.4, 0.5},
      {0.25, 0.9},
  };
  const auto success = grid.cut_piecewise_linear<ExtendType::NONE>(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 5UZ) << "Expected exactly five cut cells.";

  EXPECT_EQ(cut_idxs[0], 8UZ);
  EXPECT_EQ(cut_idxs[1], 13UZ);
  EXPECT_EQ(cut_idxs[2], 12UZ);
  EXPECT_EQ(cut_idxs[3], 11UZ);
  EXPECT_EQ(cut_idxs[4], 16UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[0]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[1]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[2]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[3]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[4]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().exit_loc, TOP);

  for (size_t i = 0; i < cut_idxs.size() - 1; ++i) {
    const auto p1 = grid[cut_idxs[i]].cut_exit<SIM_C>();
    const auto p2 = grid[cut_idxs[i + 1]].cut_entry<SIM_C>();
    EXPECT_DOUBLE_EQ((p1 - p2).norm(), 0.0);
  }
}

// =================================================================================================
TEST(CutGrid, PiecewiseLinearMultiplePointsExtendMax) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, Float> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {0.85, 0.35},
      {0.4, 0.5},
      {0.25, 0.9},
  };
  const auto success = grid.cut_piecewise_linear<ExtendType::MAX>(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 7UZ) << "Expected exactly seven cut cells.";

  EXPECT_EQ(cut_idxs[0], 9UZ);
  EXPECT_EQ(cut_idxs[1], 8UZ);
  EXPECT_EQ(cut_idxs[2], 13UZ);
  EXPECT_EQ(cut_idxs[3], 12UZ);
  EXPECT_EQ(cut_idxs[4], 11UZ);
  EXPECT_EQ(cut_idxs[5], 16UZ);
  EXPECT_EQ(cut_idxs[6], 21UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[0]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[1]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[2]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[3]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[4]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[5]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[5]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[6]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[6]].get_cut().exit_loc, TOP);

  for (size_t i = 0; i < cut_idxs.size() - 1; ++i) {
    const auto p1 = grid[cut_idxs[i]].cut_exit<SIM_C>();
    const auto p2 = grid[cut_idxs[i + 1]].cut_entry<SIM_C>();
    EXPECT_DOUBLE_EQ((p1 - p2).norm(), 0.0);
  }
}

// =================================================================================================
TEST(CutGrid, PiecewiseLinearMultiplePointsExtendNearest) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, Float> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {0.85, 0.35},
      {0.4, 0.5},
      {0.25, 0.91},
  };
  const auto success = grid.cut_piecewise_linear<ExtendType::NEAREST>(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 6UZ) << "Expected exactly six cut cells.";

  EXPECT_EQ(cut_idxs[0], 8UZ);
  EXPECT_EQ(cut_idxs[1], 13UZ);
  EXPECT_EQ(cut_idxs[2], 12UZ);
  EXPECT_EQ(cut_idxs[3], 11UZ);
  EXPECT_EQ(cut_idxs[4], 16UZ);
  EXPECT_EQ(cut_idxs[5], 21UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[0]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[1]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[2]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[3]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[4]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[5]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[5]].get_cut().exit_loc, TOP);

  for (size_t i = 0; i < cut_idxs.size() - 1; ++i) {
    const auto p1 = grid[cut_idxs[i]].cut_exit<SIM_C>();
    const auto p2 = grid[cut_idxs[i + 1]].cut_entry<SIM_C>();
    EXPECT_DOUBLE_EQ((p1 - p2).norm(), 0.0);
  }
}

// =================================================================================================
TEST(CutGrid, CutOnCorner) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, Float> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {0.8, 0.2},
      {0.4, 0.4},
      {0.0, 0.8},
  };
  const auto success = grid.cut_piecewise_linear<ExtendType::NEAREST>(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 4UZ) << "Expected exactly four cut cells.";

  EXPECT_EQ(cut_idxs[0], 8UZ);
  EXPECT_EQ(cut_idxs[1], 7UZ);
  EXPECT_EQ(cut_idxs[2], 11UZ);
  EXPECT_EQ(cut_idxs[3], 15UZ);

  // Multiple acceptable cut types
  EXPECT_GE(grid[cut_idxs[0]].get_cut().entry_loc & (BOTTOM | RIGHT), 0);
  EXPECT_GE(grid[cut_idxs[0]].get_cut().exit_loc & LEFT, 0);

  EXPECT_GE(grid[cut_idxs[1]].get_cut().entry_loc & RIGHT, 0);
  EXPECT_GE(grid[cut_idxs[1]].get_cut().exit_loc & (LEFT | TOP), 0);

  EXPECT_GE(grid[cut_idxs[2]].get_cut().entry_loc & (RIGHT | BOTTOM), 0);
  EXPECT_GE(grid[cut_idxs[2]].get_cut().exit_loc & (TOP | LEFT), 0);

  EXPECT_GE(grid[cut_idxs[3]].get_cut().entry_loc & (RIGHT | BOTTOM), 0);
  EXPECT_GE(grid[cut_idxs[3]].get_cut().exit_loc & (TOP | LEFT), 0);

  for (size_t i = 0; i < cut_idxs.size() - 1; ++i) {
    const auto p1 = grid[cut_idxs[i]].cut_exit<SIM_C>();
    const auto p2 = grid[cut_idxs[i + 1]].cut_entry<SIM_C>();
    EXPECT_DOUBLE_EQ((p1 - p2).norm(), 0.0);
  }
}

// =================================================================================================
TEST(CutGrid, PointsFullyInside) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, Float> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {0.85, 0.35},
      {0.55, 0.5},
      {0.45, 0.5},
      {0.4, 0.5},
      {0.25, 0.91},
  };
  const auto success = grid.cut_piecewise_linear<ExtendType::MAX>(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 7UZ) << "Expected exactly seven cut cells.";

  EXPECT_EQ(cut_idxs[0], 9UZ);
  EXPECT_EQ(cut_idxs[1], 8UZ);
  EXPECT_EQ(cut_idxs[2], 13UZ);
  EXPECT_EQ(cut_idxs[3], 12UZ);
  EXPECT_EQ(cut_idxs[4], 11UZ);
  EXPECT_EQ(cut_idxs[5], 16UZ);
  EXPECT_EQ(cut_idxs[6], 21UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[0]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[1]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[2]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[3]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().exit_loc, LEFT);

  EXPECT_EQ(grid[cut_idxs[4]].get_cut().entry_loc, RIGHT);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[5]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[5]].get_cut().exit_loc, TOP);

  EXPECT_EQ(grid[cut_idxs[6]].get_cut().entry_loc, BOTTOM);
  EXPECT_EQ(grid[cut_idxs[6]].get_cut().exit_loc, TOP);

  for (size_t i = 0; i < cut_idxs.size() - 1; ++i) {
    const auto p1 = grid[cut_idxs[i]].cut_exit<SIM_C>();
    const auto p2 = grid[cut_idxs[i + 1]].cut_entry<SIM_C>();
    EXPECT_DOUBLE_EQ((p1 - p2).norm(), 0.0);
  }
}

// =================================================================================================
TEST(CutGrid, HandleCutOnSingleSide) {
  // Synthetic example
  {
    UniformGrid<double, double> grid(0.0, 1.0, 2, 0.0, 1.0, 2);
    std::vector<SimCoord<double>> points = {
        {0.125, 0.0},
        {0.5 + 50 * EPS<double>(), 0.5 + 1e-2},
        {0.375, 1.0},
    };

    const auto success = grid.cut_piecewise_linear<ExtendType::MAX>(points);
    EXPECT_TRUE(success) << "Cut not cut the grid along the given points.";

    const auto cut_idxs = grid.cut_cell_idxs();
    ASSERT_EQ(cut_idxs.size(), 2UZ);

    EXPECT_EQ(cut_idxs[0], 0UZ);
    EXPECT_EQ(cut_idxs[1], 2UZ);

    EXPECT_EQ(grid[cut_idxs[0]].get_cut().entry_loc, BOTTOM);
    EXPECT_EQ(grid[cut_idxs[0]].get_cut().exit_loc, TOP);
    EXPECT_EQ(grid[cut_idxs[1]].get_cut().entry_loc, BOTTOM);
    EXPECT_EQ(grid[cut_idxs[1]].get_cut().exit_loc, TOP);

    EXPECT_DOUBLE_EQ(
        (grid[cut_idxs[0]].cut_exit<SIM_C>() - grid[cut_idxs[1]].cut_entry<SIM_C>()).norm(), 0.0);
  }

  // Real example
  {
    UniformGrid<double, double> grid(0.0, 1.0, 13, 0.0, 1.0, 13);
    std::vector<GridCoord<double>> points = {
        {8.920254698767621, 1.9924984252186384},
        {9.000887022549623, 3.0133103606505567},
        {8.817575199595618, 4.042412896103696},
    };

    const auto success = grid.cut_piecewise_linear<ExtendType::MAX>(points);
    EXPECT_TRUE(success) << "Cut not cut the grid along the given points.";

    const auto cut_idxs = grid.cut_cell_idxs();

    for (size_t i = 0; i < cut_idxs.size() - 1; ++i) {
      EXPECT_DOUBLE_EQ(
          (grid[cut_idxs[i]].cut_exit<SIM_C>() - grid[cut_idxs[i + 1]].cut_entry<SIM_C>()).norm(),
          0.0);
    }
  }

  // Entry point on corner
  {
    UniformGrid<double, double> grid(0.0, 1.0, 2, 0.0, 1.0, 2);
    std::vector<SimCoord<double>> points = {
        {0.0, 0.0},
        {0.625, 0.625},
        {0.25, 1.0},
    };

    const auto success = grid.cut_piecewise_linear<ExtendType::MAX>(points);
    EXPECT_TRUE(success) << "Cut not cut the grid along the given points.";

    const auto cut_idxs = grid.cut_cell_idxs();
    ASSERT_EQ(cut_idxs.size(), 2UZ);

    EXPECT_EQ(cut_idxs[0], 0UZ);
    EXPECT_EQ(cut_idxs[1], 2UZ);

    EXPECT_EQ(grid[cut_idxs[0]].get_cut().entry_loc, BOTTOM);
    EXPECT_EQ(grid[cut_idxs[0]].get_cut().exit_loc, TOP);
    EXPECT_EQ(grid[cut_idxs[1]].get_cut().entry_loc, BOTTOM);
    EXPECT_EQ(grid[cut_idxs[1]].get_cut().exit_loc, TOP);

    EXPECT_DOUBLE_EQ(
        (grid[cut_idxs[0]].cut_exit<SIM_C>() - grid[cut_idxs[1]].cut_entry<SIM_C>()).norm(), 0.0);
    EXPECT_DOUBLE_EQ((grid[cut_idxs[0]].cut_exit<SIM_C>() - SimCoord<double>{0.5, 0.5}).norm(),
                     0.0);
  }

  // Exit point on corner
  {
    UniformGrid<double, double> grid(0.0, 1.0, 2, 0.0, 1.0, 2);
    std::vector<SimCoord<double>> points = {
        {0.25, 0.0},
        {0.625, 0.375},
        {0.0, 1.0},
    };

    const auto success = grid.cut_piecewise_linear<ExtendType::MAX>(points);
    EXPECT_TRUE(success) << "Cut not cut the grid along the given points.";

    const auto cut_idxs = grid.cut_cell_idxs();
    ASSERT_EQ(cut_idxs.size(), 2UZ);

    EXPECT_EQ(cut_idxs[0], 0UZ);
    EXPECT_EQ(cut_idxs[1], 2UZ);

    EXPECT_EQ(grid[cut_idxs[0]].get_cut().entry_loc, BOTTOM);
    EXPECT_EQ(grid[cut_idxs[0]].get_cut().exit_loc, TOP);
    EXPECT_EQ(grid[cut_idxs[1]].get_cut().entry_loc, BOTTOM);
    EXPECT_EQ(grid[cut_idxs[1]].get_cut().exit_loc, TOP);

    EXPECT_DOUBLE_EQ(
        (grid[cut_idxs[0]].cut_exit<SIM_C>() - grid[cut_idxs[1]].cut_entry<SIM_C>()).norm(), 0.0);
    EXPECT_DOUBLE_EQ((grid[cut_idxs[0]].cut_exit<SIM_C>() - SimCoord<double>{0.5, 0.5}).norm(),
                     0.0);
  }

  // Entry and exit point on corner
  {
    UniformGrid<double, double> grid(0.0, 1.0, 3, 0.0, 1.0, 3);
    std::vector<GridCoord<double>> points = {
        {0.0, 0.0},
        {1.5, 1.5},
        {0.0, 3.0},
    };

    const auto success = grid.cut_piecewise_linear<ExtendType::MAX>(points);
    EXPECT_TRUE(success) << "Cut not cut the grid along the given points.";

    const auto cut_idxs = grid.cut_cell_idxs();
    ASSERT_EQ(cut_idxs.size(), 3UZ);

    EXPECT_EQ(cut_idxs[0], 0UZ);
    EXPECT_EQ(cut_idxs[1], 4UZ);
    EXPECT_EQ(cut_idxs[2], 6UZ);

    EXPECT_EQ(grid[cut_idxs[0]].get_cut().entry_loc, BOTTOM);
    EXPECT_EQ(grid[cut_idxs[0]].get_cut().exit_loc, TOP);
    EXPECT_EQ(grid[cut_idxs[1]].get_cut().entry_loc, BOTTOM);
    EXPECT_EQ(grid[cut_idxs[1]].get_cut().exit_loc, TOP);
    EXPECT_EQ(grid[cut_idxs[2]].get_cut().entry_loc, BOTTOM);
    EXPECT_EQ(grid[cut_idxs[2]].get_cut().exit_loc, TOP);

    for (size_t i = 0; i < cut_idxs.size() - 1; ++i) {
      EXPECT_DOUBLE_EQ(
          (grid[cut_idxs[i]].cut_exit<SIM_C>() - grid[cut_idxs[i + 1]].cut_entry<SIM_C>()).norm(),
          0.0);
    }
  }
}

// =================================================================================================
TEST(CutGrid, PiecewiseLinearRoundingBug) {
  UniformGrid<double, double> grid(0.0, 5.0, 5UZ, 0.0, 5.0, 5UZ);

  const std::vector<GridCoord<double>> points = {
      {3.2663460232270394, 0.05915286922710798},
      {3.3001116149507848, 1.0830318867396005},
      {3.1653095473707733, 2.0765210702248598},
      {3.0838033071117104, 2.267462824652602},
      {2.2576268619934536, 3.0912892589539536},
      {2.0824742005628303, 3.1510810962845373},
      {1.081172302885356, 3.279781651119647},
      {0.05849934452716325, 3.254446826278905},
  };
  const bool success = grid.cut_piecewise_linear<ExtendType::NEAREST>(points);
  ASSERT_TRUE(success) << "Could not cut the grid.";

  EXPECT_EQ(grid.cut_cell_idxs().size(), 7UZ)
      << "Expected 7 cut cells but got only " << grid.cut_cell_idxs().size();

  EXPECT_EQ(grid.cut_cell_idxs(), (std::vector<size_t>{3UZ, 8UZ, 13UZ, 12UZ, 17UZ, 16UZ, 15UZ}));
}
