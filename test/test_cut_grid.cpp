#include <gtest/gtest.h>

#include <numbers>

#include <CellBased/Grid.hpp>

using namespace Zap::CellBased;

TEST(CutGrid, CurveQuarterCircle) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, 1> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

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
  ASSERT_EQ(cut_idxs.size(), 5UZ) << "Expect exactly five cut cells.";

  EXPECT_EQ(cut_idxs[0], 2UZ);
  EXPECT_EQ(cut_idxs[1], 7UZ);
  EXPECT_EQ(cut_idxs[2], 6UZ);
  EXPECT_EQ(cut_idxs[3], 11UZ);
  EXPECT_EQ(cut_idxs[4], 10UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().type, CutType::MIDDLE_VERT);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().type, CutType::BOTTOM_LEFT);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().type, CutType::TOP_RIGHT);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().type, CutType::BOTTOM_LEFT);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().type, CutType::MIDDLE_HORI);
}

TEST(CutGrid, PiecewiseLinearSingleLine) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, 1> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {.x = 0.7, .y = 0.0},
      {.x = 0.0, .y = 0.5},
  };
  const auto success = grid.cut_piecewise_linear(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 6UZ) << "Expect exactly five cut cells.";

  EXPECT_EQ(cut_idxs[0], 3UZ);
  EXPECT_EQ(cut_idxs[1], 2UZ);
  EXPECT_EQ(cut_idxs[2], 7UZ);
  EXPECT_EQ(cut_idxs[3], 6UZ);
  EXPECT_EQ(cut_idxs[4], 5UZ);
  EXPECT_EQ(cut_idxs[5], 10UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().type, CutType::BOTTOM_LEFT);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().type, CutType::TOP_RIGHT);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().type, CutType::BOTTOM_LEFT);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().type, CutType::MIDDLE_HORI);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().type, CutType::TOP_RIGHT);
  EXPECT_EQ(grid[cut_idxs[5]].get_cut().type, CutType::BOTTOM_LEFT);
}

TEST(CutGrid, PiecewiseLinearVertLineOnGrid) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, 1> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {.x = 0.8, .y = 0.0},
      {.x = 0.8, .y = 1.0},
  };
  const auto success = grid.cut_piecewise_linear(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 5UZ) << "Expect exactly five cut cells.";

  EXPECT_EQ(cut_idxs[0], 4UZ);
  EXPECT_EQ(cut_idxs[1], 9UZ);
  EXPECT_EQ(cut_idxs[2], 14UZ);
  EXPECT_EQ(cut_idxs[3], 19UZ);
  EXPECT_EQ(cut_idxs[4], 24UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().type, CutType::MIDDLE_VERT);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().type, CutType::MIDDLE_VERT);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().type, CutType::MIDDLE_VERT);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().type, CutType::MIDDLE_VERT);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().type, CutType::MIDDLE_VERT);
}

TEST(CutGrid, PiecewiseLinearHoriLineOnGrid) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, 1> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {.x = 0.0, .y = 0.8},
      {.x = 1.0, .y = 0.8},
  };
  const auto success = grid.cut_piecewise_linear(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 5UZ) << "Expect exactly five cut cells.";

  EXPECT_EQ(cut_idxs[0], 20UZ);
  EXPECT_EQ(cut_idxs[1], 21UZ);
  EXPECT_EQ(cut_idxs[2], 22UZ);
  EXPECT_EQ(cut_idxs[3], 23UZ);
  EXPECT_EQ(cut_idxs[4], 24UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().type, CutType::MIDDLE_HORI);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().type, CutType::MIDDLE_HORI);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().type, CutType::MIDDLE_HORI);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().type, CutType::MIDDLE_HORI);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().type, CutType::MIDDLE_HORI);
}

TEST(CutGrid, PiecewiseLinearMultiplePointsNoExtend) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, 1> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {.x = 0.85, .y = 0.35},
      {.x = 0.4, .y = 0.5},
      {.x = 0.25, .y = 0.9},
  };
  const auto success = grid.cut_piecewise_linear<ExtendType::NONE>(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 5UZ) << "Expect exactly five cut cells.";

  EXPECT_EQ(cut_idxs[0], 8UZ);
  EXPECT_EQ(cut_idxs[1], 13UZ);
  EXPECT_EQ(cut_idxs[2], 12UZ);
  EXPECT_EQ(cut_idxs[3], 11UZ);
  EXPECT_EQ(cut_idxs[4], 16UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().type, CutType::TOP_RIGHT);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().type, CutType::BOTTOM_LEFT);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().type, CutType::MIDDLE_HORI);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().type, CutType::TOP_RIGHT);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().type, CutType::MIDDLE_VERT);
}

TEST(CutGrid, PiecewiseLinearMultiplePointsExtendMax) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, 1> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {.x = 0.85, .y = 0.35},
      {.x = 0.4, .y = 0.5},
      {.x = 0.25, .y = 0.9},
  };
  const auto success = grid.cut_piecewise_linear<ExtendType::MAX>(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 7UZ) << "Expect exactly five cut cells.";

  EXPECT_EQ(cut_idxs[0], 9UZ);
  EXPECT_EQ(cut_idxs[1], 8UZ);
  EXPECT_EQ(cut_idxs[2], 13UZ);
  EXPECT_EQ(cut_idxs[3], 12UZ);
  EXPECT_EQ(cut_idxs[4], 11UZ);
  EXPECT_EQ(cut_idxs[5], 16UZ);
  EXPECT_EQ(cut_idxs[6], 21UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().type, CutType::MIDDLE_HORI);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().type, CutType::TOP_RIGHT);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().type, CutType::BOTTOM_LEFT);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().type, CutType::MIDDLE_HORI);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().type, CutType::TOP_RIGHT);
  EXPECT_EQ(grid[cut_idxs[5]].get_cut().type, CutType::MIDDLE_VERT);
  EXPECT_EQ(grid[cut_idxs[6]].get_cut().type, CutType::MIDDLE_VERT);
}

TEST(CutGrid, PiecewiseLinearMultiplePointsExtendNearest) {
  using Float = double;

  const Float x_min = 0.0;
  const Float x_max = 1.0;
  const Float y_min = 0.0;
  const Float y_max = 1.0;

  UniformGrid<Float, 1> grid(x_min, x_max, 5UZ, y_min, y_max, 5UZ);

  std::vector<SimCoord<Float>> points = {
      {.x = 0.85, .y = 0.35},
      {.x = 0.4, .y = 0.5},
      {.x = 0.25, .y = 0.91},
  };
  const auto success = grid.cut_piecewise_linear<ExtendType::NEAREST>(points);
  ASSERT_TRUE(success) << "Could not cut the grid successfully.";

  const auto cut_idxs = grid.cut_cell_idxs();
  ASSERT_EQ(cut_idxs.size(), 6UZ) << "Expect exactly five cut cells.";

  EXPECT_EQ(cut_idxs[0], 8UZ);
  EXPECT_EQ(cut_idxs[1], 13UZ);
  EXPECT_EQ(cut_idxs[2], 12UZ);
  EXPECT_EQ(cut_idxs[3], 11UZ);
  EXPECT_EQ(cut_idxs[4], 16UZ);
  EXPECT_EQ(cut_idxs[5], 21UZ);

  EXPECT_EQ(grid[cut_idxs[0]].get_cut().type, CutType::TOP_RIGHT);
  EXPECT_EQ(grid[cut_idxs[1]].get_cut().type, CutType::BOTTOM_LEFT);
  EXPECT_EQ(grid[cut_idxs[2]].get_cut().type, CutType::MIDDLE_HORI);
  EXPECT_EQ(grid[cut_idxs[3]].get_cut().type, CutType::TOP_RIGHT);
  EXPECT_EQ(grid[cut_idxs[4]].get_cut().type, CutType::MIDDLE_VERT);
  EXPECT_EQ(grid[cut_idxs[5]].get_cut().type, CutType::MIDDLE_VERT);
}
