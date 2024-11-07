#include <gtest/gtest.h>

#include <CellBased/Cell.hpp>
#include <CellBased/Interface.hpp>

using namespace Zap::CellBased;

using Float = double;

TEST(OuterInterfaces, A) {
  Cell<Float, Float> cell = {
      .cell_type =
          CutValue<Float>{
              .left_value    = 42.0,
              .right_value   = 39.0,
              .rel_cut_entry = {0.0, 0.5},
              .entry_loc     = LEFT,
              .rel_cut_exit  = {0.5, 1.0},
              .exit_loc      = TOP,
          },
      .m_x_min = 0.0,
      .m_y_min = 0.0,
      .m_dx    = 1.0,
      .m_dy    = 1.0,
  };

  const auto outer_interfaces = get_outer_cell_interfaces<Float, Float, GridCoord<Float>>(cell);

  const auto left_side = outer_interfaces.left;
  ASSERT_EQ(left_side.size(), 2UZ);
  EXPECT_LE((left_side[0].begin - GridCoord<Float>{0, 0}).norm(), EPS<Float>);
  EXPECT_LE((left_side[0].end - GridCoord<Float>{0, 0.5}).norm(), EPS<Float>);
  EXPECT_DOUBLE_EQ(left_side[0].value, 39.0);
  EXPECT_LE((left_side[1].begin - GridCoord<Float>{0, 0.5}).norm(), EPS<Float>);
  EXPECT_LE((left_side[1].end - GridCoord<Float>{0.0, 1.0}).norm(), EPS<Float>);
  EXPECT_DOUBLE_EQ(left_side[1].value, 42.0);

  const auto right_side = outer_interfaces.right;
  ASSERT_EQ(right_side.size(), 1UZ);
  IGOR_DEBUG_PRINT(right_side[0].begin);
  IGOR_DEBUG_PRINT(right_side[0].end);
  EXPECT_LE((right_side[0].begin - GridCoord<Float>{1.0, 0.0}).norm(), EPS<Float>);
  EXPECT_LE((right_side[0].end - GridCoord<Float>{1.0, 1.0}).norm(), EPS<Float>);
  EXPECT_DOUBLE_EQ(right_side[0].value, 39.0);

  const auto bottom_side = outer_interfaces.bottom;
  ASSERT_EQ(bottom_side.size(), 1UZ);
  EXPECT_LE((bottom_side[0].begin - GridCoord<Float>{1.0, 0.0}).norm(), EPS<Float>);
  EXPECT_LE((bottom_side[0].end - GridCoord<Float>{0.0, 0.0}).norm(), EPS<Float>);
  EXPECT_DOUBLE_EQ(bottom_side[0].value, 39.0);

  const auto top_side = outer_interfaces.top;
  ASSERT_EQ(top_side.size(), 2UZ);
  EXPECT_LE((top_side[0].begin - GridCoord<Float>{0.5, 1.0}).norm(), EPS<Float>);
  EXPECT_LE((top_side[0].end - GridCoord<Float>{0.0, 1.0}).norm(), EPS<Float>);
  EXPECT_DOUBLE_EQ(top_side[0].value, 42.0);
  EXPECT_LE((top_side[1].begin - GridCoord<Float>{1.0, 1.0}).norm(), EPS<Float>);
  EXPECT_LE((top_side[1].end - GridCoord<Float>{0.5, 1.0}).norm(), EPS<Float>);
  EXPECT_DOUBLE_EQ(top_side[1].value, 39.0);
}
