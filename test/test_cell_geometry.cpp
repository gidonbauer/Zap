#include <gtest/gtest.h>

#include "CellBased/Cell.hpp"
#include "CellBased/Geometry.hpp"

using namespace Zap::CellBased;

using Float          = double;
constexpr size_t DIM = 2;

TEST(CellGeometry, AreaBottomLeft) {
  const Cell<Float, DIM> cell{
      .value =
          CutValue<Float, DIM>{
              .left_value  = Eigen::Vector<Float, DIM>::Zero(),
              .right_value = Eigen::Vector<Float, DIM>::Zero(),
              .type        = CutType::BOTTOM_LEFT,
              .x1_cut      = -4.75,
              .y1_cut      = 1,
              .x2_cut      = -5,
              .y2_cut      = 1.125,
          },
      .x_min = -5,
      .dx    = 0.5,
      .y_min = 1,
      .dy    = 0.25,
  };

  const auto cartesian_polygon = cell.get_cartesian_polygon();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), cell.dx * cell.dy);

  const auto left_cut_polygon  = cell.get_cut_left_polygon();
  const auto right_cut_polygon = cell.get_cut_right_polygon();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), left_cut_polygon.area() + right_cut_polygon.area());

  EXPECT_DOUBLE_EQ(left_cut_polygon.area(),
                   (cell.get_cut().x1_cut - cell.x_min) * (cell.get_cut().y2_cut - cell.y_min) / 2);
  EXPECT_DOUBLE_EQ(right_cut_polygon.area(),
                   cell.dx * cell.dy - (cell.get_cut().x1_cut - cell.x_min) *
                                           (cell.get_cut().y2_cut - cell.y_min) / 2);
}

TEST(CellGeometry, AreaBottomRight) {
  const Cell<Float, DIM> cell{
      .value =
          CutValue<Float, DIM>{
              .left_value  = Eigen::Vector<Float, DIM>::Zero(),
              .right_value = Eigen::Vector<Float, DIM>::Zero(),
              .type        = CutType::BOTTOM_RIGHT,
              .x1_cut      = -4.75,
              .y1_cut      = 1,
              .x2_cut      = -4.5,
              .y2_cut      = 1.125,
          },
      .x_min = -5,
      .dx    = 0.5,
      .y_min = 1,
      .dy    = 0.25,
  };

  const auto cartesian_polygon = cell.get_cartesian_polygon();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), cell.dx * cell.dy);

  const auto left_cut_polygon  = cell.get_cut_left_polygon();
  const auto right_cut_polygon = cell.get_cut_right_polygon();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), left_cut_polygon.area() + right_cut_polygon.area());

  EXPECT_DOUBLE_EQ(left_cut_polygon.area(),
                   cell.dx * cell.dy - (cell.get_cut().x1_cut - cell.x_min) *
                                           (cell.get_cut().y2_cut - cell.y_min) / 2);
  EXPECT_DOUBLE_EQ(right_cut_polygon.area(),
                   (cell.x_min + cell.dx - cell.get_cut().x1_cut) *
                       (cell.get_cut().y2_cut - cell.y_min) / 2);
}

TEST(CellGeometry, AreaTopRight) {
  const Cell<Float, DIM> cell{
      .value =
          CutValue<Float, DIM>{
              .left_value  = Eigen::Vector<Float, DIM>::Zero(),
              .right_value = Eigen::Vector<Float, DIM>::Zero(),
              .type        = CutType::TOP_RIGHT,
              .x1_cut      = -4.5,
              .y1_cut      = 1.125,
              .x2_cut      = -4.75,
              .y2_cut      = 1.25,
          },
      .x_min = -5,
      .dx    = 0.5,
      .y_min = 1,
      .dy    = 0.25,
  };

  const auto cartesian_polygon = cell.get_cartesian_polygon();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), cell.dx * cell.dy);

  const auto left_cut_polygon  = cell.get_cut_left_polygon();
  const auto right_cut_polygon = cell.get_cut_right_polygon();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), left_cut_polygon.area() + right_cut_polygon.area());

  EXPECT_DOUBLE_EQ(left_cut_polygon.area(),
                   cell.dx * cell.dy - (cell.x_min + cell.dx - cell.get_cut().x2_cut) *
                                           (cell.y_min + cell.dy - cell.get_cut().y1_cut) / 2);
  EXPECT_DOUBLE_EQ(right_cut_polygon.area(),
                   (cell.x_min + cell.dx - cell.get_cut().x2_cut) *
                       (cell.y_min + cell.dy - cell.get_cut().y1_cut) / 2);
}

TEST(CellGeometry, AreaTopLeft) {
  const Cell<Float, DIM> cell{
      .value =
          CutValue<Float, DIM>{
              .left_value  = Eigen::Vector<Float, DIM>::Zero(),
              .right_value = Eigen::Vector<Float, DIM>::Zero(),
              .type        = CutType::TOP_LEFT,
              .x1_cut      = -4.75,
              .y1_cut      = 1.25,
              .x2_cut      = -5,
              .y2_cut      = 1.125,
          },
      .x_min = -5,
      .dx    = 0.5,
      .y_min = 1,
      .dy    = 0.25,
  };

  const auto cartesian_polygon = cell.get_cartesian_polygon();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), cell.dx * cell.dy);

  const auto left_cut_polygon  = cell.get_cut_left_polygon();
  const auto right_cut_polygon = cell.get_cut_right_polygon();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), left_cut_polygon.area() + right_cut_polygon.area());

  EXPECT_DOUBLE_EQ(left_cut_polygon.area(),
                   cell.dx * cell.dy - (cell.x_min + cell.dx - cell.get_cut().x1_cut) *
                                           (cell.y_min + cell.dy - cell.get_cut().y2_cut) / 2);
  EXPECT_DOUBLE_EQ(right_cut_polygon.area(),
                   (cell.x_min + cell.dx - cell.get_cut().x1_cut) *
                       (cell.y_min + cell.dy - cell.get_cut().y2_cut) / 2);
}

TEST(CellGeometry, AreaMiddleHori) {
  const Cell<Float, DIM> cell{
      .value =
          CutValue<Float, DIM>{
              .left_value  = Eigen::Vector<Float, DIM>::Zero(),
              .right_value = Eigen::Vector<Float, DIM>::Zero(),
              .type        = CutType::MIDDLE_HORI,
              .x1_cut      = -4.5,
              .y1_cut      = 1.125,
              .x2_cut      = -5,
              .y2_cut      = 1.01,
          },
      .x_min = -5,
      .dx    = 0.5,
      .y_min = 1,
      .dy    = 0.25,
  };

  const auto cartesian_polygon = cell.get_cartesian_polygon();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), cell.dx * cell.dy);

  const auto left_cut_polygon  = cell.get_cut_left_polygon();
  const auto right_cut_polygon = cell.get_cut_right_polygon();
  EXPECT_NEAR(cartesian_polygon.area(), left_cut_polygon.area() + right_cut_polygon.area(), 1e-15)
      << "Expect area of two subcells to be equal to the area of the entire cell, difference is "
      << cartesian_polygon.area() - (left_cut_polygon.area() + right_cut_polygon.area()) << '.';
}

TEST(CellGeometry, AreaMiddleVert) {
  const Cell<Float, DIM> cell{
      .value =
          CutValue<Float, DIM>{
              .left_value  = Eigen::Vector<Float, DIM>::Zero(),
              .right_value = Eigen::Vector<Float, DIM>::Zero(),
              .type        = CutType::MIDDLE_VERT,
              .x1_cut      = -4.75,
              .y1_cut      = 1,
              .x2_cut      = -4.9,
              .y2_cut      = 1.25,
          },
      .x_min = -5,
      .dx    = 0.5,
      .y_min = 1,
      .dy    = 0.25,
  };

  const auto cartesian_polygon = cell.get_cartesian_polygon();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), cell.dx * cell.dy);

  const auto left_cut_polygon  = cell.get_cut_left_polygon();
  const auto right_cut_polygon = cell.get_cut_right_polygon();
  EXPECT_NEAR(cartesian_polygon.area(), left_cut_polygon.area() + right_cut_polygon.area(), 1e-15)
      << "Expect area of two subcells to be equal to the area of the entire cell, difference is "
      << cartesian_polygon.area() - (left_cut_polygon.area() + right_cut_polygon.area()) << '.';
}
