#include <gtest/gtest.h>

#include "CellBased/Cell.hpp"
#include "CellBased/Geometry.hpp"

using namespace Zap::CellBased;

using Float          = double;
constexpr size_t DIM = 2;

TEST(CellGeometry, AreaBottomLeft) {
  const Cell<Float, Float, DIM> cell{
      .value =
          CutValue<Float, DIM>{
              .left_value  = Eigen::Vector<Float, DIM>::Zero(),
              .right_value = Eigen::Vector<Float, DIM>::Zero(),
              .type        = CutType::BOTTOM_LEFT,
              .rel_cut1    = {0.5, 0.0},
              .rel_cut2    = {0.0, 0.5},
          },
      .m_x_min = -5,
      .m_y_min = 1,
      .m_dx    = 0.5,
      .m_dy    = 0.25,
  };

  const auto cartesian_polygon = cell.template get_cartesian_polygon<SIM_C>();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), cell.template dx<SIM_C>() * cell.template dy<SIM_C>());

  const auto left_cut_polygon  = cell.template get_cut_left_polygon<SIM_C>();
  const auto right_cut_polygon = cell.template get_cut_right_polygon<SIM_C>();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), left_cut_polygon.area() + right_cut_polygon.area());

  EXPECT_DOUBLE_EQ(left_cut_polygon.area(),
                   (cell.template cut1<SIM_C>().x - cell.template x_min<SIM_C>()) *
                       (cell.template cut2<SIM_C>().y - cell.template y_min<SIM_C>()) / 2);
  EXPECT_DOUBLE_EQ(right_cut_polygon.area(),
                   cell.template dx<SIM_C>() * cell.template dy<SIM_C>() -
                       (cell.template cut1<SIM_C>().x - cell.template x_min<SIM_C>()) *
                           (cell.template cut2<SIM_C>().y - cell.template y_min<SIM_C>()) / 2);

  EXPECT_DOUBLE_EQ(cartesian_polygon.area(),
                   cell.template get_cartesian_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
  EXPECT_DOUBLE_EQ(left_cut_polygon.area(),
                   cell.template get_cut_left_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
  EXPECT_DOUBLE_EQ(right_cut_polygon.area(),
                   cell.template get_cut_right_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
}

TEST(CellGeometry, AreaBottomRight) {
  const Cell<Float, Float, DIM> cell{
      .value =
          CutValue<Float, DIM>{
              .left_value  = Eigen::Vector<Float, DIM>::Zero(),
              .right_value = Eigen::Vector<Float, DIM>::Zero(),
              .type        = CutType::BOTTOM_RIGHT,
              .rel_cut1    = {0.5, 0.0},
              .rel_cut2    = {1.0, 0.5},
          },
      .m_x_min = -5,
      .m_y_min = 1,
      .m_dx    = 0.5,
      .m_dy    = 0.25,
  };

  const auto cartesian_polygon = cell.template get_cartesian_polygon<SIM_C>();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), cell.template dx<SIM_C>() * cell.template dy<SIM_C>());

  const auto left_cut_polygon  = cell.template get_cut_left_polygon<SIM_C>();
  const auto right_cut_polygon = cell.template get_cut_right_polygon<SIM_C>();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), left_cut_polygon.area() + right_cut_polygon.area());

  EXPECT_DOUBLE_EQ(left_cut_polygon.area(),
                   cell.template dx<SIM_C>() * cell.template dy<SIM_C>() -
                       (cell.template cut1<SIM_C>().x - cell.template x_min<SIM_C>()) *
                           (cell.template cut2<SIM_C>().y - cell.template y_min<SIM_C>()) / 2);
  EXPECT_DOUBLE_EQ(
      right_cut_polygon.area(),
      (cell.template x_min<SIM_C>() + cell.template dx<SIM_C>() - cell.template cut1<SIM_C>().x) *
          (cell.template cut2<SIM_C>().y - cell.template y_min<SIM_C>()) / 2);

  EXPECT_DOUBLE_EQ(cartesian_polygon.area(),
                   cell.template get_cartesian_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
  EXPECT_DOUBLE_EQ(left_cut_polygon.area(),
                   cell.template get_cut_left_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
  EXPECT_DOUBLE_EQ(right_cut_polygon.area(),
                   cell.template get_cut_right_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
}

TEST(CellGeometry, AreaTopRight) {
  const Cell<Float, Float, DIM> cell{
      .value =
          CutValue<Float, DIM>{
              .left_value  = Eigen::Vector<Float, DIM>::Zero(),
              .right_value = Eigen::Vector<Float, DIM>::Zero(),
              .type        = CutType::TOP_RIGHT,
              .rel_cut1    = {1.0, 0.5},
              .rel_cut2    = {0.5, 1.0},
          },
      .m_x_min = -5,
      .m_y_min = 1,
      .m_dx    = 0.5,
      .m_dy    = 0.25,
  };

  const auto cartesian_polygon = cell.template get_cartesian_polygon<SIM_C>();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), cell.template dx<SIM_C>() * cell.template dy<SIM_C>());

  const auto left_cut_polygon  = cell.template get_cut_left_polygon<SIM_C>();
  const auto right_cut_polygon = cell.template get_cut_right_polygon<SIM_C>();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), left_cut_polygon.area() + right_cut_polygon.area());

  EXPECT_DOUBLE_EQ(left_cut_polygon.area(),
                   cell.template dx<SIM_C>() * cell.template dy<SIM_C>() -
                       (cell.template x_min<SIM_C>() + cell.template dx<SIM_C>() -
                        cell.template cut2<SIM_C>().x) *
                           (cell.template y_min<SIM_C>() + cell.template dy<SIM_C>() -
                            cell.template cut1<SIM_C>().y) /
                           2);
  EXPECT_DOUBLE_EQ(
      right_cut_polygon.area(),
      (cell.template x_min<SIM_C>() + cell.template dx<SIM_C>() - cell.template cut2<SIM_C>().x) *
          (cell.template y_min<SIM_C>() + cell.template dy<SIM_C>() -
           cell.template cut1<SIM_C>().y) /
          2);

  EXPECT_DOUBLE_EQ(cartesian_polygon.area(),
                   cell.template get_cartesian_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
  EXPECT_DOUBLE_EQ(left_cut_polygon.area(),
                   cell.template get_cut_left_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
  EXPECT_DOUBLE_EQ(right_cut_polygon.area(),
                   cell.template get_cut_right_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
}

TEST(CellGeometry, AreaTopLeft) {
  const Cell<Float, Float, DIM> cell{
      .value =
          CutValue<Float, DIM>{
              .left_value  = Eigen::Vector<Float, DIM>::Zero(),
              .right_value = Eigen::Vector<Float, DIM>::Zero(),
              .type        = CutType::TOP_LEFT,
              .rel_cut1    = {0.5, 1.0},
              .rel_cut2    = {0.0, 0.5},
          },
      .m_x_min = -5,
      .m_y_min = 1,
      .m_dx    = 0.5,
      .m_dy    = 0.25,
  };

  const auto cartesian_polygon = cell.template get_cartesian_polygon<SIM_C>();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), cell.template dx<SIM_C>() * cell.template dy<SIM_C>());

  const auto left_cut_polygon  = cell.template get_cut_left_polygon<SIM_C>();
  const auto right_cut_polygon = cell.template get_cut_right_polygon<SIM_C>();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), left_cut_polygon.area() + right_cut_polygon.area());

  EXPECT_DOUBLE_EQ(left_cut_polygon.area(),
                   cell.template dx<SIM_C>() * cell.template dy<SIM_C>() -
                       (cell.template x_min<SIM_C>() + cell.template dx<SIM_C>() -
                        cell.template cut1<SIM_C>().x) *
                           (cell.template y_min<SIM_C>() + cell.template dy<SIM_C>() -
                            cell.template cut2<SIM_C>().y) /
                           2);
  EXPECT_DOUBLE_EQ(
      right_cut_polygon.area(),
      (cell.template x_min<SIM_C>() + cell.template dx<SIM_C>() - cell.template cut1<SIM_C>().x) *
          (cell.template y_min<SIM_C>() + cell.template dy<SIM_C>() -
           cell.template cut2<SIM_C>().y) /
          2);

  EXPECT_DOUBLE_EQ(cartesian_polygon.area(),
                   cell.template get_cartesian_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
  EXPECT_DOUBLE_EQ(left_cut_polygon.area(),
                   cell.template get_cut_left_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
  EXPECT_DOUBLE_EQ(right_cut_polygon.area(),
                   cell.template get_cut_right_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
}

TEST(CellGeometry, AreaMiddleHori) {
  const Cell<Float, Float, DIM> cell{
      .value =
          CutValue<Float, DIM>{
              .left_value  = Eigen::Vector<Float, DIM>::Zero(),
              .right_value = Eigen::Vector<Float, DIM>::Zero(),
              .type        = CutType::MIDDLE_HORI,
              .rel_cut1    = {1.0, 0.5},
              .rel_cut2    = {0.0, 0.04},
          },
      .m_x_min = -5,
      .m_y_min = 1,
      .m_dx    = 0.5,
      .m_dy    = 0.25,
  };

  const auto cartesian_polygon = cell.template get_cartesian_polygon<SIM_C>();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), cell.template dx<SIM_C>() * cell.template dy<SIM_C>());

  const auto left_cut_polygon  = cell.template get_cut_left_polygon<SIM_C>();
  const auto right_cut_polygon = cell.template get_cut_right_polygon<SIM_C>();
  EXPECT_NEAR(cartesian_polygon.area(), left_cut_polygon.area() + right_cut_polygon.area(), 1e-15)
      << "Expect area of two subcells to be equal to the area of the entire cell, difference is "
      << cartesian_polygon.area() - (left_cut_polygon.area() + right_cut_polygon.area()) << '.';

  EXPECT_DOUBLE_EQ(cartesian_polygon.area(),
                   cell.template get_cartesian_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
  EXPECT_NEAR(left_cut_polygon.area(),
              cell.template get_cut_left_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy,
              1e-15);
  EXPECT_NEAR(right_cut_polygon.area(),
              cell.template get_cut_right_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy,
              1e-15);
}

TEST(CellGeometry, AreaMiddleVert) {
  const Cell<Float, Float, DIM> cell{
      .value =
          CutValue<Float, DIM>{
              .left_value  = Eigen::Vector<Float, DIM>::Zero(),
              .right_value = Eigen::Vector<Float, DIM>::Zero(),
              .type        = CutType::MIDDLE_VERT,
              .rel_cut1    = {0.5, 0.0},
              .rel_cut2    = {0.2, 1.0},
          },
      .m_x_min = -5,
      .m_y_min = 1,
      .m_dx    = 0.5,
      .m_dy    = 0.25,
  };

  const auto cartesian_polygon = cell.template get_cartesian_polygon<SIM_C>();
  EXPECT_DOUBLE_EQ(cartesian_polygon.area(), cell.template dx<SIM_C>() * cell.template dy<SIM_C>());

  const auto left_cut_polygon  = cell.template get_cut_left_polygon<SIM_C>();
  const auto right_cut_polygon = cell.template get_cut_right_polygon<SIM_C>();
  EXPECT_NEAR(cartesian_polygon.area(), left_cut_polygon.area() + right_cut_polygon.area(), 1e-15)
      << "Expect area of two subcells to be equal to the area of the entire cell, difference is "
      << cartesian_polygon.area() - (left_cut_polygon.area() + right_cut_polygon.area()) << '.';

  EXPECT_DOUBLE_EQ(cartesian_polygon.area(),
                   cell.template get_cartesian_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy);
  EXPECT_NEAR(left_cut_polygon.area(),
              cell.template get_cut_left_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy,
              1e-15);
  EXPECT_NEAR(right_cut_polygon.area(),
              cell.template get_cut_right_polygon<GRID_C>().area() * cell.m_dx * cell.m_dy,
              1e-15);
}
