#include <gtest/gtest.h>

#include "CellBased/Geometry.hpp"
namespace Geo = Zap::CellBased::Geometry;
using Zap::CellBased::SmallVector;

TEST(GeometryPolygon, AddPoint) {
  using Float = double;

  Geo::Polygon<Float> polygon{};
  EXPECT_EQ(polygon.size(), 0);
  EXPECT_DOUBLE_EQ(polygon.area(), 0.0);

  polygon.add_point(Eigen::Vector<Float, 2>{0.0, 0.0});
  EXPECT_EQ(polygon.size(), 1);
  EXPECT_DOUBLE_EQ(polygon.area(), 0.0);

  polygon.add_point(Eigen::Vector<Float, 2>{0.0, 1.0});
  EXPECT_EQ(polygon.size(), 2);
  EXPECT_DOUBLE_EQ(polygon.area(), 0.0);

  polygon.add_point(Eigen::Vector<Float, 2>{1.0, 0.0});
  EXPECT_EQ(polygon.size(), 3);
  EXPECT_DOUBLE_EQ(polygon.area(), 0.5);

  polygon.add_point(Eigen::Vector<Float, 2>{1.0, 1.0});
  EXPECT_EQ(polygon.size(), 4);
  EXPECT_DOUBLE_EQ(polygon.area(), 1.0);

  polygon.add_point(Eigen::Vector<Float, 2>{1.0, 1.0});
  EXPECT_EQ(polygon.size(), 4);
  EXPECT_DOUBLE_EQ(polygon.area(), 1.0);
}

TEST(GeometryPolygon, Area) {
  {
    using Float = float;
    const SmallVector<Eigen::Vector<Float, 2>> points{
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {1.0, 1.0},
    };

    const Geo::Polygon polygon{points};

    EXPECT_EQ(polygon.size(), 4);
    EXPECT_DOUBLE_EQ(polygon.area(), 1.0);
  }

  {
    using Float = double;
    const SmallVector<Eigen::Vector<Float, 2>> points{
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {1.0, 1.0},
    };

    const Geo::Polygon polygon{points};

    EXPECT_EQ(polygon.size(), 4);
    EXPECT_DOUBLE_EQ(polygon.area(), 1.0);
  }

  {
    using Float = double;
    const SmallVector<Eigen::Vector<Float, 2>> points{
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 1.0},
    };

    const Geo::Polygon polygon{points};

    EXPECT_DOUBLE_EQ(polygon.area(), 0.5);
  }

  {
    using Float = double;
    const SmallVector<Eigen::Vector<Float, 2>> points{
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {1.5, 0.5},
    };

    const Geo::Polygon polygon{points};

    EXPECT_DOUBLE_EQ(polygon.area(), 1.25);
  }
}

TEST(GeometryPolygon, Intersect) {
  {
    const Geo::Polygon<double> poly1{{
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
    }};

    const Geo::Polygon<double> poly2{{}};

    const auto intersect1 = Geo::intersection(poly1, poly2);
    const auto intersect2 = Geo::intersection(poly2, poly1);

    EXPECT_EQ(intersect1.size(), intersect2.size());
    EXPECT_DOUBLE_EQ(intersect1.area(), intersect2.area());

    EXPECT_EQ(intersect1.size(), 0);
    EXPECT_DOUBLE_EQ(intersect1.area(), 0.0);
  }

  {
    const Geo::Polygon<double> poly1{{
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
    }};

    const Geo::Polygon<double> poly2{{
        {2.0, 0.0},
        {2.0, 1.0},
        {3.0, 0.0},
    }};

    const auto intersect1 = Geo::intersection(poly1, poly2);
    const auto intersect2 = Geo::intersection(poly2, poly1);

    EXPECT_EQ(intersect1.size(), intersect2.size());
    EXPECT_DOUBLE_EQ(intersect1.area(), intersect2.area());

    EXPECT_EQ(intersect1.size(), 0);
    EXPECT_DOUBLE_EQ(intersect1.area(), 0.0);
  }

  {
    const Geo::Polygon<double> poly1{{
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
    }};

    const Geo::Polygon<double> poly2{{
        {1.0, 0.0},
        {1.0, 1.0},
        {3.0, 0.0},
    }};

    const auto intersect1 = Geo::intersection(poly1, poly2);
    const auto intersect2 = Geo::intersection(poly2, poly1);

    EXPECT_EQ(intersect1.size(), intersect2.size());
    EXPECT_DOUBLE_EQ(intersect1.area(), intersect2.area());

    EXPECT_EQ(intersect1.size(), 2);
    EXPECT_DOUBLE_EQ(intersect1.area(), 0.0);
  }

  {
    const Geo::Polygon<double> poly1{{
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
    }};

    const Geo::Polygon<double> poly2{{
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
    }};

    const auto intersect1 = Geo::intersection(poly1, poly2);
    const auto intersect2 = Geo::intersection(poly2, poly1);

    EXPECT_EQ(intersect1.size(), intersect2.size());
    EXPECT_DOUBLE_EQ(intersect1.area(), intersect2.area());

    ASSERT_EQ(intersect1.size(), poly2.size());
    for (size_t i = 0; i < intersect1.size(); ++i) {
      EXPECT_LE((intersect1[i] - poly2[i]).norm(), 1e-8);
      EXPECT_LE((intersect2[i] - poly2[i]).norm(), 1e-8);
    }

    EXPECT_DOUBLE_EQ(intersect1.area(), poly2.area());
  }

  {
    const Geo::Polygon<double> poly1{{
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {1.5, 0.5},
    }};

    const Geo::Polygon<double> poly2{{
        {0.0, 0.0},
        {0.0, 1.0},
        {1.0, 0.0},
        {1.0, 1.0},
    }};

    const auto intersect1 = Geo::intersection(poly1, poly2);
    const auto intersect2 = Geo::intersection(poly2, poly1);

    EXPECT_EQ(intersect1.size(), intersect2.size());
    EXPECT_DOUBLE_EQ(intersect1.area(), intersect2.area());

    ASSERT_EQ(intersect1.size(), poly2.size());
    for (size_t i = 0; i < intersect1.size(); ++i) {
      EXPECT_LE((intersect1[i] - poly2[i]).norm(), 1e-8);
      EXPECT_LE((intersect2[i] - poly2[i]).norm(), 1e-8);
    }

    EXPECT_DOUBLE_EQ(intersect1.area(), poly2.area());
  }
}

TEST(GeometryPolygon, IntersectRoundingBug) {
  {
    const Geo::Polygon<float> p1{{
        {1.8633899774085183f, 1.6666666709613775f},
        {1.6666666709613778f, 1.8633899774085183f},
        {1.6666666666666665f, 1.6666666666666665f},
    }};

    const Geo::Polygon<float> p2{{
        {1.8251895692700746f, 1.6666666666666665f},
        {1.8251895692700746f, 1.8633899774085183f},
        {1.6666666666666665f, 1.8633899774085183f},
        {1.6666666666666665f, 1.6666666666666665f},
    }};

    const auto intersect_12 = Geo::intersection(p1, p2);
    const auto intersect_21 = Geo::intersection(p2, p1);

#define EXPECT_EPS_EQ(a, b) EXPECT_NEAR(a, b, 1e-6f)  // NOLINT

    EXPECT_EPS_EQ(intersect_12.area(), intersect_21.area());

    ASSERT_EQ(intersect_12.size(), intersect_21.size());
    for (size_t i = 0; i < intersect_12.size(); ++i) {
      EXPECT_EPS_EQ(intersect_12[i](0), intersect_21[i](0));
      EXPECT_EPS_EQ(intersect_12[i](1), intersect_21[i](1));
    }

#undef EXPECT_EPS_EQ
  }

  {
    const Geo::Polygon<double> p1{{
        {1.8633899774085183, 1.6666666709613775},
        {1.6666666709613778, 1.8633899774085183},
        {1.6666666666666665, 1.6666666666666665},
    }};

    const Geo::Polygon<double> p2{{
        {1.8251895692700746, 1.6666666666666665},
        {1.8251895692700746, 1.8633899774085183},
        {1.6666666666666665, 1.8633899774085183},
        {1.6666666666666665, 1.6666666666666665},
    }};

    const auto intersect_12 = Geo::intersection(p1, p2);
    const auto intersect_21 = Geo::intersection(p2, p1);

#define EXPECT_EPS_EQ(a, b) EXPECT_NEAR(a, b, 1e-8)  // NOLINT

    EXPECT_EPS_EQ(intersect_12.area(), intersect_21.area());

    ASSERT_EQ(intersect_12.size(), intersect_21.size());
    for (size_t i = 0; i < intersect_12.size(); ++i) {
      EXPECT_EPS_EQ(intersect_12[i](0), intersect_21[i](0));
      EXPECT_EPS_EQ(intersect_12[i](1), intersect_21[i](1));
    }

#undef EXPECT_EPS_EQ
  }
}

TEST(GeometryPolygon, PointInPolygon) {
  Geo::Polygon<double> polygon({
      {0.0, 0.0},
      {1.0, 0.0},
      {0.0, 1.0},
      {1.0, 1.0},
  });

  EXPECT_TRUE(polygon.point_in_polygon(Zap::CellBased::Point<double>{0.5, 0.5}));
  EXPECT_FALSE(polygon.point_in_polygon(Zap::CellBased::Point<double>{1.5, 0.5}));
  EXPECT_FALSE(polygon.point_in_polygon(Zap::CellBased::Point<double>{0.5, 1.5}));
  EXPECT_FALSE(polygon.point_in_polygon(Zap::CellBased::Point<double>{1.5, 1.5}));

  EXPECT_TRUE(polygon.point_in_polygon(Zap::CellBased::Point<double>{0.0, 0.0}));
  EXPECT_TRUE(polygon.point_in_polygon(Zap::CellBased::Point<double>{0.1, 0.0}));
  EXPECT_TRUE(polygon.point_in_polygon(Zap::CellBased::Point<double>{0.0, 0.1}));
  EXPECT_TRUE(polygon.point_in_polygon(Zap::CellBased::Point<double>{0.1, 0.1}));
}
