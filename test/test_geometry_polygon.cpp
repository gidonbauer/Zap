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
