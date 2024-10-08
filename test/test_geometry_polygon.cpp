#include <gtest/gtest.h>

#include "CellBased/Geometry.hpp"
namespace Geo = Zap::CellBased::Geometry;
using Zap::CellBased::SmallVector;

template <typename T>
struct P2D {
  T x;
  T y;

  [[nodiscard]] static constexpr auto Zero() noexcept -> P2D { return {.x = 0, .y = 0}; }

  [[nodiscard]] constexpr auto norm() const noexcept -> T { return std::sqrt(x * x + y * y); }
  [[nodiscard]] constexpr auto normalized() const noexcept -> P2D { return *this / norm(); }

  [[nodiscard]] constexpr auto dot(const P2D& other) const noexcept -> T {
    return x * other.x + y * other.y;
  }

  [[nodiscard]] friend constexpr auto operator+(const P2D& lhs, const P2D& rhs) noexcept -> P2D {
    return {.x = lhs.x + rhs.x, .y = lhs.y + rhs.y};
  }

  [[nodiscard]] friend constexpr auto operator-(const P2D& lhs, const P2D& rhs) noexcept -> P2D {
    return {.x = lhs.x - rhs.x, .y = lhs.y - rhs.y};
  }

  [[nodiscard]] friend constexpr auto operator*(const P2D& lhs, T scalar) noexcept -> P2D {
    return {.x = lhs.x * scalar, .y = lhs.y * scalar};
  }

  [[nodiscard]] friend constexpr auto operator*(T scalar, const P2D& rhs) noexcept -> P2D {
    return {.x = scalar * rhs.x, .y = scalar * rhs.y};
  }

  [[nodiscard]] friend constexpr auto operator/(const P2D& lhs, T scalar) noexcept -> P2D {
    return {.x = lhs.x / scalar, .y = lhs.y / scalar};
  }

  constexpr auto operator+=(const P2D& other) noexcept -> P2D& {
    x += other.x;
    y += other.y;
    return *this;
  }

  constexpr auto operator-=(const P2D& other) noexcept -> P2D& {
    x -= other.x;
    y -= other.y;
    return *this;
  }

  constexpr auto operator*=(T scalar) noexcept -> P2D& {
    x *= scalar;
    y *= scalar;
    return *this;
  }

  constexpr auto operator/=(T scalar) noexcept -> P2D& {
    x /= scalar;
    y /= scalar;
    return *this;
  }
};

TEST(GeometryPolygon, AddPoint) {
  using Float = double;

  Geo::Polygon<P2D<Float>> polygon{};
  EXPECT_EQ(polygon.size(), 0);
  EXPECT_DOUBLE_EQ(polygon.area(), 0.0);

  polygon.add_point(P2D<Float>{.x = 0.0, .y = 0.0});
  EXPECT_EQ(polygon.size(), 1);
  EXPECT_DOUBLE_EQ(polygon.area(), 0.0);

  polygon.add_point(P2D<Float>{.x = 0.0, .y = 1.0});
  EXPECT_EQ(polygon.size(), 2);
  EXPECT_DOUBLE_EQ(polygon.area(), 0.0);

  polygon.add_point(P2D<Float>{.x = 1.0, .y = 0.0});
  EXPECT_EQ(polygon.size(), 3);
  EXPECT_DOUBLE_EQ(polygon.area(), 0.5);

  polygon.add_point(P2D<Float>{.x = 1.0, .y = 1.0});
  EXPECT_EQ(polygon.size(), 4);
  EXPECT_DOUBLE_EQ(polygon.area(), 1.0);

  polygon.add_point(P2D<Float>{.x = 1.0, .y = 1.0});
  EXPECT_EQ(polygon.size(), 4);
  EXPECT_DOUBLE_EQ(polygon.area(), 1.0);
}

TEST(GeometryPolygon, Area) {
  {
    using Float = float;
    const SmallVector<P2D<Float>> points{
        {.x = 0.0, .y = 0.0},
        {.x = 0.0, .y = 1.0},
        {.x = 1.0, .y = 0.0},
        {.x = 1.0, .y = 1.0},
        {.x = 1.0, .y = 1.0},
    };

    const Geo::Polygon polygon{points};

    EXPECT_EQ(polygon.size(), 4);
    EXPECT_DOUBLE_EQ(polygon.area(), 1.0);
  }

  {
    using Float = double;
    const SmallVector<P2D<Float>> points{
        {.x = 0.0, .y = 0.0},
        {.x = 0.0, .y = 1.0},
        {.x = 1.0, .y = 0.0},
        {.x = 1.0, .y = 1.0},
        {.x = 1.0, .y = 1.0},
    };

    const Geo::Polygon polygon{points};

    EXPECT_EQ(polygon.size(), 4);
    EXPECT_DOUBLE_EQ(polygon.area(), 1.0);
  }

  {
    using Float = double;
    const SmallVector<P2D<Float>> points{
        {.x = 0.0, .y = 0.0},
        {.x = 0.0, .y = 1.0},
        {.x = 1.0, .y = 1.0},
    };

    const Geo::Polygon polygon{points};

    EXPECT_DOUBLE_EQ(polygon.area(), 0.5);
  }

  {
    using Float = double;
    const SmallVector<P2D<Float>> points{
        {.x = 0.0, .y = 0.0},
        {.x = 0.0, .y = 1.0},
        {.x = 1.0, .y = 0.0},
        {.x = 1.0, .y = 1.0},
        {.x = 1.5, .y = 0.5},
    };

    const Geo::Polygon polygon{points};

    EXPECT_DOUBLE_EQ(polygon.area(), 1.25);
  }
}

TEST(GeometryPolygon, Intersect) {
  {
    const Geo::Polygon<P2D<double>> poly1{{
        {.x = 0.0, .y = 0.0},
        {.x = 0.0, .y = 1.0},
        {.x = 1.0, .y = 0.0},
        {.x = 1.0, .y = 1.0},
    }};

    const Geo::Polygon<P2D<double>> poly2{{}};

    const auto intersect1 = Geo::intersection(poly1, poly2);
    const auto intersect2 = Geo::intersection(poly2, poly1);

    EXPECT_EQ(intersect1.size(), intersect2.size());
    EXPECT_DOUBLE_EQ(intersect1.area(), intersect2.area());

    EXPECT_EQ(intersect1.size(), 0);
    EXPECT_DOUBLE_EQ(intersect1.area(), 0.0);
  }

  {
    const Geo::Polygon<P2D<double>> poly1{{
        {.x = 0.0, .y = 0.0},
        {.x = 0.0, .y = 1.0},
        {.x = 1.0, .y = 0.0},
        {.x = 1.0, .y = 1.0},
    }};

    const Geo::Polygon<P2D<double>> poly2{{
        {.x = 2.0, .y = 0.0},
        {.x = 2.0, .y = 1.0},
        {.x = 3.0, .y = 0.0},
    }};

    const auto intersect1 = Geo::intersection(poly1, poly2);
    const auto intersect2 = Geo::intersection(poly2, poly1);

    EXPECT_EQ(intersect1.size(), intersect2.size());
    EXPECT_DOUBLE_EQ(intersect1.area(), intersect2.area());

    EXPECT_EQ(intersect1.size(), 0);
    EXPECT_DOUBLE_EQ(intersect1.area(), 0.0);
  }

  {
    const Geo::Polygon<P2D<double>> poly1{{
        {.x = 0.0, .y = 0.0},
        {.x = 0.0, .y = 1.0},
        {.x = 1.0, .y = 0.0},
        {.x = 1.0, .y = 1.0},
    }};

    const Geo::Polygon<P2D<double>> poly2{{
        {.x = 1.0, .y = 0.0},
        {.x = 1.0, .y = 1.0},
        {.x = 3.0, .y = 0.0},
    }};

    const auto intersect1 = Geo::intersection(poly1, poly2);
    const auto intersect2 = Geo::intersection(poly2, poly1);

    EXPECT_EQ(intersect1.size(), intersect2.size());
    EXPECT_DOUBLE_EQ(intersect1.area(), intersect2.area());

    EXPECT_EQ(intersect1.size(), 2);
    EXPECT_DOUBLE_EQ(intersect1.area(), 0.0);
  }

  {
    const Geo::Polygon<P2D<double>> poly1{{
        {.x = 0.0, .y = 0.0},
        {.x = 0.0, .y = 1.0},
        {.x = 1.0, .y = 0.0},
        {.x = 1.0, .y = 1.0},
    }};

    const Geo::Polygon<P2D<double>> poly2{{
        {.x = 0.0, .y = 0.0},
        {.x = 0.0, .y = 1.0},
        {.x = 1.0, .y = 0.0},
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
    const Geo::Polygon<P2D<double>> poly1{{
        {.x = 0.0, .y = 0.0},
        {.x = 0.0, .y = 1.0},
        {.x = 1.0, .y = 0.0},
        {.x = 1.0, .y = 1.0},
        {.x = 1.5, .y = 0.5},
    }};

    const Geo::Polygon<P2D<double>> poly2{{
        {.x = 0.0, .y = 0.0},
        {.x = 0.0, .y = 1.0},
        {.x = 1.0, .y = 0.0},
        {.x = 1.0, .y = 1.0},
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

// TEST(GeometryPolygon, IntersectRoundingBugFloat) {
//   // TODO: float tolerance seems to be wrong for this bug
//   const Geo::Polygon<P2D<float>> p1{{
//       {.x = 1.8633899774085183f, .y = 1.6666666709613775f},
//       {.x = 1.6666666709613778f, .y = 1.8633899774085183f},
//       {.x = 1.6666666666666665f, .y = 1.6666666666666665f},
//   }};
//
//   const Geo::Polygon<P2D<float>> p2{{
//       {.x = 1.8251895692700746f, .y = 1.6666666666666665f},
//       {.x = 1.8251895692700746f, .y = 1.8633899774085183f},
//       {.x = 1.6666666666666665f, .y = 1.8633899774085183f},
//       {.x = 1.6666666666666665f, .y = 1.6666666666666665f},
//   }};
//
//   const auto intersect_12 = Geo::intersection(p1, p2);
//   const auto intersect_21 = Geo::intersection(p2, p1);
//
// #define EXPECT_EPS_EQ(a, b) EXPECT_NEAR(a, b, 1e-6f)  // NOLINT
//
//   EXPECT_EPS_EQ(intersect_12.area(), intersect_21.area());
//
//   ASSERT_EQ(intersect_12.size(), intersect_21.size());
//   for (size_t i = 0; i < intersect_12.size(); ++i) {
//     EXPECT_EPS_EQ(intersect_12[i].x, intersect_21[i].x);
//     EXPECT_EPS_EQ(intersect_12[i].y, intersect_21[i].y);
//   }
//
// #undef EXPECT_EPS_EQ
// }

TEST(GeometryPolygon, IntersectRoundingBugDouble) {
  const Geo::Polygon<P2D<double>> p1{{
      {.x = 1.8633899774085183, .y = 1.6666666709613775},
      {.x = 1.6666666709613778, .y = 1.8633899774085183},
      {.x = 1.6666666666666665, .y = 1.6666666666666665},
  }};

  const Geo::Polygon<P2D<double>> p2{{
      {.x = 1.8251895692700746, .y = 1.6666666666666665},
      {.x = 1.8251895692700746, .y = 1.8633899774085183},
      {.x = 1.6666666666666665, .y = 1.8633899774085183},
      {.x = 1.6666666666666665, .y = 1.6666666666666665},
  }};

  const auto intersect_12 = Geo::intersection(p1, p2);
  const auto intersect_21 = Geo::intersection(p2, p1);

#define EXPECT_EPS_EQ(a, b) EXPECT_NEAR(a, b, 1e-8)  // NOLINT

  EXPECT_EPS_EQ(intersect_12.area(), intersect_21.area());

  ASSERT_EQ(intersect_12.size(), intersect_21.size());
  for (size_t i = 0; i < intersect_12.size(); ++i) {
    EXPECT_EPS_EQ(intersect_12[i].x, intersect_21[i].x);
    EXPECT_EPS_EQ(intersect_12[i].y, intersect_21[i].y);
  }

#undef EXPECT_EPS_EQ
}

TEST(GeometryPolygon, PointInPolygon) {
  Geo::Polygon<P2D<double>> polygon({
      {.x = 0.0, .y = 0.0},
      {.x = 1.0, .y = 0.0},
      {.x = 0.0, .y = 1.0},
      {.x = 1.0, .y = 1.0},
  });

  EXPECT_TRUE(polygon.point_in_polygon(P2D<double>{0.5, 0.5}));
  EXPECT_FALSE(polygon.point_in_polygon(P2D<double>{1.5, 0.5}));
  EXPECT_FALSE(polygon.point_in_polygon(P2D<double>{0.5, 1.5}));
  EXPECT_FALSE(polygon.point_in_polygon(P2D<double>{1.5, 1.5}));

  EXPECT_TRUE(polygon.point_in_polygon(P2D<double>{0.0, 0.0}));
  EXPECT_TRUE(polygon.point_in_polygon(P2D<double>{0.1, 0.0}));
  EXPECT_TRUE(polygon.point_in_polygon(P2D<double>{0.0, 0.1}));
  EXPECT_TRUE(polygon.point_in_polygon(P2D<double>{0.1, 0.1}));
}

TEST(GeometryPolygon, ScalingAndRelativeSizes) {
  Geo::Polygon<P2D<double>> polygon_A({
      {.x = 0.2, .y = 0.4},
      {.x = 0.4, .y = 0.2},
      {.x = 0.55, .y = 0.6},
  });

  Geo::Polygon<P2D<double>> polygon_B({
      {.x = 0.8, .y = 0.35},
      {.x = 0.75, .y = 0.8},
      {.x = 0.25, .y = 0.8},
      {.x = 0.4, .y = 0.45},
  });

  const auto intersect_AB = Geo::intersection(polygon_A, polygon_B);

  const double scale_x = 1.75;
  const double scale_y = 0.55;

  Geo::Polygon<P2D<double>> scaled_polygon_A;
  for (const auto& p : polygon_A.points()) {
    scaled_polygon_A.add_point({.x = scale_x * p.x, .y = scale_y * p.y});
  }

  Geo::Polygon<P2D<double>> scaled_polygon_B;
  for (const auto& p : polygon_B.points()) {
    scaled_polygon_B.add_point({.x = scale_x * p.x, .y = scale_y * p.y});
  }

  const auto scaled_intersect_AB = Geo::intersection(scaled_polygon_A, scaled_polygon_B);

  EXPECT_DOUBLE_EQ(scaled_polygon_A.area(), scale_x * scale_y * polygon_A.area());

  EXPECT_DOUBLE_EQ(scaled_polygon_B.area(), scale_x * scale_y * polygon_B.area());

  EXPECT_DOUBLE_EQ(polygon_A.area() / polygon_B.area(),
                   scaled_polygon_A.area() / scaled_polygon_B.area());

  EXPECT_DOUBLE_EQ(scaled_polygon_A.area() / polygon_A.area(),
                   scaled_polygon_B.area() / polygon_B.area());

  EXPECT_NEAR(scaled_intersect_AB.area(),
              scale_x * scale_y * intersect_AB.area(),
              Zap::CellBased::EPS<double>);

  EXPECT_NEAR(intersect_AB.area() / polygon_A.area(),
              scaled_intersect_AB.area() / scaled_polygon_A.area(),
              Zap::CellBased::EPS<double>);

  EXPECT_NEAR(intersect_AB.area() / polygon_B.area(),
              scaled_intersect_AB.area() / scaled_polygon_B.area(),
              Zap::CellBased::EPS<double>);
}
