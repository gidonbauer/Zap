#include <gtest/gtest.h>

#include <numbers>

// #define ZAP_QUAD_2
#include "CellBased/Definitions.hpp"
#include "CellBased/Grid.hpp"
#include "CellBased/Quadrature.hpp"

using namespace Zap::CellBased;

using Float     = double;
using PointType = GenCoord<Float>;

// -------------------------------------------------------------------------------------------------
TEST(Quadrature, UnitSquare) {
  Geometry::Polygon<PointType> polygon({{-1.0, -1.0}, {-1.0, 1.0}, {1.0, 1.0}, {1.0, -1.0}});
  auto f = [](Float x, Float y) { return x * x * y * y; };

  const auto integral = quadrature(f, polygon);
  static_assert(std::is_same_v<decltype(integral), const double>);

  EXPECT_NEAR(integral, 4.0 / 9.0, EPS<Float>());
}

// -------------------------------------------------------------------------------------------------
TEST(Quadrature, TranslatedUnitSquare) {
  Geometry::Polygon<PointType> polygon({{1.0, 1.0}, {1.0, 3.0}, {3.0, 3.0}, {3.0, 1.0}});
  auto f = [](Float x, Float y) {
    x -= 2;
    y -= 2;
    return x * x * y * y;
  };

  const auto integral = quadrature(f, polygon);
  static_assert(std::is_same_v<decltype(integral), const double>);

  EXPECT_NEAR(integral, 4.0 / 9.0, EPS<Float>());
}

// -------------------------------------------------------------------------------------------------
TEST(Quadrature, FourPointsArea) {
  {
    Geometry::Polygon<PointType> polygon({
        {-1.5, -1.5},
        {2.0, -3.0},
        {2.5, 2.5},
        {-2.0, 2.0},
    });
    auto f = [](Float /*x*/, Float /*y*/) { return static_cast<Float>(1); };

    const auto integral = quadrature(f, polygon);
    static_assert(std::is_same_v<decltype(integral), const double>);

    EXPECT_NEAR(integral, polygon.area(), EPS<Float>());
  }

  {
    Geometry::Polygon<PointType> polygon({
        {0.5, 1.5},
        {4.0, 0.0},
        {4.5, 5.5},
        {0.0, 5.0},
    });
    auto f = [](Float /*x*/, Float /*y*/) { return static_cast<Float>(1); };

    const auto integral = quadrature(f, polygon);
    static_assert(std::is_same_v<decltype(integral), const double>);

    EXPECT_NEAR(integral, polygon.area(), EPS<Float>());
  }
}

// -------------------------------------------------------------------------------------------------
TEST(Quadrature, TriangleArea) {
  {
    Geometry::Polygon<PointType> polygon({
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 1.0},
    });
    auto f = [](Float /*x*/, Float /*y*/) { return static_cast<Float>(1); };

    const auto integral = quadrature(f, polygon);
    static_assert(std::is_same_v<decltype(integral), const double>);

    EXPECT_NEAR(integral, polygon.area(), EPS<Float>());
  }
}

// -------------------------------------------------------------------------------------------------
TEST(Quadrature, PentagonArea) {
  {
    Geometry::Polygon<PointType> polygon({
        {-1.0, -1.0},
        {1.0, -1.0},
        {-1.0, 1.0},
        {1.0, 1.0},
        {0.0, 2.0},
    });
    auto f = [](Float /*x*/, Float /*y*/) { return static_cast<Float>(1); };

    const auto integral = quadrature(f, polygon);
    static_assert(std::is_same_v<decltype(integral), const double>);

    EXPECT_NEAR(integral, polygon.area(), EPS<Float>());
  }
}

// -------------------------------------------------------------------------------------------------
TEST(Quadrature, GridCellArea) {
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

  auto f = [](Float /*x*/, Float /*y*/) -> Float { return 1; };
  for (const auto& cell : grid) {
    if (cell.is_cartesian()) {
      const auto polygon  = cell.get_cartesian_polygon<SIM_C>();
      const auto integral = quadrature(f, polygon);
      const auto area     = polygon.area();
      EXPECT_NEAR(integral, area, EPS<Float>());
    } else {
      {
        const auto polygon  = cell.get_cut_left_polygon<SIM_C>();
        const auto integral = quadrature(f, polygon);
        const auto area     = polygon.area();
        EXPECT_NEAR(integral, area, EPS<Float>());
      }
      {
        const auto polygon  = cell.get_cut_right_polygon<SIM_C>();
        const auto integral = quadrature(f, polygon);
        const auto area     = polygon.area();
        EXPECT_NEAR(integral, area, EPS<Float>());
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------
TEST(Quadrature, IntegrateOverPentagon) {
  {
    Geometry::Polygon<PointType> polygon({
        {-2.5, -1.0},
        {-2.0, -1.3333333333333333},
        {-2.0, 2.0},
    });

    auto f = [](Float x, Float y) -> Float { return std::sin(x) + std::exp(-y); };

    const auto integral = quadrature(f, polygon);
    // From scipy dblquad
    constexpr auto expected_integral = static_cast<Float>(7.229971054944182);

    EXPECT_NEAR(integral, expected_integral, EPS<Float>());
  }

  {
    Geometry::Polygon<PointType> polygon({
        {-2.5, -1.0},
        {0.5, -3.0},
        {0.5, 2.5},
        {-2.0, 2.0},
    });

    auto f = [](Float x, Float y) -> Float { return std::sin(x) + std::exp(-y); };

    const auto integral = quadrature(f, polygon);
    // From scipy dblquad
    constexpr auto expected_integral = static_cast<Float>(34.710327710602606);

    EXPECT_NEAR(integral, expected_integral, EPS<Float>());
  }

  {
    Geometry::Polygon<PointType> polygon({
        {-2.5, -1.0},
        {0.5, -3.0},
        {4.0, -1.5},
        {3.0, 3.0},
        {-2.0, 2.0},
    });

    auto f = [](Float x, Float y) -> Float { return std::sin(x) + std::exp(-y); };

    const auto integral = quadrature(f, polygon);
    // From scipy dblquad
    constexpr auto expected_integral = static_cast<Float>(37.7994083154342064);

    EXPECT_NEAR(integral, expected_integral, EPS<Float>());
  }
}
