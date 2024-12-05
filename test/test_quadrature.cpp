#include <gtest/gtest.h>

#include <numbers>

#include "CellBased/Definitions.hpp"
#include "CellBased/Grid.hpp"
#include "CellBased/Quadrature.hpp"

using namespace Zap::CellBased;

using Float        = double;
using PointType    = GenCoord<Float>;
constexpr size_t N = 15UZ;

// -------------------------------------------------------------------------------------------------
TEST(Quadrature, UnitSquare) {
  Geometry::Polygon<PointType> polygon({{-1.0, -1.0}, {-1.0, 1.0}, {1.0, 1.0}, {1.0, -1.0}});
  auto f = [](Float x, Float y) { return x * x * y * y; };

  const auto integral = quadrature<N>(f, polygon);
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

  const auto integral = quadrature<N>(f, polygon);
  static_assert(std::is_same_v<decltype(integral), const double>);

  EXPECT_NEAR(integral, 4.0 / 9.0, EPS<Float>());
}

// -------------------------------------------------------------------------------------------------
TEST(Quadrature, TranslatedScaledSquare) {
  Geometry::Polygon<PointType> polygon({{1.0, 3.0}, {1.0, 7.0}, {10.0, 7.0}, {10.0, 3.0}});
  auto f = [](Float x, Float y) -> Float { return std::sin(x) + std::exp(-y); };

  const auto integral = quadrature<N>(f, polygon);
  static_assert(std::is_same_v<decltype(integral), const double>);
  const Float expected_integral =
      -4 * std::cos(10.0) + 4 * std::cos(1.0) - 9 * (std::exp(-7.0) - std::exp(-3.0));

  EXPECT_NEAR(integral, expected_integral, EPS<Float>());
}

// -------------------------------------------------------------------------------------------------
TEST(Quadrature, Quadrilateral) {
  const Geometry::Polygon<PointType> polygon({{1.0, 0.5}, {5.0, 0.5}, {5.0, 2.0}, {1.0, 3.0}});
  auto f = [](Float x, Float y) -> Float { return std::sin(x) + std::exp(-y); };

  const auto integral = quadrature<N>(f, polygon);
  static_assert(std::is_same_v<decltype(integral), const double>);
  // From Wolfram Alpha
  const Float expected_integral =
      4.0 / std::pow(std::numbers::e_v<Float>, 3.0) -
      4.0 / std::pow(std::numbers::e_v<Float>, 2.0) + 4.0 / std::sqrt(std::numbers::e_v<Float>) +
      1.0 / 4.0 * (std::sin(1.0) - std::sin(5.0) + 10.0 * std::cos(1.0) - 6.0 * std::cos(5.0));

  EXPECT_NEAR(integral, expected_integral, EPS<Float>());
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

    const auto integral = quadrature<N>(f, polygon);
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

    const auto integral = quadrature<N>(f, polygon);
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

    const auto integral = quadrature<N>(f, polygon);
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

    const auto integral = quadrature<N>(f, polygon);
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
      const auto integral = quadrature<N>(f, polygon);
      const auto area     = polygon.area();
      EXPECT_NEAR(integral, area, EPS<Float>());
    } else {
      {
        const auto polygon  = cell.get_cut_left_polygon<SIM_C>();
        const auto integral = quadrature<N>(f, polygon);
        const auto area     = polygon.area();
        EXPECT_NEAR(integral, area, EPS<Float>());
      }
      {
        const auto polygon  = cell.get_cut_right_polygon<SIM_C>();
        const auto integral = quadrature<N>(f, polygon);
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

    const auto integral = quadrature<N>(f, polygon);
    // From scipy dblquad
    constexpr auto expected_integral = static_cast<Float>(0.4975755788114332);

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

    const auto integral = quadrature<N>(f, polygon);
    // From scipy dblquad
    constexpr auto expected_integral = static_cast<Float>(19.6588273896563095);

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

    const auto integral = quadrature<N>(f, polygon);
    // From scipy dblquad
    constexpr auto expected_integral = static_cast<Float>(64.2509552141199833);

    EXPECT_NEAR(integral, expected_integral, EPS<Float>());
  }
}
