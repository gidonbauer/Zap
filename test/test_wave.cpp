#include <gtest/gtest.h>

#include <numbers>

#include "CellBased/Wave.hpp"

using namespace Zap::CellBased;

// - Flux function for Burgers Equation ------------------------------------------------------------
template <typename Float>
[[nodiscard]] constexpr auto f(Float u) noexcept -> Float {
  return u * u / 2;
}

// -------------------------------------------------------------------------------------------------
TEST(Wave, XDirection) {
  const FullInterface<double, GridCoord<double>> interface{
      .left_value  = 1.0,
      .right_value = 0.0,
      .begin       = {0.0, 0.0},
      .end         = {0.0, 1.0},
  };
  const double dx = 0.1;
  const double dy = 0.2;
  const double dt = 0.03;

  const auto wave = normal_wave<X>(interface, dx, dy, dt);
  static_assert(std::is_same_v<std::remove_const_t<decltype(wave)>,
                               AxisAlignedWave<double, GridCoord<double>, X>>,
                "Expected an axis-aligned wave in x-direction.");

  EXPECT_LE((wave.begin - interface.begin).norm(), EPS<double>);
  EXPECT_LE((wave.end - interface.end).norm(), EPS<double>);

  EXPECT_TRUE(wave.is_right_going) << "Expected wave to move to the right but is not.";
  EXPECT_DOUBLE_EQ(wave.speed,
                   (f(interface.right_value) - f(interface.left_value)) /
                       (interface.right_value - interface.left_value))
      << "Wave speed does not match with Rankine-Hugoniot condition";

  // Igor::Info("wave.first_order_update  = {}", wave.first_order_update);
  // Igor::Info("wave.second_order_update = {}", wave.second_order_update);
}

// -------------------------------------------------------------------------------------------------
TEST(Wave, YDirection) {
  const FullInterface<double, GridCoord<double>> interface{
      .left_value  = 1.0,
      .right_value = 0.0,
      .begin       = {1.0, 0.0},
      .end         = {0.0, 0.0},
  };
  const double dx = 0.1;
  const double dy = 0.2;
  const double dt = 0.03;

  const auto wave = normal_wave<Y>(interface, dx, dy, dt);
  static_assert(std::is_same_v<std::remove_const_t<decltype(wave)>,
                               AxisAlignedWave<double, GridCoord<double>, Y>>,
                "Expected an axis-aligned wave in x-direction.");

  EXPECT_LE((wave.begin - interface.begin).norm(), EPS<double>);
  EXPECT_LE((wave.end - interface.end).norm(), EPS<double>);

  EXPECT_TRUE(wave.is_right_going) << "Expected wave to move to the right but is not.";
  EXPECT_DOUBLE_EQ(wave.speed,
                   (f(interface.right_value) - f(interface.left_value)) /
                       (interface.right_value - interface.left_value))
      << "Wave speed does not match with Rankine-Hugoniot condition";

  // Igor::Info("wave.first_order_update  = {}", wave.first_order_update);
  // Igor::Info("wave.second_order_update = {}", wave.second_order_update);
}

// -------------------------------------------------------------------------------------------------
TEST(Wave, FourtyFiveDegree) {
  const FullInterface<double, GridCoord<double>> interface{
      .left_value  = 1.0,
      .right_value = 0.0,
      .begin       = {1.0, 0.0},
      .end         = {0.0, 1.0},
  };
  const double dx = 0.1;
  const double dy = 0.2;
  const double dt = 0.03;

  const auto wave = normal_wave<FREE>(interface, dx, dy, dt);
  static_assert(
      std::is_same_v<std::remove_const_t<decltype(wave)>, FreeWave<double, GridCoord<double>>>,
      "Expected an axis-aligned wave in x-direction.");

  EXPECT_LE((wave.begin - interface.begin).norm(), EPS<double>);
  EXPECT_LE((wave.end - interface.end).norm(), EPS<double>);
  EXPECT_LE((wave.normal -
             GridCoord<double>{std::numbers::sqrt2 / (2 * dx), std::numbers::sqrt2 / (2 * dy)})
                .norm(),
            EPS<double>);

  EXPECT_TRUE(wave.is_right_going) << "Expected wave to move to the right but is not.";
  EXPECT_DOUBLE_EQ(wave.speed,
                   (std::cos(std::numbers::pi / 4) + std::sin(std::numbers::pi / 4)) *
                       (f(interface.right_value) - f(interface.left_value)) /
                       (interface.right_value - interface.left_value))
      << "Wave speed does not match with Rankine-Hugoniot condition";

  // Igor::Info("wave.first_order_update  = {}", wave.first_order_update);
  // Igor::Info("wave.second_order_update = {}", wave.second_order_update);
}

// -------------------------------------------------------------------------------------------------
TEST(Wave, TwohundredTwentyFiveDegree) {
  const FullInterface<double, GridCoord<double>> interface{
      .left_value  = 0.0,
      .right_value = 1.0,
      .begin       = {0.0, 1.0},
      .end         = {1.0, 0.0},
  };
  const double dx = 0.1;
  const double dy = 0.2;
  const double dt = 0.03;

  const auto wave = normal_wave<FREE>(interface, dx, dy, dt);
  static_assert(
      std::is_same_v<std::remove_const_t<decltype(wave)>, FreeWave<double, GridCoord<double>>>,
      "Expected an axis-aligned wave in x-direction.");

  EXPECT_LE((wave.begin - interface.begin).norm(), EPS<double>);
  EXPECT_LE((wave.end - interface.end).norm(), EPS<double>);
  EXPECT_LE((wave.normal -
             GridCoord<double>{-std::numbers::sqrt2 / (2 * dx), -std::numbers::sqrt2 / (2 * dy)})
                .norm(),
            EPS<double>);

  EXPECT_FALSE(wave.is_right_going) << "Expected wave to move to the left but is not.";
  EXPECT_DOUBLE_EQ(wave.speed,
                   (std::cos(5 * std::numbers::pi / 4) + std::sin(5 * std::numbers::pi / 4)) *
                       (f(interface.right_value) - f(interface.left_value)) /
                       (interface.right_value - interface.left_value))
      << "Wave speed does not match with Rankine-Hugoniot condition";
  EXPECT_LT(wave.speed, 0.0) << "Expected negative wave speed";

  // Igor::Info("wave.first_order_update  = {}", wave.first_order_update);
  // Igor::Info("wave.second_order_update = {}", wave.second_order_update);
}
