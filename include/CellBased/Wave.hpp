#ifndef ZAP_CELL_BASED_WAVE_HPP
#define ZAP_CELL_BASED_WAVE_HPP

#ifndef ZAP_NO_TANGENTIAL_CORRECTION
#define ZAP_TANGENTIAL_CORRECTION
#endif  //  ZAP_NO_TANGENTIAL_CORRECTION

#include <optional>

#include <fmt/format.h>

#include "CellBased/Definitions.hpp"
#include "CellBased/Geometry.hpp"
#include "CellBased/Interface.hpp"

namespace Zap::CellBased {

// - Scale vector when operating in grid coordinates -----------------------------------------------
template <typename PassiveFloat, Point2D_c PointType>
[[nodiscard]] constexpr auto
scale_if_grid_coord(const PointType& p, PassiveFloat dx, PassiveFloat dy) noexcept -> PointType {
  if constexpr (is_GridCoord_v<PointType>) {
    return PointType{p.x / dx, p.y / dy};
  } else {
    return p;
  }
}

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, Point2D_c PointType, Orientation orientation>
struct AxisAlignedWave {
  static_assert(orientation == X || orientation == Y, "Orientation must be `X` or `Y`.");

  ActiveFloat first_order_update;
  ActiveFloat second_order_update;

  ActiveFloat normal_speed;

  ActiveFloat tangent_speed;

  PointType begin;
  PointType end;
};

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, Point2D_c PointType>
struct FreeWave {
  ActiveFloat first_order_update;
  ActiveFloat second_order_update;

  ActiveFloat normal_speed;
  PointType normal;

  ActiveFloat tangent_speed;
  PointType tangent;

  PointType begin;
  PointType end;
};

// -------------------------------------------------------------------------------------------------
// Reference:
// - normal direction:     https://github.com/clawpack/riemann/blob/master/src/rpn2_burgers.f90
// - tangential direction: https://github.com/clawpack/riemann/blob/master/src/rpt2_burgers.f90
template <Orientation orientation, typename ActiveFloat, typename PassiveFloat, Point2D_c PointType>
requires(orientation == X || orientation == Y || orientation == FREE)
[[nodiscard]] constexpr auto calc_wave(const FullInterface<ActiveFloat, PointType>& interface,
                                       PassiveFloat dx,
                                       PassiveFloat dy,
                                       [[maybe_unused]] ActiveFloat dt) noexcept
    -> std::optional<std::conditional_t<orientation == FREE,
                                        FreeWave<ActiveFloat, PointType>,
                                        AxisAlignedWave<ActiveFloat, PointType, orientation>>> {
  if ((interface.end - interface.begin).norm() < EPS<PassiveFloat>) { return std::nullopt; }

  const ActiveFloat wave = interface.right_value - interface.left_value;  // As in LeVeque Book

  // - Axis aligned wave ---------------------------------------------------------------------------
  if constexpr (orientation == X || orientation == Y) {
    const ActiveFloat wave_speed = (interface.left_value + interface.right_value) / 2;

    // TODO: wave in second_order_update can be replaced by limited version
    // TODO: how is the second_order_update distributed when multiple cells are swept over?
    // TODO: maybe add entropy fix for transsonic
    if constexpr (orientation == X) {
      IGOR_ASSERT(approx_eq(interface.begin.x, interface.end.x),
                  "Interface is not aligned to x-axis, begin={}, end={}",
                  interface.begin,
                  interface.end);

      return std::make_optional(AxisAlignedWave<ActiveFloat, PointType, X>{
          .first_order_update = sign(wave_speed) * wave,
#ifdef ZAP_2ND_ORDER_CORRECTION
          .second_order_update =
              0.5 * std::abs(wave_speed) * (1 - dt / dx * std::abs(wave_speed)) * wave,
#else
          .second_order_update = 0,
#endif  //  ZAP_2ND_ORDER_CORRECTION

          .normal_speed = wave_speed,

#ifdef ZAP_TANGENTIAL_CORRECTION
          .tangent_speed = wave_speed,
#else
          .tangent_speed = 0,
#endif  //  ZAP_TANGENTIAL_CORRECTION

          .begin = interface.begin,
          .end   = interface.end,
      });
    } else {
      IGOR_ASSERT(approx_eq(interface.begin.y, interface.end.y),
                  "Interface is not aligned to y-axis, begin={}, end={}",
                  interface.begin,
                  interface.end);

      return std::make_optional(AxisAlignedWave<ActiveFloat, PointType, Y>{
          .first_order_update = sign(wave_speed) * wave,
#ifdef ZAP_2ND_ORDER_CORRECTION
          .second_order_update =
              0.5 * std::abs(wave_speed) * (1 - dt / dy * std::abs(wave_speed)) * wave,
#else
          .second_order_update = 0,
#endif  //  ZAP_2ND_ORDER_CORRECTION

          .normal_speed = wave_speed,

#ifdef ZAP_TANGENTIAL_CORRECTION
          .tangent_speed = wave_speed,
#else
          .tangent_speed = 0,
#endif  //  ZAP_TANGENTIAL_CORRECTION

          .begin = interface.begin,
          .end   = interface.end,
      });
    }
  }
  // - Free wave -----------------------------------------------------------------------------------
  else {
    const PointType tangent_vector = (interface.end - interface.begin).normalized();
    const PointType normal_vector  = {tangent_vector.y, -tangent_vector.x};

    // Rotate PDE if interface does not align with an axis:
    //  tangent.y                                     = cos(interface_angle)
    //  -sign(tangent.x) * std::sqrt(1 - tangent.y^2) = sin(interface_angle)
    const ActiveFloat cos_angle = tangent_vector.y;
    const ActiveFloat sin_angle =
        -sign(tangent_vector.x) * std::sqrt(1 - tangent_vector.y * tangent_vector.y);
    const ActiveFloat u_mid         = (interface.left_value + interface.right_value) / 2;
    const ActiveFloat normal_speed  = u_mid * (cos_angle + sin_angle);
    const ActiveFloat tangent_speed = u_mid * (-sin_angle + cos_angle);

    return std::make_optional(FreeWave<ActiveFloat, PointType>{
        .first_order_update =
            sign(normal_speed) * wave,  // TODO: Double check if sign(normal_speed) is correct
        .second_order_update = 0,       // TODO: How does the second_order_update work here?

        .normal_speed = normal_speed,
        .normal       = scale_if_grid_coord(normal_vector, dx, dy),

#ifdef ZAP_TANGENTIAL_CORRECTION
        .tangent_speed = tangent_speed,
#else
        .tangent_speed = 0,
#endif  //  ZAP_TANGENTIAL_CORRECTION
        .tangent = scale_if_grid_coord(tangent_vector, dx, dy),

        .begin = interface.begin,
        .end   = interface.end,
    });
  }
}

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, typename PointType>
[[nodiscard]] constexpr auto calc_wave_polygon(const FreeWave<ActiveFloat, PointType>& wave,
                                               const ActiveFloat& dt) noexcept
    -> Geometry::Polygon<PointType> {
  return Geometry::Polygon<PointType>({
      wave.begin,
      wave.begin + dt * wave.normal_speed * wave.normal + dt * wave.tangent_speed * wave.tangent,
      wave.end,
      wave.end + dt * wave.normal_speed * wave.normal + dt * wave.tangent_speed * wave.tangent,
  });
}

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, typename PassiveFloat, typename PointType, Orientation orientation>
[[nodiscard]] constexpr auto
calc_wave_polygon(const AxisAlignedWave<ActiveFloat, PointType, orientation>& wave,
                  PassiveFloat dx,
                  PassiveFloat dy,
                  const ActiveFloat& dt) noexcept -> Geometry::Polygon<PointType> {
  const PointType normal =
      scale_if_grid_coord(PointType{static_cast<PassiveFloat>(orientation == X),
                                    static_cast<PassiveFloat>(orientation == Y)},
                          dx,
                          dy);

  const PointType tangent =
      scale_if_grid_coord(PointType{-static_cast<PassiveFloat>(orientation == Y),
                                    static_cast<PassiveFloat>(orientation == X)},
                          dx,
                          dy);

  return Geometry::Polygon<PointType>({
      wave.begin,
      wave.begin + dt * wave.normal_speed * normal + dt * wave.tangent_speed * tangent,
      wave.end,
      wave.end + dt * wave.normal_speed * normal + dt * wave.tangent_speed * tangent,
  });
}

}  // namespace Zap::CellBased

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat,
          Zap::CellBased::Point2D_c PointType,
          Zap::CellBased::Orientation orientation>
struct fmt::formatter<Zap::CellBased::AxisAlignedWave<ActiveFloat, PointType, orientation>> {
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) const noexcept {
    return ctx.begin();
  }
  template <typename FormatContext>
  constexpr auto
  format(const Zap::CellBased::AxisAlignedWave<ActiveFloat, PointType, orientation>& wave,
         FormatContext& ctx) const noexcept {
    return fmt::format_to(ctx.out(),
                          "{{\n"
                          "  .first_order_update = {}\n"
                          "  .second_order_update = {}\n"
                          "  .normal_speed = {}\n"
                          "  .tangent_speed = {}\n"
                          "  .begin = {}\n"
                          "  .end = {}\n"
                          "  .orientation = {}\n"
                          "}}",
                          wave.first_order_update,
                          wave.second_order_update,
                          wave.normal_speed,
                          wave.tangent_speed,
                          wave.begin,
                          wave.end,
                          orientation == Zap::CellBased::X ? 'X' : 'Y');
  }
};

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, Zap::CellBased::Point2D_c PointType>
struct fmt::formatter<Zap::CellBased::FreeWave<ActiveFloat, PointType>> {
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) const noexcept {
    return ctx.begin();
  }
  template <typename FormatContext>
  constexpr auto format(const Zap::CellBased::FreeWave<ActiveFloat, PointType>& wave,
                        FormatContext& ctx) const noexcept {
    return fmt::format_to(ctx.out(),
                          "{{\n"
                          "  .first_order_update = {}\n"
                          "  .second_order_update = {}\n"
                          "  .normal_speed = {}\n"
                          "  .normal = {}\n"
                          "  .tangent_speed = {}\n"
                          "  .tangent = {}\n"
                          "  .begin = {}\n"
                          "  .end = {}\n"
                          "}}",
                          wave.first_order_update,
                          wave.second_order_update,
                          wave.normal_speed,
                          wave.normal,
                          wave.tangent_speed,
                          wave.tangent,
                          wave.begin,
                          wave.end,
                          wave.normal);
  }
};

#endif  // ZAP_CELL_BASED_WAVE_HPP
