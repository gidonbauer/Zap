#ifndef ZAP_CELL_BASED_WAVE_HPP
#define ZAP_CELL_BASED_WAVE_HPP

#include <fmt/format.h>

#include "CellBased/Definitions.hpp"
#include "CellBased/Geometry.hpp"
#include "CellBased/Interface.hpp"

namespace Zap::CellBased {

// template <typename ActiveFloat>
// [[nodiscard]] constexpr auto minmod(const ActiveFloat& a, const ActiveFloat& b) noexcept
//     -> ActiveFloat {
//   if (a * b <= 0.0) { return 0.0; }
//   if (std::abs(a) < std::abs(b)) { return a; }
//   return b;
// }

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
requires(orientation == X || orientation == Y)
struct AxisAlignedWave {
  ActiveFloat first_order_update;
  ActiveFloat second_order_update;
  ActiveFloat speed;
  bool is_right_going;
  PointType begin;
  PointType end;
};

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, Point2D_c PointType>
struct FreeWave {
  ActiveFloat first_order_update;
  ActiveFloat second_order_update;
  ActiveFloat speed;
  bool is_right_going;
  PointType begin;
  PointType end;
  PointType normal;
};

// -------------------------------------------------------------------------------------------------
// Reference:
// - normal direction:     https://github.com/clawpack/riemann/blob/master/src/rpn2_burgers.f90
// - tangential direction: https://github.com/clawpack/riemann/blob/master/src/rpt2_burgers.f90
template <Orientation orientation, typename ActiveFloat, typename PassiveFloat, Point2D_c PointType>
requires(orientation == X || orientation == Y || orientation == FREE)
[[nodiscard]] constexpr auto normal_wave(const FullInterface<ActiveFloat, PointType>& interface,
                                         PassiveFloat dx,
                                         PassiveFloat dy,
                                         [[maybe_unused]] ActiveFloat dt) noexcept {
  IGOR_ASSERT((interface.end - interface.begin).norm() >= EPS<PassiveFloat>,
              "Empty interface with begin={} and end={}",
              interface.begin,
              interface.end);

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

      return AxisAlignedWave<ActiveFloat, PointType, X>{
          .first_order_update = sign(wave_speed) * wave,
#ifdef ZAP_2ND_ORDER_CORRECTION
          .second_order_update =
              0.5 * std::abs(wave_speed) * (1 - dt / dx * std::abs(wave_speed)) * wave,
#else
          .second_order_update = 0,
#endif  //  ZAP_2ND_ORDER_CORRECTION
          .speed          = wave_speed,
          .is_right_going = wave_speed >= 0,
          .begin          = interface.begin,
          .end            = interface.end,
      };
    } else {
      IGOR_ASSERT(approx_eq(interface.begin.y, interface.end.y),
                  "Interface is not aligned to y-axis, begin={}, end={}",
                  interface.begin,
                  interface.end);

      return AxisAlignedWave<ActiveFloat, PointType, Y>{
          .first_order_update = sign(wave_speed) * wave,
#ifdef ZAP_2ND_ORDER_CORRECTION
          .second_order_update =
              0.5 * std::abs(wave_speed) * (1 - dt / dy * std::abs(wave_speed)) * wave,
#else
          .second_order_update = 0,
#endif  //  ZAP_2ND_ORDER_CORRECTION
          .speed          = wave_speed,
          .is_right_going = wave_speed >= 0,
          .begin          = interface.begin,
          .end            = interface.end,
      };
    }
  }
  // - Free wave -----------------------------------------------------------------------------------
  else {
    const PointType tangent_vector = (interface.end - interface.begin).normalized();
    const PointType normal_vector  = {tangent_vector.y, -tangent_vector.x};

    // Rotate PDE if interface does not align with an axis:
    //  tangent.y                                     = cos(interface_angle)
    //  -sign(tangent.x) * std::sqrt(1 - tangent.y^2) = sin(interface_angle)
    const ActiveFloat u_mid = (interface.left_value + interface.right_value) / 2;
    const ActiveFloat wave_speed =
        u_mid * (tangent_vector.y +
                 -sign(tangent_vector.x) * std::sqrt(1 - tangent_vector.y * tangent_vector.y));

    return FreeWave<ActiveFloat, PointType>{
        .first_order_update  = sign(wave_speed) * wave,
        .second_order_update = 0,  // TODO: How does the second_order_update work here?
        .speed               = wave_speed,
        .is_right_going      = wave_speed >= 0,
        .begin               = interface.begin,
        .end                 = interface.end,
        .normal              = scale_if_grid_coord(normal_vector, dx, dy),
    };
  }
}

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, typename PointType>
[[nodiscard]] constexpr auto calc_wave_polygon(const FreeWave<ActiveFloat, PointType>& wave,
                                               const ActiveFloat& dt) noexcept
    -> Geometry::Polygon<PointType> {
  return Geometry::Polygon<PointType>({
      wave.begin,
      wave.begin + dt * wave.speed * wave.normal,
      wave.end,
      wave.end + dt * wave.speed * wave.normal,
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

  return Geometry::Polygon<PointType>({
      wave.begin,
      wave.begin + dt * wave.speed * normal,
      wave.end,
      wave.end + dt * wave.speed * normal,
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
                          "  .speed = {}\n"
                          "  .is_right_going = {}\n"
                          "  .begin = {}\n"
                          "  .end = {}\n"
                          "  .orientation = {}\n"
                          "}}",
                          wave.first_order_update,
                          wave.second_order_update,
                          wave.speed,
                          wave.is_right_going,
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
                          "  .speed = {}\n"
                          "  .is_right_going = {}\n"
                          "  .begin = {}\n"
                          "  .end = {}\n"
                          "  .normal = {}\n"
                          "}}",
                          wave.first_order_update,
                          wave.second_order_update,
                          wave.speed,
                          wave.is_right_going,
                          wave.begin,
                          wave.end,
                          wave.normal);
  }
};

#endif  // ZAP_CELL_BASED_WAVE_HPP
