#ifndef ZAP_CELL_BASED_INTERFACE_HPP_
#define ZAP_CELL_BASED_INTERFACE_HPP_

#include "CellBased/Cell.hpp"
#include "CellBased/Definitions.hpp"
#include "CellBased/Grid.hpp"

namespace Zap::CellBased {

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, Point2D_c PointType>
struct HalfInterface {
  ActiveFloat value;
  PointType begin;
  PointType end;
};

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, Point2D_c PointType>
struct FullInterface {
  ActiveFloat left_value;
  ActiveFloat right_value;
  PointType begin;
  PointType end;
};

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, Point2D_c PointType>
struct CellHalfInterfaces {
  SmallVector<HalfInterface<ActiveFloat, PointType>> left;
  SmallVector<HalfInterface<ActiveFloat, PointType>> right;
  SmallVector<HalfInterface<ActiveFloat, PointType>> bottom;
  SmallVector<HalfInterface<ActiveFloat, PointType>> top;
};

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, Point2D_c PointType>
[[nodiscard]] constexpr auto
same_interface(const HalfInterface<ActiveFloat, PointType>& i1,
               const HalfInterface<ActiveFloat, PointType>& i2) noexcept -> bool {
  return (i1.begin - i2.begin).norm() <= EPS<ActiveFloat> &&
         (i1.end - i2.end).norm() <= EPS<ActiveFloat>;
}

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, typename PassiveFloat, Point2D_c PointType>
[[nodiscard]] constexpr auto
get_outer_cell_interfaces(const Cell<ActiveFloat, PassiveFloat>& cell) noexcept
    -> CellHalfInterfaces<ActiveFloat, PointType> {
  constexpr CoordType coord_type = PointType2CoordType<PointType>;

  if (cell.is_cartesian()) {
    return CellHalfInterfaces<ActiveFloat, PointType>{
        .left =
            {
                {
                    .value = cell.get_cartesian().value,
                    .begin = {cell.template x_min<coord_type>(), cell.template y_min<coord_type>()},
                    .end   = {cell.template x_min<coord_type>(),
                              cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                },
            },
        .right =
            {
                {
                    .value = cell.get_cartesian().value,
                    .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                              cell.template y_min<coord_type>()},
                    .end   = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                              cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                },
            },
        .bottom =
            {
                {
                    .value = cell.get_cartesian().value,
                    .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                              cell.template y_min<coord_type>()},
                    .end   = {cell.template x_min<coord_type>(), cell.template y_min<coord_type>()},
                },
            },
        .top =
            {
                {
                    .value = cell.get_cartesian().value,
                    .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                              cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    .end   = {cell.template x_min<coord_type>(),
                              cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                },
            },
    };
  } else if (cell.is_cut()) {
    CellHalfInterfaces<ActiveFloat, PointType> outer_interfaces;

    // - Left-side ---------------------------------------------------------------------------------
    if ((cell.get_cut().entry_loc & LEFT) > 0) {
      outer_interfaces.left = {
          {
              .value = cell.get_cut().right_value,
              .begin = {cell.template x_min<coord_type>(), cell.template y_min<coord_type>()},
              .end   = {cell.template x_min<coord_type>(), cell.template cut_entry<coord_type>().y},
          },
          {
              .value = cell.get_cut().left_value,
              .begin = {cell.template x_min<coord_type>(), cell.template cut_entry<coord_type>().y},
              .end   = {cell.template x_min<coord_type>(),
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
          },
      };
    } else if ((cell.get_cut().exit_loc & LEFT) > 0) {
      outer_interfaces.left = {
          {
              .value = cell.get_cut().left_value,
              .begin = {cell.template x_min<coord_type>(), cell.template y_min<coord_type>()},
              .end   = {cell.template x_min<coord_type>(), cell.template cut_exit<coord_type>().y},
          },
          {
              .value = cell.get_cut().right_value,
              .begin = {cell.template x_min<coord_type>(), cell.template cut_exit<coord_type>().y},
              .end   = {cell.template x_min<coord_type>(),
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
          },
      };
    } else {
      outer_interfaces.left = {
          {
              .value = cell.get_cut().entry_loc < cell.get_cut().exit_loc
                           ? cell.get_cut().left_value
                           : cell.get_cut().right_value,
              .begin = {cell.template x_min<coord_type>(), cell.template y_min<coord_type>()},
              .end   = {cell.template x_min<coord_type>(),
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
          },
      };
    }

    // - Right-side --------------------------------------------------------------------------------
    if ((cell.get_cut().entry_loc & RIGHT) > 0) {
      outer_interfaces.right = {
          {
              .value = cell.get_cut().left_value,
              .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template y_min<coord_type>()},
              .end   = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template cut_entry<coord_type>().y},
          },
          {
              .value = cell.get_cut().right_value,
              .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template cut_entry<coord_type>().y},
              .end   = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
          },
      };
    } else if ((cell.get_cut().exit_loc & RIGHT) > 0) {
      outer_interfaces.right = {
          {
              .value = cell.get_cut().right_value,
              .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template y_min<coord_type>()},
              .end   = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template cut_exit<coord_type>().y},
          },
          {
              .value = cell.get_cut().left_value,
              .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template cut_exit<coord_type>().y},
              .end   = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
          },
      };
    } else {
      ActiveFloat value;
      const auto cut_type = make_cut_type(cell.get_cut().entry_loc, cell.get_cut().exit_loc);
      switch (cut_type) {
        case make_cut_type(BOTTOM, TOP):
        case make_cut_type(LEFT, TOP):
        case make_cut_type(BOTTOM, LEFT): value = cell.get_cut().right_value; break;
        case make_cut_type(TOP, BOTTOM):
        case make_cut_type(TOP, LEFT):
        case make_cut_type(LEFT, BOTTOM): value = cell.get_cut().left_value; break;
        default:
          Igor::Panic("Invalid combination of entry_loc = {} and exit_loc = {}",
                      cell.get_cut().entry_loc,
                      cell.get_cut().exit_loc);
          std::unreachable();
      }

      outer_interfaces.right = {
          {
              .value = value,
              .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template y_min<coord_type>()},
              .end   = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
          },
      };
    }

    // - Bottom-side -------------------------------------------------------------------------------
    if ((cell.get_cut().entry_loc & BOTTOM) > 0) {
      outer_interfaces.bottom = {
          {
              .value = cell.get_cut().left_value,
              .begin = {cell.template cut_entry<coord_type>().x, cell.template y_min<coord_type>()},
              .end   = {cell.template x_min<coord_type>(), cell.template y_min<coord_type>()},
          },
          {
              .value = cell.get_cut().right_value,
              .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template y_min<coord_type>()},
              .end   = {cell.template cut_entry<coord_type>().x, cell.template y_min<coord_type>()},
          },
      };
    } else if ((cell.get_cut().exit_loc & BOTTOM) > 0) {
      outer_interfaces.bottom = {
          {
              .value = cell.get_cut().right_value,
              .begin = {cell.template cut_exit<coord_type>().x, cell.template y_min<coord_type>()},
              .end   = {cell.template x_min<coord_type>(), cell.template y_min<coord_type>()},
          },
          {
              .value = cell.get_cut().left_value,
              .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template y_min<coord_type>()},
              .end   = {cell.template cut_exit<coord_type>().x, cell.template y_min<coord_type>()},
          },
      };
    } else {
      outer_interfaces.bottom = {{
          .value = cell.get_cut().entry_loc < cell.get_cut().exit_loc ? cell.get_cut().left_value
                                                                      : cell.get_cut().right_value,
          .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                    cell.template y_min<coord_type>()},
          .end   = {cell.template x_min<coord_type>(), cell.template y_min<coord_type>()},
      }};
    }

    // - Top-side ----------------------------------------------------------------------------------
    if ((cell.get_cut().entry_loc & TOP) > 0) {
      outer_interfaces.top = {
          {
              .value = cell.get_cut().right_value,
              .begin = {cell.template cut_entry<coord_type>().x,
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
              .end   = {cell.template x_min<coord_type>(),
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
          },
          {
              .value = cell.get_cut().left_value,
              .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
              .end   = {cell.template cut_entry<coord_type>().x,
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
          },
      };
    } else if ((cell.get_cut().exit_loc & TOP) > 0) {
      outer_interfaces.top = {
          {
              .value = cell.get_cut().left_value,
              .begin = {cell.template cut_exit<coord_type>().x,
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
              .end   = {cell.template x_min<coord_type>(),
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
          },
          {
              .value = cell.get_cut().right_value,
              .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
              .end   = {cell.template cut_exit<coord_type>().x,
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
          },
      };
    } else {
      ActiveFloat value;
      const auto cut_type = make_cut_type(cell.get_cut().entry_loc, cell.get_cut().exit_loc);
      switch (cut_type) {
        case make_cut_type(RIGHT, LEFT):
        case make_cut_type(BOTTOM, LEFT):
        case make_cut_type(RIGHT, BOTTOM): value = cell.get_cut().right_value; break;
        case make_cut_type(LEFT, RIGHT):
        case make_cut_type(LEFT, BOTTOM):
        case make_cut_type(BOTTOM, RIGHT): value = cell.get_cut().left_value; break;
        default:
          Igor::Panic("Invalid combination of entry_loc = {} and exit_loc = {}",
                      cell.get_cut().entry_loc,
                      cell.get_cut().exit_loc);
          std::unreachable();
      }

      outer_interfaces.top = {
          {
              .value = value,
              .begin = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
              .end   = {cell.template x_min<coord_type>(),
                        cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
          },
      };
    }

    return outer_interfaces;
  } else {
    Igor::Panic("Unknown cell type with variant index {}", cell.cell_type.index());
    std::unreachable();
  }
}

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, typename PassiveFloat, Point2D_c PointType>
[[nodiscard]] constexpr auto
get_shared_interfaces(const Cell<ActiveFloat, PassiveFloat>& center_cell,
                      const Cell<ActiveFloat, PassiveFloat>& other_cell,
                      Side side) noexcept -> SmallVector<FullInterface<ActiveFloat, PointType>> {
  const auto center_interfaces =
      get_outer_cell_interfaces<ActiveFloat, PassiveFloat, PointType>(center_cell);
  const auto other_interfaces =
      get_outer_cell_interfaces<ActiveFloat, PassiveFloat, PointType>(other_cell);

  SmallVector<HalfInterface<ActiveFloat, PointType>> left_side{};
  SmallVector<HalfInterface<ActiveFloat, PointType>> right_side{};
  bool left_is_center;
  switch (side) {
    case LEFT:
      {
        left_side      = other_interfaces.right;
        right_side     = center_interfaces.left;
        left_is_center = false;
      }
      break;
    case RIGHT:
      {
        left_side      = center_interfaces.right;
        right_side     = other_interfaces.left;
        left_is_center = true;
      }
      break;
    case BOTTOM:
      {
        left_side      = other_interfaces.top;
        right_side     = center_interfaces.bottom;
        left_is_center = false;
      }
      break;
    case TOP:
      {
        left_side      = center_interfaces.top;
        right_side     = other_interfaces.bottom;
        left_is_center = true;
      }
      break;
    default:
      Igor::Panic("side must be LEFT, RIGHT, BOTTOM, or TOP, but has value {}", side);
      std::unreachable();
  }

  // = Same amount of interfaces on left and right side ============================================
  if (left_side.size() == right_side.size()) {
    const auto n = left_side.size();
    SmallVector<FullInterface<ActiveFloat, PointType>> interfaces(n);
    for (size_t i = 0; i < n; ++i) {
      // TODO: Do this check for periodic boundary conditions.
      // IGOR_ASSERT(same_interface(left_side[i], right_side[i]),
      //             "Incompatible interface: left = {} -> {} and right = {} -> {}",
      //             left_side[i].begin,
      //             left_side[i].end,
      //             right_side[i].begin,
      //             right_side[i].end);
      interfaces[i] = FullInterface<ActiveFloat, PointType>{
          .left_value  = left_side[i].value,
          .right_value = right_side[i].value,
          .begin       = left_is_center ? left_side[i].begin : right_side[i].begin,
          .end         = left_is_center ? left_side[i].end : right_side[i].end,
      };
    }
    return interfaces;
  }
  // = Different amount of interfaces on left and right side =======================================
  else {
    // = Right side is resonsible for dividing the interface =======================================
    if (left_side.size() == 1) {
      SmallVector<FullInterface<ActiveFloat, PointType>> interfaces(right_side.size());
      for (size_t i = 0; i < right_side.size(); ++i) {
        PointType begin = right_side[i].begin;
        PointType end   = right_side[i].end;
        if (left_is_center) {
          if (side == LEFT || side == RIGHT) {
            begin.x = left_side[0].begin.x;
            end.x   = left_side[0].end.x;
          } else if (side == BOTTOM || side == TOP) {
            begin.y = left_side[0].begin.y;
            end.y   = left_side[0].end.y;
          }
        }

        interfaces[i] = FullInterface<ActiveFloat, PointType>{
            .left_value  = left_side[0].value,
            .right_value = right_side[i].value,
            .begin       = begin,
            .end         = end,
        };
      }
      return interfaces;
    }
    // = Left side is resonsible for dividing the interface ========================================
    else if (right_side.size() == 1) {
      SmallVector<FullInterface<ActiveFloat, PointType>> interfaces(left_side.size());
      for (size_t i = 0; i < left_side.size(); ++i) {
        PointType begin = left_side[i].begin;
        PointType end   = left_side[i].end;
        if (!left_is_center) {
          if (side == LEFT || side == RIGHT) {
            begin.x = right_side[0].begin.x;
            end.x   = right_side[0].end.x;
          } else if (side == BOTTOM || side == TOP) {
            begin.y = right_side[0].begin.y;
            end.y   = right_side[0].end.y;
          }
        }

        interfaces[i] = FullInterface<ActiveFloat, PointType>{
            .left_value  = left_side[i].value,
            .right_value = right_side[0].value,
            .begin       = begin,
            .end         = end,
        };
      }
      return interfaces;
    }
    // = Multiple incompatible interfaces ==========================================================
    Igor::Panic("Only simple case when either left_side or right side have only one interface is "
                "implemented.");
    std::unreachable();
  }
}

}  // namespace Zap::CellBased

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, Zap::CellBased::Point2D_c PointType>
struct fmt::formatter<Zap::CellBased::HalfInterface<ActiveFloat, PointType>> {
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) const noexcept {
    return ctx.begin();
  }
  template <typename FormatContext>
  constexpr auto format(const Zap::CellBased::HalfInterface<ActiveFloat, PointType>& half_interface,
                        FormatContext& ctx) const noexcept {
    return fmt::format_to(ctx.out(),
                          "{{\n"
                          "  .value = {}\n"
                          "  .begin = {}\n"
                          "  .end   = {}\n"
                          "}}",
                          half_interface.value,
                          half_interface.begin,
                          half_interface.end);
  }
};

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, Zap::CellBased::Point2D_c PointType>
struct fmt::formatter<Zap::CellBased::FullInterface<ActiveFloat, PointType>> {
  template <typename ParseContext>
  constexpr auto parse(ParseContext& ctx) const noexcept {
    return ctx.begin();
  }
  template <typename FormatContext>
  constexpr auto format(const Zap::CellBased::FullInterface<ActiveFloat, PointType>& full_interface,
                        FormatContext& ctx) const noexcept {
    return fmt::format_to(ctx.out(),
                          "{{\n"
                          "  .left_value  = {}\n"
                          "  .right_value = {}\n"
                          "  .begin       = {}\n"
                          "  .end         = {}\n"
                          "}}",
                          full_interface.left_value,
                          full_interface.right_value,
                          full_interface.begin,
                          full_interface.end);
  }
};

#endif  // ZAP_CELL_BASED_INTERFACE_HPP_
