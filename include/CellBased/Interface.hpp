#ifndef ZAP_CELL_BASED_INTERFACE_HPP_
#define ZAP_CELL_BASED_INTERFACE_HPP_

#include "CellBased/Cell.hpp"
#include "CellBased/Definitions.hpp"
#include "CellBased/Grid.hpp"

namespace Zap::CellBased {

// TODO: Fix interfaces

template <typename Float, size_t DIM, Point2D_c PointType>
struct HalfInterface {
  Eigen::Vector<Float, DIM> value;
  PointType begin;
  PointType end;
};

template <typename Float, size_t DIM, Point2D_c PointType>
struct FullInterface {
  Eigen::Vector<Float, DIM> left_value;
  Eigen::Vector<Float, DIM> right_value;
  PointType begin;
  PointType end;
};

template <typename Float, size_t DIM, Point2D_c PointType>
struct CellHalfInterfaces {
  SmallVector<HalfInterface<Float, DIM, PointType>> left;
  SmallVector<HalfInterface<Float, DIM, PointType>> right;
  SmallVector<HalfInterface<Float, DIM, PointType>> bottom;
  SmallVector<HalfInterface<Float, DIM, PointType>> top;
};

template <typename Float, size_t DIM, Point2D_c PointType>
[[nodiscard]] constexpr auto same_interface(const HalfInterface<Float, DIM, PointType>& i1,
                                            const HalfInterface<Float, DIM, PointType>& i2) noexcept
    -> bool {
  constexpr Float tol = 1e-6;
  return (i1.begin - i2.begin).norm() <= tol && (i1.end - i2.end).norm() <= tol;
}

template <typename Float, size_t DIM, Point2D_c PointType>
[[nodiscard]] constexpr auto get_outer_cell_interfaces(const Cell<Float, DIM>& cell) noexcept
    -> CellHalfInterfaces<Float, DIM, PointType> {
  constexpr CoordType coord_type = std::is_same_v<std::remove_cvref_t<PointType>, SimCoord<Float>>
                                       ? CoordType::SIM
                                       : CoordType::GRID;

  if (cell.is_cartesian()) {
    return CellHalfInterfaces<Float, DIM, PointType>{
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
                    .begin = {cell.template x_min<coord_type>(), cell.template y_min<coord_type>()},
                    .end   = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                              cell.template y_min<coord_type>()},
                },
            },
        .top =
            {
                {
                    .value = cell.get_cartesian().value,
                    .begin = {cell.template x_min<coord_type>(),
                              cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    .end   = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                              cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                },
            },
    };
  } else if (cell.is_cut()) {
    switch (cell.get_cut().type) {
      case CutType::BOTTOM_LEFT:
        return CellHalfInterfaces<Float, DIM, PointType>{
            .left =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end   = {cell.template x_min<coord_type>(),
                                  cell.template cut2<coord_type>().y},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template cut2<coord_type>().y},
                        .end   = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
            .right =
                {
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template x_min<coord_type>() +
                                      cell.template dx<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
            .bottom =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end   = {cell.template cut1<coord_type>().x,
                                  cell.template y_min<coord_type>()},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template cut1<coord_type>().x,
                                  cell.template y_min<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>()},
                    },
                },
            .top =
                {
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>() +
                                      cell.template dy<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
        };

      case CutType::BOTTOM_RIGHT:
        return CellHalfInterfaces<Float, DIM, PointType>{
            .left =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end   = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
            .right =
                {
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template x_min<coord_type>() +
                                      cell.template dx<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template cut2<coord_type>().y},
                    },
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>() +
                                      cell.template dx<coord_type>(),
                                  cell.template cut2<coord_type>().y},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
            .bottom =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end   = {cell.template cut1<coord_type>().x,
                                  cell.template y_min<coord_type>()},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template cut1<coord_type>().x,
                                  cell.template y_min<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>()},
                    },
                },
            .top =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>() +
                                      cell.template dy<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
        };

      case CutType::TOP_RIGHT:
        return CellHalfInterfaces<Float, DIM, PointType>{
            .left =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end   = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
            .right =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>() +
                                      cell.template dx<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template cut1<coord_type>().y},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template x_min<coord_type>() +
                                      cell.template dx<coord_type>(),
                                  cell.template cut1<coord_type>().y},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
            .bottom =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>()},
                    },
                },
            .top =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>() +
                                      cell.template dy<coord_type>()},
                        .end   = {cell.template cut2<coord_type>().x,
                                  cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template cut2<coord_type>().x,
                                  cell.template y_min<coord_type>() +
                                      cell.template dy<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
        };

      case CutType::TOP_LEFT:
        return CellHalfInterfaces<Float, DIM, PointType>{
            .left =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end   = {cell.template x_min<coord_type>(),
                                  cell.template cut2<coord_type>().y},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template cut2<coord_type>().y},
                        .end   = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
            .right =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>() +
                                      cell.template dx<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
            .bottom =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>()},
                    },
                },
            .top =
                {
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>() +
                                      cell.template dy<coord_type>()},
                        .end   = {cell.template cut1<coord_type>().x,
                                  cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template cut1<coord_type>().x,
                                  cell.template y_min<coord_type>() +
                                      cell.template dy<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
        };

      case CutType::MIDDLE_HORI:
        return CellHalfInterfaces<Float, DIM, PointType>{
            .left =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end   = {cell.template x_min<coord_type>(),
                                  cell.template cut2<coord_type>().y},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template cut2<coord_type>().y},
                        .end   = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
            .right =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>() +
                                      cell.template dx<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template cut1<coord_type>().y},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template x_min<coord_type>() +
                                      cell.template dx<coord_type>(),
                                  cell.template cut1<coord_type>().y},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
            .bottom =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>()},
                    },
                },
            .top =
                {
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>() +
                                      cell.template dy<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
        };

      case CutType::MIDDLE_VERT:
        return CellHalfInterfaces<Float, DIM, PointType>{
            .left =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end   = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
            .right =
                {
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template x_min<coord_type>() +
                                      cell.template dx<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
            .bottom =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>()},
                        .end   = {cell.template cut1<coord_type>().x,
                                  cell.template y_min<coord_type>()},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template cut1<coord_type>().x,
                                  cell.template y_min<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>()},
                    },
                },
            .top =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.template x_min<coord_type>(),
                                  cell.template y_min<coord_type>() +
                                      cell.template dy<coord_type>()},
                        .end   = {cell.template cut2<coord_type>().x,
                                  cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.template cut2<coord_type>().x,
                                  cell.template y_min<coord_type>() +
                                      cell.template dy<coord_type>()},
                        .end = {cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                                cell.template y_min<coord_type>() + cell.template dy<coord_type>()},
                    },
                },
        };

      default:
        Igor::Panic("Unknown cut type with value {}", std::to_underlying(cell.get_cut().type));
        std::unreachable();
    }
    Igor::Panic("Unreachable");
    std::unreachable();
  } else {
    Igor::Panic("Unknown cell type with variant index {}", cell.value.index());
    std::unreachable();
  }
}

template <typename Float, size_t DIM, Point2D_c PointType>
[[nodiscard]] constexpr auto get_shared_interfaces(const Cell<Float, DIM>& center_cell,
                                                   const Cell<Float, DIM>& other_cell,
                                                   Side side) noexcept
    -> SmallVector<FullInterface<Float, DIM, PointType>> {
  const auto center_interfaces = get_outer_cell_interfaces<Float, DIM, PointType>(center_cell);
  const auto other_interfaces  = get_outer_cell_interfaces<Float, DIM, PointType>(other_cell);

  SmallVector<HalfInterface<Float, DIM, PointType>> left_side{};
  SmallVector<HalfInterface<Float, DIM, PointType>> right_side{};
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
      Igor::Panic("side must be LEFT, RIGHT, BOTTOM, or TOP, but has value {}",
                  std::to_underlying(side));
      std::unreachable();
  }

  if (left_side.size() == right_side.size()) {
    const auto n = left_side.size();
    SmallVector<FullInterface<Float, DIM, PointType>> interfaces(n);
    for (size_t i = 0; i < n; ++i) {
      // TODO: Is this necessary?
      // assert(same_interface(left_side[i], right_side[i]));
      interfaces[i] = FullInterface<Float, DIM, PointType>{
          .left_value  = left_side[i].value,
          .right_value = right_side[i].value,
          // TODO: Is this necessary?
          .begin = left_is_center ? left_side[i].begin : right_side[i].begin,
          .end   = left_is_center ? left_side[i].end : right_side[i].end,
      };
    }
    return interfaces;
  } else {
    if (left_side.size() == 1) {
      SmallVector<FullInterface<Float, DIM, PointType>> interfaces(right_side.size());
      for (size_t i = 0; i < right_side.size(); ++i) {
        interfaces[i] = FullInterface<Float, DIM, PointType>{
            .left_value  = left_side[0].value,
            .right_value = right_side[i].value,
            .begin       = right_side[i].begin,
            .end         = right_side[i].end,
        };
      }
      return interfaces;
    } else if (right_side.size() == 1) {
      SmallVector<FullInterface<Float, DIM, PointType>> interfaces(left_side.size());
      for (size_t i = 0; i < left_side.size(); ++i) {
        interfaces[i] = FullInterface<Float, DIM, PointType>{
            .left_value  = left_side[i].value,
            .right_value = right_side[0].value,
            .begin       = left_side[i].begin,
            .end         = left_side[i].end,
        };
      }
      return interfaces;
    } else {
      Igor::Panic("Only simple case when either left_side or right side have only one interface is "
                  "implemented.");
      std::unreachable();
    }

    Igor::Debug("left_side:");
    for (size_t i = 0; i < left_side.size(); ++i) {
      Igor::Debug("[{}] = {{ .value = {}, .begin = {}, .end = {} }}",
                  i,
                  left_side[i].value,
                  left_side[i].begin,
                  left_side[i].end);
    }

    Igor::Debug("right_side:");
    for (size_t i = 0; i < right_side.size(); ++i) {
      Igor::Debug("[{}] = {{ .value = {}, .begin = {}, .end = {} }}",
                  i,
                  right_side[i].value,
                  right_side[i].begin,
                  right_side[i].end);
    }

    {
      std::stringstream s;
      s << center_cell;
      Igor::Debug("center_cell = {}", s.str());
    }
    {
      std::stringstream s;
      s << other_cell;
      Igor::Debug("other_cell = {}", s.str());
    }

    Igor::Todo("Handle incompatible interfaces.");
    std::unreachable();
  }
}

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_INTERFACE_HPP_
