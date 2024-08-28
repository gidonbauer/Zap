#ifndef ZAP_CELL_BASED_INTERFACE_HPP_
#define ZAP_CELL_BASED_INTERFACE_HPP_

#include "CellBased/Cell.hpp"
#include "CellBased/Grid.hpp"

namespace Zap::CellBased {

template <typename Float, size_t DIM>
struct HalfInterface {
  Eigen::Vector<Float, DIM> value;
  Eigen::Vector<Float, 2> begin;
  Eigen::Vector<Float, 2> end;
};

template <typename Float, size_t DIM>
struct FullInterface {
  Eigen::Vector<Float, DIM> left_value;
  Eigen::Vector<Float, DIM> right_value;
  Eigen::Vector<Float, 2> begin;
  Eigen::Vector<Float, 2> end;
};

template <typename Float, size_t DIM>
struct CellHalfInterfaces {
  SmallVector<HalfInterface<Float, DIM>> left;
  SmallVector<HalfInterface<Float, DIM>> right;
  SmallVector<HalfInterface<Float, DIM>> bottom;
  SmallVector<HalfInterface<Float, DIM>> top;
};

template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto same_interface(const HalfInterface<Float, DIM>& i1,
                                            const HalfInterface<Float, DIM>& i2) noexcept -> bool {
  constexpr Float tol = 1e-6;
  return (i1.begin - i2.begin).norm() <= tol && (i1.end - i2.end).norm() <= tol;
}

template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto
get_outer_cell_interfaces(const Cell<Float, DIM>& cell) noexcept -> CellHalfInterfaces<Float, DIM> {
  if (cell.is_cartesian()) {
    return CellHalfInterfaces<Float, DIM>{
        .left =
            {
                {
                    .value = cell.get_cartesian().value,
                    .begin = {cell.x_min, cell.y_min},
                    .end   = {cell.x_min, cell.y_min + cell.dy},
                },
            },
        .right =
            {
                {
                    .value = cell.get_cartesian().value,
                    .begin = {cell.x_min + cell.dx, cell.y_min},
                    .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                },
            },
        .bottom =
            {
                {
                    .value = cell.get_cartesian().value,
                    .begin = {cell.x_min, cell.y_min},
                    .end   = {cell.x_min + cell.dx, cell.y_min},
                },
            },
        .top =
            {
                {
                    .value = cell.get_cartesian().value,
                    .begin = {cell.x_min, cell.y_min + cell.dy},
                    .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                },
            },
    };
  } else if (cell.is_cut()) {
    switch (cell.get_cut().type) {
      case CutType::BOTTOM_LEFT:
        return CellHalfInterfaces<Float, DIM>{
            .left =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min},
                        .end   = {cell.x_min, cell.get_cut().y2_cut},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.x_min, cell.get_cut().y2_cut},
                        .end   = {cell.x_min, cell.y_min + cell.dy},
                    },
                },
            .right =
                {
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.x_min + cell.dx, cell.y_min},
                        .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                    },
                },
            .bottom =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min},
                        .end   = {cell.get_cut().x1_cut, cell.y_min},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.get_cut().x1_cut, cell.y_min},
                        .end   = {cell.x_min + cell.dx, cell.y_min},
                    },
                },
            .top =
                {
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.x_min, cell.y_min + cell.dy},
                        .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                    },
                },
        };

      case CutType::BOTTOM_RIGHT:
        return CellHalfInterfaces<Float, DIM>{
            .left =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min},
                        .end   = {cell.x_min, cell.y_min + cell.dy},
                    },
                },
            .right =
                {
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.x_min + cell.dx, cell.y_min},
                        .end   = {cell.x_min + cell.dx, cell.get_cut().y2_cut},
                    },
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min + cell.dx, cell.get_cut().y2_cut},
                        .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                    },
                },
            .bottom =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min},
                        .end   = {cell.get_cut().x1_cut, cell.y_min},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.get_cut().x1_cut, cell.y_min},
                        .end   = {cell.x_min + cell.dx, cell.y_min},
                    },
                },
            .top =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min + cell.dy},
                        .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                    },
                },
        };

      case CutType::TOP_RIGHT:
        return CellHalfInterfaces<Float, DIM>{
            .left =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min},
                        .end   = {cell.x_min, cell.y_min + cell.dy},
                    },
                },
            .right =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min + cell.dx, cell.y_min},
                        .end   = {cell.x_min + cell.dx, cell.get_cut().y1_cut},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.x_min + cell.dx, cell.get_cut().y1_cut},
                        .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                    },
                },
            .bottom =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min},
                        .end   = {cell.x_min + cell.dx, cell.y_min},
                    },
                },
            .top =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min + cell.dy},
                        .end   = {cell.get_cut().x2_cut, cell.y_min + cell.dy},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.get_cut().x2_cut, cell.y_min + cell.dy},
                        .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                    },
                },
        };

      case CutType::TOP_LEFT:
        return CellHalfInterfaces<Float, DIM>{
            .left =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min},
                        .end   = {cell.x_min, cell.get_cut().y2_cut},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.x_min, cell.get_cut().y2_cut},
                        .end   = {cell.x_min, cell.y_min + cell.dy},
                    },
                },
            .right =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min + cell.dx, cell.y_min},
                        .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                    },
                },
            .bottom =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min},
                        .end   = {cell.x_min + cell.dx, cell.y_min},
                    },
                },
            .top =
                {
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.x_min, cell.y_min + cell.dy},
                        .end   = {cell.get_cut().x1_cut, cell.y_min + cell.dy},
                    },
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.get_cut().x1_cut, cell.y_min + cell.dy},
                        .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                    },
                },
        };

      case CutType::MIDDLE_HORI:
        return CellHalfInterfaces<Float, DIM>{
            .left =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min},
                        .end   = {cell.x_min, cell.get_cut().y2_cut},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.x_min, cell.get_cut().y2_cut},
                        .end   = {cell.x_min, cell.y_min + cell.dy},
                    },
                },
            .right =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min + cell.dx, cell.y_min},
                        .end   = {cell.x_min + cell.dx, cell.get_cut().y1_cut},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.x_min + cell.dx, cell.get_cut().y1_cut},
                        .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                    },
                },
            .bottom =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min},
                        .end   = {cell.x_min + cell.dx, cell.y_min},
                    },
                },
            .top =
                {
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.x_min, cell.y_min + cell.dy},
                        .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                    },
                },
        };

      case CutType::MIDDLE_VERT:
        return CellHalfInterfaces<Float, DIM>{
            .left =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min},
                        .end   = {cell.x_min, cell.y_min + cell.dy},
                    },
                },
            .right =
                {
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.x_min + cell.dx, cell.y_min},
                        .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
                    },
                },
            .bottom =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min},
                        .end   = {cell.get_cut().x1_cut, cell.y_min},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.get_cut().x1_cut, cell.y_min},
                        .end   = {cell.x_min + cell.dx, cell.y_min},
                    },
                },
            .top =
                {
                    {
                        .value = cell.get_cut().left_value,
                        .begin = {cell.x_min, cell.y_min + cell.dy},
                        .end   = {cell.get_cut().x2_cut, cell.y_min + cell.dy},
                    },
                    {
                        .value = cell.get_cut().right_value,
                        .begin = {cell.get_cut().x2_cut, cell.y_min + cell.dy},
                        .end   = {cell.x_min + cell.dx, cell.y_min + cell.dy},
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

template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto
get_shared_interfaces(const Cell<Float, DIM>& center_cell,
                      const Cell<Float, DIM>& other_cell,
                      Side side) noexcept -> SmallVector<FullInterface<Float, DIM>> {
  const auto center_interfaces = get_outer_cell_interfaces(center_cell);
  const auto other_interfaces  = get_outer_cell_interfaces(other_cell);

  SmallVector<HalfInterface<Float, DIM>> left_side{};
  SmallVector<HalfInterface<Float, DIM>> right_side{};
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
    SmallVector<FullInterface<Float, DIM>> interfaces(n);
    for (size_t i = 0; i < n; ++i) {
      assert(same_interface(left_side[i], right_side[i]));
      interfaces[i] = FullInterface<Float, DIM>{
          .left_value  = left_side[i].value,
          .right_value = right_side[i].value,
          .begin       = left_is_center ? left_side[i].begin : right_side[i].begin,
          .end         = left_is_center ? left_side[i].end : right_side[i].end,
      };
    }
    return interfaces;
  } else {
    // if (left_side.size() == 1) {
    //   Igor::Todo("left_side.size() == 1");
    // } else if (right_side.size() == 1) {
    //   Igor::Todo("right_side.size() == 1");
    // } else {
    //   Igor::Panic("Only simple case when either left_side or right side have only one interface
    //   is "
    //               "implemented.");
    //   std::unreachable();
    // }

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
