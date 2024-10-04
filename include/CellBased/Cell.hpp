#ifndef ZAP_CELL_BASED_CELL_HPP_
#define ZAP_CELL_BASED_CELL_HPP_

#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <variant>

#include <Eigen/Dense>

#include "CellBased/Definitions.hpp"
#include "CellBased/Geometry.hpp"
#include "CellBased/SmallVector.hpp"
#include "Igor.hpp"

namespace Zap::CellBased {

// -------------------------------------------------------------------------------------------------
enum : size_t {
  SAME_VALUE_INDEX = std::numeric_limits<size_t>::max() - 2,
  ZERO_FLUX_INDEX  = std::numeric_limits<size_t>::max() - 1,
  NULL_INDEX       = std::numeric_limits<size_t>::max() - 0,
};

// -------------------------------------------------------------------------------------------------
enum class CutType : char {
  BOTTOM_LEFT,
  BOTTOM_RIGHT,
  TOP_RIGHT,
  TOP_LEFT,
  MIDDLE_HORI,
  MIDDLE_VERT,
};

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
struct CartesianValue {
  static_assert(DIM > 0, "Dimension must be at least 1.");
  Eigen::Vector<Float, DIM> value = Eigen::Vector<Float, DIM>::Zero();
};

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
struct CutValue {
  static_assert(DIM > 0, "Dimension must be at least 1.");
  Eigen::Vector<Float, DIM> left_value  = Eigen::Vector<Float, DIM>::Zero();
  Eigen::Vector<Float, DIM> right_value = Eigen::Vector<Float, DIM>::Zero();

  // Linear cut
  CutType type{};
  Point<Float> cut1{};  // Entry point
  Point<Float> cut2{};  // Exit point
};

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
class UniformGrid;

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
struct Cell {
  static_assert(DIM > 0, "Dimension must be at least 1.");

  // Cell value
  std::variant<CartesianValue<Float, DIM>, CutValue<Float, DIM>> value =
      CartesianValue<Float, DIM>{};

  UniformGrid<Float, DIM>* parent;

  // Grid index
  size_t x_idx{};
  size_t y_idx{};

  // Connectivity
  size_t left_idx   = NULL_INDEX;
  size_t right_idx  = NULL_INDEX;
  size_t bottom_idx = NULL_INDEX;
  size_t top_idx    = NULL_INDEX;

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto x_min() const noexcept -> Float {
    return parent->m_x_min + static_cast<Float>(x_idx) * parent->m_dx;
  }
  [[nodiscard]] constexpr auto y_min() const noexcept -> Float {
    return parent->m_y_min + static_cast<Float>(y_idx) * parent->m_dy;
  }
  [[nodiscard]] constexpr auto dx() const noexcept -> Float { return parent->m_dx; }
  [[nodiscard]] constexpr auto dy() const noexcept -> Float { return parent->m_dy; }

  [[nodiscard]] constexpr auto cut1() const noexcept -> Point<Float> {
    assert(is_cut());
    return Point<Float>{
        get_cut().cut1(X) * dx() + x_min(),
        get_cut().cut1(Y) * dy() + y_min(),
    };
  }

  [[nodiscard]] constexpr auto cut2() const noexcept -> Point<Float> {
    assert(is_cut());
    return Point<Float>{
        get_cut().cut2(X) * dx() + x_min(),
        get_cut().cut2(Y) * dy() + y_min(),
    };
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto is_cartesian() const noexcept -> bool {
    return std::holds_alternative<CartesianValue<Float, DIM>>(value);
  }
  [[nodiscard]] constexpr auto is_cut() const noexcept -> bool {
    return std::holds_alternative<CutValue<Float, DIM>>(value);
  }
  [[nodiscard]] constexpr auto get_cartesian() noexcept -> CartesianValue<Float, DIM>& {
    assert(is_cartesian());
    return std::get<CartesianValue<Float, DIM>>(value);
  }
  [[nodiscard]] constexpr auto get_cartesian() const noexcept -> const CartesianValue<Float, DIM>& {
    assert(is_cartesian());
    return std::get<CartesianValue<Float, DIM>>(value);
  }
  [[nodiscard]] constexpr auto get_cut() noexcept -> CutValue<Float, DIM>& {
    assert(is_cut());
    return std::get<CutValue<Float, DIM>>(value);
  }
  [[nodiscard]] constexpr auto get_cut() const noexcept -> const CutValue<Float, DIM>& {
    assert(is_cut());
    return std::get<CutValue<Float, DIM>>(value);
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto get_cartesian_polygon() const noexcept -> Geometry::Polygon<Float> {
    return Geometry::Polygon<Float>{{
        Point<Float>{x_min(), y_min()},
        Point<Float>{x_min() + dx(), y_min()},
        Point<Float>{x_min(), y_min() + dy()},
        Point<Float>{x_min() + dx(), y_min() + dy()},
    }};
  }
  [[nodiscard]] constexpr auto get_cut_left_polygon() const noexcept -> Geometry::Polygon<Float> {
    assert(is_cut());
    return Geometry::Polygon<Float>{get_left_points()};
  }
  [[nodiscard]] constexpr auto get_cut_right_polygon() const noexcept -> Geometry::Polygon<Float> {
    assert(is_cut());
    return Geometry::Polygon<Float>{get_right_points()};
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto get_left_points() const noexcept -> SmallVector<Point<Float>> {
    const auto& cell_value = get_cut();
    switch (cell_value.type) {
      case CutType::BOTTOM_LEFT:
        return {
            Point<Float>{x_min(), y_min()},
            cut1(),
            cut2(),
        };
      case CutType::BOTTOM_RIGHT:
        return {
            Point<Float>{x_min(), y_min()},
            cut1(),
            cut2(),
            Point<Float>{x_min() + dx(), y_min() + dy()},
            Point<Float>{x_min(), y_min() + dy()},
        };
      case CutType::TOP_RIGHT:
        return {
            Point<Float>{x_min(), y_min()},
            Point<Float>{x_min() + dx(), y_min()},
            cut1(),
            cut2(),
            Point<Float>{x_min(), y_min() + dy()},
        };
      case CutType::TOP_LEFT:
        return {
            Point<Float>{x_min(), y_min()},
            Point<Float>{x_min() + dx(), y_min()},
            Point<Float>{x_min() + dx(), y_min() + dy()},
            cut1(),
            cut2(),
        };
      case CutType::MIDDLE_HORI:
        return {
            Point<Float>{x_min(), y_min()},
            Point<Float>{x_min() + dx(), y_min()},
            cut1(),
            cut2(),
        };
      case CutType::MIDDLE_VERT:
        return {
            Point<Float>{x_min(), y_min()},
            Point<Float>{x_min(), y_min() + dy()},
            cut1(),
            cut2(),
        };
      default:
        Igor::Panic("Unknown cut type with value {}", static_cast<int>(cell_value.type));
        std::unreachable();
    }
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto get_right_points() const noexcept -> SmallVector<Point<Float>> {
    const auto& cell_value = get_cut();
    switch (cell_value.type) {
      case CutType::BOTTOM_LEFT:
        return {
            cut1(),
            cut2(),
            Point<Float>{x_min(), y_min() + dy()},
            Point<Float>{x_min() + dx(), y_min()},
            Point<Float>{x_min() + dx(), y_min() + dy()},
        };
      case CutType::BOTTOM_RIGHT:
        return {
            cut1(),
            cut2(),
            Point<Float>{x_min() + dx(), y_min()},
        };
      case CutType::TOP_RIGHT:
        return {
            cut1(),
            cut2(),
            Point<Float>{x_min() + dx(), y_min() + dy()},
        };
      case CutType::TOP_LEFT:
        return {
            cut1(),
            cut2(),
            Point<Float>{x_min(), y_min() + dy()},
        };
      case CutType::MIDDLE_HORI:
        return {
            cut1(),
            cut2(),
            Point<Float>{x_min(), y_min() + dy()},
            Point<Float>{x_min() + dx(), y_min() + dy()},
        };
      case CutType::MIDDLE_VERT:
        return {
            cut1(),
            Point<Float>{x_min() + dx(), y_min()},
            cut2(),
            Point<Float>{x_min() + dx(), y_min() + dy()},
        };
      default:
        Igor::Panic("Unknown cut type with value {}", static_cast<int>(cell_value.type));
        std::unreachable();
    }
  }
};

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
auto operator<<(std::ostream& out, const CartesianValue<Float, DIM>& cart_value) noexcept
    -> std::ostream& {
  constexpr auto end_char      = '\n';
  constexpr auto single_indent = "  ";
  constexpr auto double_indent = "    ";

  out << '{' << end_char;

  out << double_indent << ".value = [ ";
  out << cart_value.value(0);
  for (Eigen::Index i = 1; i < static_cast<Eigen::Index>(DIM); ++i) {
    out << ", " << cart_value.value(i);
  }
  out << " ]," << end_char;

  out << single_indent << '}';

  return out;
}

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
auto operator<<(std::ostream& out, const CutValue<Float, DIM>& cut_value) noexcept
    -> std::ostream& {
  constexpr auto end_char      = '\n';
  constexpr auto single_indent = "  ";
  constexpr auto double_indent = "    ";

  out << '{' << end_char;

  out << double_indent << ".left_value = [ ";
  out << cut_value.left_value(0);
  for (Eigen::Index i = 1; i < static_cast<Eigen::Index>(DIM); ++i) {
    out << ", " << cut_value.left_value(i);
  }
  out << " ]," << end_char;

  out << double_indent << ".right_value = [ ";
  out << cut_value.right_value(0);
  for (Eigen::Index i = 1; i < static_cast<Eigen::Index>(DIM); ++i) {
    out << ", " << cut_value.right_value(i);
  }
  out << " ]," << end_char;

  out << double_indent << ".type = ";
  switch (cut_value.type) {
    case CutType::BOTTOM_LEFT:  out << "BOTTOM_LEFT"; break;
    case CutType::BOTTOM_RIGHT: out << "BOTTOM_RIGHT"; break;
    case CutType::TOP_RIGHT:    out << "TOP_RIGHT"; break;
    case CutType::TOP_LEFT:     out << "TOP_LEFT"; break;
    case CutType::MIDDLE_HORI:  out << "MIDDLE_HORI"; break;
    case CutType::MIDDLE_VERT:  out << "MIDDLE_VERT"; break;
  }
  out << ',' << end_char;

  out << double_indent << ".cut1 = [" << cut_value.cut1(X) << ", " << cut_value.cut1(Y) << ']'
      << end_char;
  out << double_indent << ".cut2 = [" << cut_value.cut2(X) << ", " << cut_value.cut2(Y) << ']'
      << end_char;

  out << single_indent << '}';

  return out;
}

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
auto operator<<(std::ostream& out, const Cell<Float, DIM>& cell) noexcept -> std::ostream& {
  constexpr auto end_char = '\n';
  constexpr auto indent   = "  ";

  out << "{" << end_char;

  out << indent << ".value = ";
  if (cell.is_cartesian()) {
    out << cell.get_cartesian();
  } else if (cell.is_cut()) {
    out << cell.get_cut();
  } else {
    Igor::Panic("Unknown cell type with index {}.", cell.value.index());
  }
  out << ',' << end_char;

  out << indent << ".x_idx = " << cell.x_idx << ',' << end_char;
  out << indent << ".y_idx = " << cell.y_idx << ',' << end_char;

  out << indent << ".x_min = " << cell.x_min() << ',' << end_char;
  out << indent << ".dx = " << cell.dx() << ',' << end_char;
  out << indent << ".y_min = " << cell.y_min() << ',' << end_char;
  out << indent << ".dy = " << cell.dy() << ',' << end_char;

  out << indent << ".left_idx = ";
  if (cell.left_idx == NULL_INDEX) {
    out << "NULL";
  } else if (cell.left_idx == ZERO_FLUX_INDEX) {
    out << "ZERO_FLUX";
  } else if (cell.left_idx == SAME_VALUE_INDEX) {
    out << "SAME_VALUE";
  } else {
    out << cell.left_idx;
  }
  out << ',' << end_char;

  out << indent << ".right_idx = ";
  if (cell.right_idx == NULL_INDEX) {
    out << "NULL";
  } else if (cell.right_idx == ZERO_FLUX_INDEX) {
    out << "ZERO_FLUX";
  } else if (cell.right_idx == SAME_VALUE_INDEX) {
    out << "SAME_VALUE";
  } else {
    out << cell.right_idx;
  }
  out << ',' << end_char;

  out << indent << ".bottom_idx = ";
  if (cell.bottom_idx == NULL_INDEX) {
    out << "NULL";
  } else if (cell.bottom_idx == ZERO_FLUX_INDEX) {
    out << "ZERO_FLUX";
  } else if (cell.bottom_idx == SAME_VALUE_INDEX) {
    out << "SAME_VALUE";
  } else {
    out << cell.bottom_idx;
  }
  out << ',' << end_char;

  out << indent << ".top_idx = ";
  if (cell.top_idx == NULL_INDEX) {
    out << "NULL";
  } else if (cell.top_idx == ZERO_FLUX_INDEX) {
    out << "ZERO_FLUX";
  } else if (cell.top_idx == SAME_VALUE_INDEX) {
    out << "SAME_VALUE";
  } else {
    out << cell.top_idx;
  }
  out << ',' << end_char;

  out << "}";

  return out;
}

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_CELL_HPP_
