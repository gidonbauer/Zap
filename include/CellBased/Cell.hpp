#ifndef ZAP_CELL_BASED_CELL_HPP_
#define ZAP_CELL_BASED_CELL_HPP_

#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <variant>

#include <Eigen/Dense>

#include "Igor.hpp"

namespace Zap::CellBased {

// -------------------------------------------------------------------------------------------------
enum : size_t {
  SAME_VALUE_INDEX = std::numeric_limits<size_t>::max() - 2,
  ZERO_FLUX_INDEX  = std::numeric_limits<size_t>::max() - 1,
  NULL_INDEX       = std::numeric_limits<size_t>::max() - 0,
};

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
struct CartesianValue {
  static_assert(DIM > 0, "Dimension must be at least 1.");
  Eigen::Vector<Float, DIM> value = Eigen::Vector<Float, DIM>::Zero();
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
struct CutValue {
  static_assert(DIM > 0, "Dimension must be at least 1.");
  Eigen::Vector<Float, DIM> left_value  = Eigen::Vector<Float, DIM>::Zero();
  Eigen::Vector<Float, DIM> right_value = Eigen::Vector<Float, DIM>::Zero();

  // Linear cut
  CutType type{};
  Float x1_cut{};  // Entry point
  Float y1_cut{};  // Entry point
  Float x2_cut{};  // Exit point
  Float y2_cut{};  // Exit point
};

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
struct Cell {
  static_assert(DIM > 0, "Dimension must be at least 1.");

  // Cell value
  std::variant<CartesianValue<Float, DIM>, CutValue<Float, DIM>> value =
      CartesianValue<Float, DIM>{};

  // Extend
  Float x_min{};
  Float dx{};
  Float y_min{};
  Float dy{};

  // Connectivity
  size_t left_idx   = NULL_INDEX;
  size_t right_idx  = NULL_INDEX;
  size_t bottom_idx = NULL_INDEX;
  size_t top_idx    = NULL_INDEX;

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
};

// -------------------------------------------------------------------------------------------------
template <typename CellType, typename Float>
[[nodiscard]] constexpr auto
get_left_points(const CellType& cell) noexcept -> std::vector<Eigen::Vector<Float, 2>> {
  const auto& cell_value = cell.get_cut();
  switch (cell_value.type) {
    case CutType::BOTTOM_LEFT:
      return {
          Eigen::Vector<Float, 2>{cell.x_min, cell.y_min},
          Eigen::Vector<Float, 2>{cell_value.x1_cut, cell_value.y1_cut},
          Eigen::Vector<Float, 2>{cell_value.x2_cut, cell_value.y2_cut},
      };
    case CutType::BOTTOM_RIGHT:
      return {
          Eigen::Vector<Float, 2>{cell.x_min, cell.y_min},
          Eigen::Vector<Float, 2>{cell_value.x1_cut, cell_value.y1_cut},
          Eigen::Vector<Float, 2>{cell_value.x2_cut, cell_value.y2_cut},
          Eigen::Vector<Float, 2>{cell.x_min + cell.dx, cell.y_min + cell.dy},
          Eigen::Vector<Float, 2>{cell.x_min, cell.y_min + cell.dy},
      };
    case CutType::TOP_RIGHT:
      return {
          Eigen::Vector<Float, 2>{cell.x_min, cell.y_min},
          Eigen::Vector<Float, 2>{cell.x_min + cell.dx, cell.y_min},
          Eigen::Vector<Float, 2>{cell_value.x1_cut, cell_value.y1_cut},
          Eigen::Vector<Float, 2>{cell_value.x2_cut, cell_value.y2_cut},
          Eigen::Vector<Float, 2>{cell.x_min, cell.y_min + cell.dy},
      };
    case CutType::TOP_LEFT:
      return {
          Eigen::Vector<Float, 2>{cell.x_min, cell.y_min},
          Eigen::Vector<Float, 2>{cell.x_min + cell.dx, cell.y_min},
          Eigen::Vector<Float, 2>{cell.x_min + cell.dx, cell.y_min + cell.dy},
          Eigen::Vector<Float, 2>{cell_value.x1_cut, cell_value.y1_cut},
          Eigen::Vector<Float, 2>{cell_value.x2_cut, cell_value.y2_cut},
      };
    case CutType::MIDDLE_HORI:
      return {
          Eigen::Vector<Float, 2>{cell.x_min, cell.y_min},
          Eigen::Vector<Float, 2>{cell.x_min + cell.dx, cell.y_min},
          Eigen::Vector<Float, 2>{cell_value.x1_cut, cell_value.y1_cut},
          Eigen::Vector<Float, 2>{cell_value.x2_cut, cell_value.y2_cut},
      };
    case CutType::MIDDLE_VERT:
      return {
          Eigen::Vector<Float, 2>{cell.x_min, cell.y_min},
          Eigen::Vector<Float, 2>{cell.x_min, cell.y_min + cell.dy},
          Eigen::Vector<Float, 2>{cell_value.x1_cut, cell_value.y1_cut},
          Eigen::Vector<Float, 2>{cell_value.x2_cut, cell_value.y2_cut},
      };
    default:
      Igor::Panic("Unknown cut type with value {}", static_cast<int>(cell_value.type));
      std::unreachable();
  }
}

// -------------------------------------------------------------------------------------------------
template <typename CellType, typename Float>
[[nodiscard]] constexpr auto
get_right_points(const CellType& cell) noexcept -> std::vector<Eigen::Vector<Float, 2>> {
  const auto& cell_value = cell.get_cut();
  switch (cell_value.type) {
    case CutType::BOTTOM_LEFT:
      return {
          Eigen::Vector<Float, 2>{cell_value.x1_cut, cell_value.y1_cut},
          Eigen::Vector<Float, 2>{cell_value.x2_cut, cell_value.y2_cut},
          Eigen::Vector<Float, 2>{cell.x_min, cell.y_min + cell.dy},
          Eigen::Vector<Float, 2>{cell.x_min + cell.dx, cell.y_min},
          Eigen::Vector<Float, 2>{cell.x_min + cell.dx, cell.y_min + cell.dy},
      };
    case CutType::BOTTOM_RIGHT:
      return {
          Eigen::Vector<Float, 2>{cell_value.x1_cut, cell_value.y1_cut},
          Eigen::Vector<Float, 2>{cell_value.x2_cut, cell_value.y2_cut},
          Eigen::Vector<Float, 2>{cell.x_min + cell.dx, cell.y_min},
      };
    case CutType::TOP_RIGHT:
      return {
          Eigen::Vector<Float, 2>{cell_value.x1_cut, cell_value.y1_cut},
          Eigen::Vector<Float, 2>{cell_value.x2_cut, cell_value.y2_cut},
          Eigen::Vector<Float, 2>{cell.x_min + cell.dx, cell.y_min + cell.dy},
      };
    case CutType::TOP_LEFT:
      return {
          Eigen::Vector<Float, 2>{cell_value.x1_cut, cell_value.y1_cut},
          Eigen::Vector<Float, 2>{cell_value.x2_cut, cell_value.y2_cut},
          Eigen::Vector<Float, 2>{cell.x_min, cell.y_min + cell.dy},
      };
    case CutType::MIDDLE_HORI:
      return {
          Eigen::Vector<Float, 2>{cell_value.x1_cut, cell_value.y1_cut},
          Eigen::Vector<Float, 2>{cell_value.x2_cut, cell_value.y2_cut},
          Eigen::Vector<Float, 2>{cell.x_min, cell.y_min + cell.dy},
          Eigen::Vector<Float, 2>{cell.x_min + cell.dx, cell.y_min + cell.dy},
      };
    case CutType::MIDDLE_VERT:
      return {
          Eigen::Vector<Float, 2>{cell_value.x1_cut, cell_value.y1_cut},
          Eigen::Vector<Float, 2>{cell.x_min + cell.dx, cell.y_min},
          Eigen::Vector<Float, 2>{cell_value.x2_cut, cell_value.y2_cut},
          Eigen::Vector<Float, 2>{cell.x_min + cell.dx, cell.y_min + cell.dy},
      };
    default:
      Igor::Panic("Unknown cut type with value {}", static_cast<int>(cell_value.type));
      std::unreachable();
  }
}

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
auto operator<<(std::ostream& out,
                const CartesianValue<Float, DIM>& cart_value) noexcept -> std::ostream& {
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
auto operator<<(std::ostream& out,
                const CutValue<Float, DIM>& cut_value) noexcept -> std::ostream& {
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

  out << double_indent << ".x1_cut = " << cut_value.x1_cut << ',' << end_char;
  out << double_indent << ".y1_cut = " << cut_value.y1_cut << ',' << end_char;
  out << double_indent << ".x2_cut = " << cut_value.x2_cut << ',' << end_char;
  out << double_indent << ".y2_cut = " << cut_value.y2_cut << ',' << end_char;

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

  out << indent << ".x_min = " << cell.x_min << ',' << end_char;
  out << indent << ".dx = " << cell.dx << ',' << end_char;
  out << indent << ".y_min = " << cell.y_min << ',' << end_char;
  out << indent << ".dy = " << cell.dy << ',' << end_char;

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
