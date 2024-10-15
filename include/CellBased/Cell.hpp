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
#include "Igor/Logging.hpp"

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
template <typename ActiveFloat, size_t DIM>
struct CartesianValue {
  static_assert(DIM > 0, "Dimension must be at least 1.");
  Eigen::Vector<ActiveFloat, DIM> value = Eigen::Vector<ActiveFloat, DIM>::Zero();
};

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, size_t DIM>
struct CutValue {
  static_assert(DIM > 0, "Dimension must be at least 1.");
  Eigen::Vector<ActiveFloat, DIM> left_value  = Eigen::Vector<ActiveFloat, DIM>::Zero();
  Eigen::Vector<ActiveFloat, DIM> right_value = Eigen::Vector<ActiveFloat, DIM>::Zero();

  // Linear cut
  CutType type{};
  Point<ActiveFloat> rel_cut1{};  // Entry point
  Point<ActiveFloat> rel_cut2{};  // Exit point
};

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, typename PassiveFloat, size_t DIM>
struct Cell {
  static_assert(DIM > 0, "Dimension must be at least 1.");

  // Cell value
  std::variant<CartesianValue<ActiveFloat, DIM>, CutValue<ActiveFloat, DIM>> value =
      CartesianValue<ActiveFloat, DIM>{};

  // Copy the cell dimension from grid; same for all cells
  PassiveFloat m_x_min;
  PassiveFloat m_y_min;
  PassiveFloat m_dx;
  PassiveFloat m_dy;

  // Grid index
  size_t x_idx{};
  size_t y_idx{};

  // Connectivity
  size_t left_idx   = NULL_INDEX;
  size_t right_idx  = NULL_INDEX;
  size_t bottom_idx = NULL_INDEX;
  size_t top_idx    = NULL_INDEX;

  // -----------------------------------------------------------------------------------------------
  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto x_min() const noexcept -> PassiveFloat {
    if constexpr (COORD_TYPE == SIM_C) {
      return m_x_min + static_cast<PassiveFloat>(x_idx) * m_dx;
    } else {
      return static_cast<PassiveFloat>(x_idx);
    }
  }

  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto y_min() const noexcept -> PassiveFloat {
    if constexpr (COORD_TYPE == SIM_C) {
      return m_y_min + static_cast<PassiveFloat>(y_idx) * m_dy;
    } else {
      return static_cast<PassiveFloat>(y_idx);
    }
  }

  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto dx() const noexcept -> PassiveFloat {
    if constexpr (COORD_TYPE == SIM_C) {
      return m_dx;
    } else {
      return 1;
    }
  }

  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto dy() const noexcept -> PassiveFloat {
    if constexpr (COORD_TYPE == SIM_C) {
      return m_dy;
    } else {
      return 1;
    }
  }

  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto cut1() const noexcept
      -> CoordType2PointType<ActiveFloat, COORD_TYPE> {
    assert(is_cut());

    return CoordType2PointType<ActiveFloat, COORD_TYPE>{
        .x = get_cut().rel_cut1(X) * dx<COORD_TYPE>() + x_min<COORD_TYPE>(),
        .y = get_cut().rel_cut1(Y) * dy<COORD_TYPE>() + y_min<COORD_TYPE>(),
    };
  }

  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto cut2() const noexcept
      -> CoordType2PointType<ActiveFloat, COORD_TYPE> {
    assert(is_cut());

    return CoordType2PointType<ActiveFloat, COORD_TYPE>{
        .x = get_cut().rel_cut2(X) * dx<COORD_TYPE>() + x_min<COORD_TYPE>(),
        .y = get_cut().rel_cut2(Y) * dy<COORD_TYPE>() + y_min<COORD_TYPE>(),
    };
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto is_cartesian() const noexcept -> bool {
    return std::holds_alternative<CartesianValue<ActiveFloat, DIM>>(value);
  }
  [[nodiscard]] constexpr auto is_cut() const noexcept -> bool {
    return std::holds_alternative<CutValue<ActiveFloat, DIM>>(value);
  }
  [[nodiscard]] constexpr auto get_cartesian() noexcept -> CartesianValue<ActiveFloat, DIM>& {
    assert(is_cartesian());
    return std::get<CartesianValue<ActiveFloat, DIM>>(value);
  }
  [[nodiscard]] constexpr auto get_cartesian() const noexcept
      -> const CartesianValue<ActiveFloat, DIM>& {
    assert(is_cartesian());
    return std::get<CartesianValue<ActiveFloat, DIM>>(value);
  }
  [[nodiscard]] constexpr auto get_cut() noexcept -> CutValue<ActiveFloat, DIM>& {
    assert(is_cut());
    return std::get<CutValue<ActiveFloat, DIM>>(value);
  }
  [[nodiscard]] constexpr auto get_cut() const noexcept -> const CutValue<ActiveFloat, DIM>& {
    assert(is_cut());
    return std::get<CutValue<ActiveFloat, DIM>>(value);
  }

  // -----------------------------------------------------------------------------------------------
  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto get_cartesian_polygon() const noexcept
      -> Geometry::Polygon<CoordType2PointType<PassiveFloat, COORD_TYPE>> {
    return Geometry::Polygon<CoordType2PointType<PassiveFloat, COORD_TYPE>>{{
        {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
        {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
        {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
    }};
  }

  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto get_cut_left_polygon() const noexcept
      -> Geometry::Polygon<CoordType2PointType<ActiveFloat, COORD_TYPE>> {
    assert(is_cut());
    return Geometry::Polygon<CoordType2PointType<ActiveFloat, COORD_TYPE>>{
        get_left_points<COORD_TYPE>()};
  }

  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto get_cut_right_polygon() const noexcept
      -> Geometry::Polygon<CoordType2PointType<ActiveFloat, COORD_TYPE>> {
    assert(is_cut());
    return Geometry::Polygon<CoordType2PointType<ActiveFloat, COORD_TYPE>>{
        get_right_points<COORD_TYPE>()};
  }

  // -----------------------------------------------------------------------------------------------
  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto get_left_points() const noexcept
      -> SmallVector<CoordType2PointType<ActiveFloat, COORD_TYPE>> {
    const auto& cell_value = get_cut();
    switch (cell_value.type) {
      case CutType::BOTTOM_LEFT:
        return {
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
            cut1<COORD_TYPE>(),
            cut2<COORD_TYPE>(),
        };
      case CutType::BOTTOM_RIGHT:
        return {
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
            cut1<COORD_TYPE>(),
            cut2<COORD_TYPE>(),
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };
      case CutType::TOP_RIGHT:
        return {
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
            cut1<COORD_TYPE>(),
            cut2<COORD_TYPE>(),
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };
      case CutType::TOP_LEFT:
        return {
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
            cut1<COORD_TYPE>(),
            cut2<COORD_TYPE>(),
        };
      case CutType::MIDDLE_HORI:
        return {
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
            cut1<COORD_TYPE>(),
            cut2<COORD_TYPE>(),
        };
      case CutType::MIDDLE_VERT:
        return {
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
            cut1<COORD_TYPE>(),
            cut2<COORD_TYPE>(),
        };
      default:
        Igor::Panic("Unknown cut type with value {}", static_cast<int>(cell_value.type));
        std::unreachable();
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto get_right_points() const noexcept
      -> SmallVector<CoordType2PointType<ActiveFloat, COORD_TYPE>> {
    const auto& cell_value = get_cut();
    switch (cell_value.type) {
      case CutType::BOTTOM_LEFT:
        return {
            cut1<COORD_TYPE>(),
            cut2<COORD_TYPE>(),
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };
      case CutType::BOTTOM_RIGHT:
        return {
            cut1<COORD_TYPE>(),
            cut2<COORD_TYPE>(),
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
        };
      case CutType::TOP_RIGHT:
        return {
            cut1<COORD_TYPE>(),
            cut2<COORD_TYPE>(),
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };
      case CutType::TOP_LEFT:
        return {
            cut1<COORD_TYPE>(),
            cut2<COORD_TYPE>(),
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };
      case CutType::MIDDLE_HORI:
        return {
            cut1<COORD_TYPE>(),
            cut2<COORD_TYPE>(),
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };
      case CutType::MIDDLE_VERT:
        return {
            cut1<COORD_TYPE>(),
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
            cut2<COORD_TYPE>(),
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };
      default:
        Igor::Panic("Unknown cut type with value {}", static_cast<int>(cell_value.type));
        std::unreachable();
    }
  }
};

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, size_t DIM>
auto operator<<(std::ostream& out, const CartesianValue<ActiveFloat, DIM>& cart_value) noexcept
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
template <typename ActiveFloat, size_t DIM>
auto operator<<(std::ostream& out, const CutValue<ActiveFloat, DIM>& cut_value) noexcept
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

  out << double_indent << ".rel_cut1 = [" << cut_value.rel_cut1(X) << ", " << cut_value.rel_cut1(Y)
      << ']' << end_char;
  out << double_indent << ".rel_cut2 = [" << cut_value.rel_cut2(X) << ", " << cut_value.rel_cut2(Y)
      << ']' << end_char;

  out << single_indent << '}';

  return out;
}

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, typename PassiveFloat, size_t DIM>
auto operator<<(std::ostream& out, const Cell<ActiveFloat, PassiveFloat, DIM>& cell) noexcept
    -> std::ostream& {
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

  out << indent << ".x_min = " << cell.template x_min<SIM_C>() << ',' << end_char;
  out << indent << ".dx = " << cell.template dx<SIM_C>() << ',' << end_char;
  out << indent << ".y_min = " << cell.template y_min<SIM_C>() << ',' << end_char;
  out << indent << ".dy = " << cell.template dy<SIM_C>() << ',' << end_char;

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

namespace std {

template <typename ActiveFloat, typename PassiveFloat, size_t DIM>
struct formatter<Zap::CellBased::Cell<ActiveFloat, PassiveFloat, DIM>> {
  template <typename ParseContext>
  static constexpr auto parse(ParseContext& ctx) noexcept {
    return ctx.begin();
  }
  template <typename FormatContext>
  static constexpr auto format(const Zap::CellBased::Cell<ActiveFloat, PassiveFloat, DIM>& cell,
                               FormatContext& ctx) noexcept {
    std::stringstream s;
    s << cell;
    return std::format_to(ctx.out(), "{}", s.str());
  }
};

}  // namespace std

#endif  // ZAP_CELL_BASED_CELL_HPP_
