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
  GenCoord<ActiveFloat> rel_cut_entry{};
  Side entry_loc{};
  GenCoord<ActiveFloat> rel_cut_exit{};
  Side exit_loc{};
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
  [[nodiscard]] constexpr auto cut_entry() const noexcept
      -> CoordType2PointType<ActiveFloat, COORD_TYPE> {
    assert(is_cut());

    return CoordType2PointType<ActiveFloat, COORD_TYPE>{
        get_cut().rel_cut_entry.x * dx<COORD_TYPE>() + x_min<COORD_TYPE>(),
        get_cut().rel_cut_entry.y * dy<COORD_TYPE>() + y_min<COORD_TYPE>(),
    };
  }

  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto cut_exit() const noexcept
      -> CoordType2PointType<ActiveFloat, COORD_TYPE> {
    assert(is_cut());

    return CoordType2PointType<ActiveFloat, COORD_TYPE>{
        get_cut().rel_cut_exit.x * dx<COORD_TYPE>() + x_min<COORD_TYPE>(),
        get_cut().rel_cut_exit.y * dy<COORD_TYPE>() + y_min<COORD_TYPE>(),
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
      -> Geometry::Polygon<CoordType2PointType<ActiveFloat, COORD_TYPE>> {
    return Geometry::Polygon<CoordType2PointType<ActiveFloat, COORD_TYPE>>{{
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
    // Same order for entry_loc and exit_loc as is used in the implementation of get_subcell_points
    const auto cut_type = static_cast<uint8_t>((get_cut().exit_loc << 4) | get_cut().entry_loc);
    return get_subcell_points<COORD_TYPE>(cut_type);
  }

  // -----------------------------------------------------------------------------------------------
  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto get_right_points() const noexcept
      -> SmallVector<CoordType2PointType<ActiveFloat, COORD_TYPE>> {
    // Inverse order for entry_loc and exit_loc as is used in the implementation of
    // get_subcell_points
    const auto cut_type = static_cast<uint8_t>((get_cut().entry_loc << 4) | get_cut().exit_loc);
    return get_subcell_points<COORD_TYPE>(cut_type);
  }

  // -----------------------------------------------------------------------------------------------
  template <CoordType COORD_TYPE>
  [[nodiscard]] constexpr auto get_subcell_points(uint8_t cut_type) const noexcept
      -> SmallVector<CoordType2PointType<ActiveFloat, COORD_TYPE>> {
    switch (cut_type) {
      case make_cut_type(BOTTOM, LEFT):
        return {
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
            cut_entry<COORD_TYPE>(),
            cut_exit<COORD_TYPE>(),
        };
      case make_cut_type(LEFT, BOTTOM):
        return {
            cut_entry<COORD_TYPE>(),
            cut_exit<COORD_TYPE>(),
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };

      case make_cut_type(BOTTOM, RIGHT):
        return {
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
            cut_entry<COORD_TYPE>(),
            cut_exit<COORD_TYPE>(),
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };
      case make_cut_type(RIGHT, BOTTOM):
        return {
            cut_entry<COORD_TYPE>(),
            cut_exit<COORD_TYPE>(),
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
        };

      case make_cut_type(RIGHT, TOP):
        return {
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
            cut_entry<COORD_TYPE>(),
            cut_exit<COORD_TYPE>(),
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };
      case make_cut_type(TOP, RIGHT):
        return {
            cut_entry<COORD_TYPE>(),
            cut_exit<COORD_TYPE>(),
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };

      case make_cut_type(TOP, LEFT):
        return {
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
            cut_entry<COORD_TYPE>(),
            cut_exit<COORD_TYPE>(),
        };
      case make_cut_type(LEFT, TOP):
        return {
            cut_entry<COORD_TYPE>(),
            cut_exit<COORD_TYPE>(),
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };

      case make_cut_type(RIGHT, LEFT):
        return {
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
            cut_entry<COORD_TYPE>(),
            cut_exit<COORD_TYPE>(),
        };
      case make_cut_type(LEFT, RIGHT):
        return {
            cut_entry<COORD_TYPE>(),
            cut_exit<COORD_TYPE>(),
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };

      case make_cut_type(BOTTOM, TOP):
        return {
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>()},
            {x_min<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
            cut_entry<COORD_TYPE>(),
            cut_exit<COORD_TYPE>(),
        };

      case make_cut_type(TOP, BOTTOM):
        return {
            cut_entry<COORD_TYPE>(),
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>()},
            cut_exit<COORD_TYPE>(),
            {x_min<COORD_TYPE>() + dx<COORD_TYPE>(), y_min<COORD_TYPE>() + dy<COORD_TYPE>()},
        };

      default:
        Igor::Panic(
            "Invalid cut type with value {}. Lower bits are side `{}`, upper bits are side `{}`.",
            cut_type,
            static_cast<Side>(cut_type & 0b1111),
            static_cast<Side>(cut_type >> 4 & 0b1111));
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

  out << double_indent << ".rel_cut_entry = [" << cut_value.rel_cut_entry.x << ", "
      << cut_value.rel_cut_entry.y << ']' << end_char;
  out << double_indent << ".entry_loc = " << fmt::format("{}", cut_value.entry_loc) << end_char;
  out << double_indent << ".rel_cut_exit = [" << cut_value.rel_cut_exit.x << ", "
      << cut_value.rel_cut_exit.y << ']' << end_char;
  out << double_indent << ".exit_loc = " << fmt::format("{}", cut_value.exit_loc) << end_char;

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

template <typename ActiveFloat, typename PassiveFloat, size_t DIM>
struct fmt::formatter<Zap::CellBased::Cell<ActiveFloat, PassiveFloat, DIM>> {
  template <typename ParseContext>
  static constexpr auto parse(ParseContext& ctx) noexcept {
    return ctx.begin();
  }
  template <typename FormatContext>
  static constexpr auto format(const Zap::CellBased::Cell<ActiveFloat, PassiveFloat, DIM>& cell,
                               FormatContext& ctx) noexcept {
    std::stringstream s;
    s << cell;
    return fmt::format_to(ctx.out(), "{}", s.str());
  }
};

#endif  // ZAP_CELL_BASED_CELL_HPP_
