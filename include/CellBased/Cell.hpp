#ifndef ZAP_CELL_BASED_CELL_HPP_
#define ZAP_CELL_BASED_CELL_HPP_

#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <vector>

#include <Eigen/Dense>

#include "Igor.hpp"

namespace Zap::IO {

enum struct VTKFormat : std::uint8_t { STRUCTURED_POINTS, STRUCTURED_GRID, UNSTRUCTURED_GRID };

template <VTKFormat>
class VTKWriter;

template <typename, size_t>
class IncCellWriter;

template <typename, size_t>
class IncCellReader;

}  // namespace Zap::IO

namespace Zap::CellBased {

template <typename Float, size_t DIM>
struct CartesianCell {
  static_assert(DIM > 0, "Dimension must be at least 1.");
  enum : size_t { NULL_INDEX = std::numeric_limits<size_t>::max() };

  Eigen::Vector<Float, DIM> value{};
  Float x_min{};
  Float dx{};
  Float y_min{};
  Float dy{};
  size_t left_idx   = NULL_INDEX;
  size_t right_idx  = NULL_INDEX;
  size_t bottom_idx = NULL_INDEX;
  size_t top_idx    = NULL_INDEX;
};

template <typename Float, size_t DIM>
auto operator<<(std::ostream& out,
                const CartesianCell<Float, DIM>& cell) noexcept -> std::ostream& {
  constexpr auto end_char = '\n';
  constexpr auto indent   = "  ";

  out << "{" << end_char;

  out << indent << ".value = [ ";
  out << cell.value[0];
  for (size_t i = 1; i < DIM; ++i) {
    out << ", " << cell.value[i];
  }
  out << " ]," << end_char;

  out << indent << ".x_min = " << cell.x_min << ',' << end_char;
  out << indent << ".dx = " << cell.dx << ',' << end_char;
  out << indent << ".y_min = " << cell.y_min << ',' << end_char;
  out << indent << ".dy = " << cell.dy << ',' << end_char;

  out << indent << ".left_idx = ";
  if (cell.left_idx != cell.NULL_INDEX) {
    out << cell.left_idx;
  } else {
    out << "NULL";
  }
  out << ',' << end_char;

  out << indent << ".right_idx = ";
  if (cell.right_idx != cell.NULL_INDEX) {
    out << cell.right_idx;
  } else {
    out << "NULL";
  }
  out << ',' << end_char;

  out << indent << ".bottom_idx = ";
  if (cell.bottom_idx != cell.NULL_INDEX) {
    out << cell.bottom_idx;
  } else {
    out << "NULL";
  }
  out << ',' << end_char;

  out << indent << ".top_idx = ";
  if (cell.top_idx != cell.NULL_INDEX) {
    out << cell.top_idx;
  } else {
    out << "NULL";
  }
  out << ',' << end_char;

  out << "}";

  return out;
}

template <typename Float, size_t DIM>
class Grid {
  struct m_CornerIndices {
    size_t top_left;
    size_t top_right;
    size_t bottom_left;
    size_t bottom_right;
  };

  using m_Cell = CartesianCell<Float, DIM>;
  std::vector<m_Cell> m_cells;
  size_t m_nx;
  size_t m_ny;

  mutable std::optional<Float> m_min_delta = std::nullopt;
  mutable std::optional<Float> m_x_min     = std::nullopt;
  mutable std::optional<Float> m_x_max     = std::nullopt;
  mutable std::optional<Float> m_y_min     = std::nullopt;
  mutable std::optional<Float> m_y_max     = std::nullopt;

  // -----------------------------------------------------------------------------------------------
  constexpr Grid(size_t nx, size_t ny)
      : m_cells(nx * ny),
        m_nx(nx),
        m_ny(ny) {}

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto to_vec_idx(size_t xi, size_t yi) const noexcept -> size_t {
    assert(xi < m_nx);
    assert(yi < m_ny);
    return yi * m_nx + xi;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto to_corner_idx(size_t xi,
                                             size_t yi) const noexcept -> m_CornerIndices {
    assert(xi < m_nx);
    assert(yi < m_ny);
    return {
        .top_left     = (xi + 0) + (m_nx + 1) * (yi + 0),
        .top_right    = (xi + 1) + (m_nx + 1) * (yi + 0),
        .bottom_left  = (xi + 0) + (m_nx + 1) * (yi + 1),
        .bottom_right = (xi + 1) + (m_nx + 1) * (yi + 1),
    };
  }

 public:
  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] static constexpr auto
  Uniform(Float x_min, Float x_max, size_t nx, Float y_min, Float y_max, size_t ny) -> Grid {
    assert(nx > 0);
    assert(ny > 0);
    Grid grid(nx, ny);
    const auto dx    = (x_max - x_min) / static_cast<Float>(nx);
    const auto dy    = (y_max - y_min) / static_cast<Float>(ny);
    grid.m_min_delta = std::min(dx, dy);
    grid.m_x_min     = x_min;
    grid.m_x_max     = x_max;
    grid.m_y_min     = y_min;
    grid.m_y_max     = y_max;

    for (size_t yi = 0; yi < ny; ++yi) {
      for (size_t xi = 0; xi < nx; ++xi) {
        grid.m_cells[grid.to_vec_idx(xi, yi)] = m_Cell{
            .x_min      = x_min + static_cast<Float>(xi) * dx,
            .dx         = dx,
            .y_min      = y_min + static_cast<Float>(yi) * dy,
            .dy         = dy,
            .left_idx   = xi == 0 ? m_Cell::NULL_INDEX : grid.to_vec_idx(xi - 1, yi),
            .right_idx  = xi == (nx - 1) ? m_Cell::NULL_INDEX : grid.to_vec_idx(xi + 1, yi),
            .bottom_idx = yi == 0 ? m_Cell::NULL_INDEX : grid.to_vec_idx(xi, yi - 1),
            .top_idx    = yi == (ny - 1) ? m_Cell::NULL_INDEX : grid.to_vec_idx(xi, yi + 1),
        };
      }
    }
    return grid;
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void make_periodic() noexcept {
    for (size_t xi = 0; xi < m_nx; ++xi) {
      const auto bottom_idx          = to_vec_idx(xi, 0);
      const auto top_idx             = to_vec_idx(xi, m_ny - 1);
      m_cells[bottom_idx].bottom_idx = top_idx;
      m_cells[top_idx].top_idx       = bottom_idx;
    }
    for (size_t yi = 0; yi < m_ny; ++yi) {
      const auto left_idx          = to_vec_idx(0, yi);
      const auto right_idx         = to_vec_idx(m_nx - 1, yi);
      m_cells[left_idx].left_idx   = right_idx;
      m_cells[right_idx].right_idx = left_idx;
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <typename FUNC>
  constexpr void fill_center(FUNC f) noexcept {
    for (auto& cell : m_cells) {
      if constexpr (DIM == 1) {
        cell.value[0] = f(cell.x_min + static_cast<Float>(0.5) * cell.dx,
                          cell.y_min + static_cast<Float>(0.5) * cell.dy);
      } else {
        cell.value = f(cell.x_min + static_cast<Float>(0.5) * cell.dx,
                       cell.y_min + static_cast<Float>(0.5) * cell.dy);
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <typename FUNC>
  constexpr void fill_four_point(FUNC f) noexcept {
    for (auto& cell : m_cells) {
      if constexpr (DIM == 1) {
        cell.value(0) = (f(cell.x_min + cell.dx / 3, cell.y_min + cell.dy / 3) +
                         f(cell.x_min + 2 * cell.dx / 3, cell.y_min + cell.dy / 3) +
                         f(cell.x_min + cell.dx / 3, cell.y_min + 2 * cell.dy / 3) +
                         f(cell.x_min + 2 * cell.dx / 3, cell.y_min + 2 * cell.dy / 3)) /
                        4;
      } else {
        cell.value = (f(cell.x_min + cell.dx / 3, cell.y_min + cell.dy / 3) +
                      f(cell.x_min + 2 * cell.dx / 3, cell.y_min + cell.dy / 3) +
                      f(cell.x_min + cell.dx / 3, cell.y_min + 2 * cell.dy / 3) +
                      f(cell.x_min + 2 * cell.dx / 3, cell.y_min + 2 * cell.dy / 3)) /
                     4;
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto abs_max_value() const noexcept -> Float {
    if constexpr (DIM == 1) {
      assert(m_cells.size() > 0);
      return std::transform_reduce(
          std::cbegin(m_cells),
          std::cend(m_cells),
          m_cells[0].value[0],
          [](const auto& a, const auto& b) { return std::max(a, b); },
          [](const m_Cell& cell) { return cell.value[0]; });
    } else {
      Igor::Todo("Not implemented for DIM={}", DIM);
      return 0;
    }
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto nx() const noexcept -> size_t { return m_nx; }
  [[nodiscard]] constexpr auto ny() const noexcept -> size_t { return m_ny; }

  [[nodiscard]] constexpr auto min_delta() const noexcept -> Float {
    if (!m_min_delta.has_value()) {
      assert(m_cells.size() > 0);
      m_min_delta = std::transform_reduce(
          std::cbegin(m_cells),
          std::cend(m_cells),
          std::min(m_cells[0].dx, m_cells[0].dy),
          [](const auto& a, const auto& b) { return std::min(a, b); },
          [](const m_Cell& cell) { return std::min(cell.dx, cell.dy); });
    }
    return *m_min_delta;
  }

  [[nodiscard]] constexpr auto x_min() const noexcept -> Float {
    if (!m_x_min.has_value()) {
      Igor::Todo("Implement finding x_min.");
      m_x_min = 0;
    }
    return *m_x_min;
  }

  [[nodiscard]] constexpr auto x_max() const noexcept -> Float {
    if (!m_x_max.has_value()) {
      Igor::Todo("Implement finding x_max.");
      m_x_max = 0;
    }
    return *m_x_max;
  }

  [[nodiscard]] constexpr auto y_min() const noexcept -> Float {
    if (!m_y_min.has_value()) {
      Igor::Todo("Implement finding y_min.");
      m_y_min = 0;
    }
    return *m_y_min;
  }

  [[nodiscard]] constexpr auto y_max() const noexcept -> Float {
    if (!m_y_max.has_value()) {
      Igor::Todo("Implement finding y_max.");
      m_y_max = 0;
    }
    return *m_y_max;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto operator[](size_t idx) const noexcept -> const m_Cell& {
    assert(idx < m_nx * m_ny);
    return m_cells[idx];
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void dump_cells(std::ostream& out) const noexcept {
    for (const auto& cell : m_cells) {
      out << cell << ",\n";
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <typename F, typename G>
  friend class Solver;

  template <Zap::IO::VTKFormat F>
  friend class Zap::IO::VTKWriter;

  template <typename F, size_t D>
  friend class Zap::IO::IncCellWriter;

  template <typename F, size_t D>
  friend class Zap::IO::IncCellReader;
};

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_CELL_HPP_
