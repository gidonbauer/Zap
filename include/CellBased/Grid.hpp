#ifndef ZAP_CELL_BASED_GRID_HPP_
#define ZAP_CELL_BASED_GRID_HPP_

#include "CellBased/Cell.hpp"
#include "CellBased/Definitions.hpp"
#include "CellBased/GridHelper.hpp"
#include "CellBased/Quadrature.hpp"

#include "Igor/TypeName.hpp"

namespace Zap::CellBased {

template <typename ActiveFloat, typename PassiveFloat>
class UniformGrid {
  enum { ON_MIN = -1, NOT_ON = 0, ON_MAX = 1 };

  std::vector<Cell<ActiveFloat, PassiveFloat>> m_cells;
  size_t m_nx;
  size_t m_ny;

  PassiveFloat m_x_min;
  PassiveFloat m_x_max;
  PassiveFloat m_y_min;
  PassiveFloat m_y_max;
  PassiveFloat m_dx;
  PassiveFloat m_dy;

  std::vector<size_t> m_cut_cell_idxs;

 public:
  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto to_vec_idx(size_t xi, size_t yi) const noexcept -> size_t {
    IGOR_ASSERT(xi < m_nx, "Index {} is out of bounds with max {}.", xi, m_nx);
    IGOR_ASSERT(yi < m_ny, "Index {} is out of bounds with max {}.", yi, m_ny);
    return yi * m_nx + xi;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto to_mat_idx(size_t i) const noexcept -> std::pair<size_t, size_t> {
    IGOR_ASSERT(i < m_cells.size(), "Index {} is out of bounds with max {}.", i, m_cells.size());
    return std::make_pair(i % m_nx, i / m_nx);
  }

  // -----------------------------------------------------------------------------------------------
  template <typename T>
  [[nodiscard]] constexpr auto to_grid_coord(const SimCoord<T>& simulation_coord) const noexcept
      -> GridCoord<T> {
    return {
        (simulation_coord.x - m_x_min) / m_dx,
        (simulation_coord.y - m_y_min) / m_dy,
    };
  }

  // -----------------------------------------------------------------------------------------------
  template <typename T>
  [[nodiscard]] constexpr auto to_sim_coord(const GridCoord<T>& grid_coord) const noexcept
      -> SimCoord<T> {
    return {
        grid_coord.x * m_dx + m_x_min,
        grid_coord.y * m_dy + m_y_min,
    };
  }

  // -----------------------------------------------------------------------------------------------
  constexpr UniformGrid(PassiveFloat x_min,
                        PassiveFloat x_max,
                        size_t nx,
                        PassiveFloat y_min,
                        PassiveFloat y_max,
                        size_t ny) noexcept
      : m_cells(nx * ny),
        m_nx(nx),
        m_ny(ny),
        m_x_min(x_min),
        m_x_max(x_max),
        m_y_min(y_min),
        m_y_max(y_max),
        m_dx((x_max - x_min) / static_cast<PassiveFloat>(nx)),
        m_dy((y_max - y_min) / static_cast<PassiveFloat>(ny)) {
    IGOR_ASSERT(nx > 0, "nx must be greater than zero.");
    IGOR_ASSERT(ny > 0, "ny must be greater than zero.");

    for (size_t yi = 0; yi < ny; ++yi) {
      for (size_t xi = 0; xi < nx; ++xi) {
        // clang-format off
        m_cells[to_vec_idx(xi, yi)] = Cell<ActiveFloat, PassiveFloat>{
            .m_x_min    = m_x_min,
            .m_y_min    = m_y_min,
            .m_dx       = m_dx,
            .m_dy       = m_dy,
            .x_idx      = xi,
            .y_idx      = yi,
            .left_idx   = xi == 0        ? NULL_INDEX : to_vec_idx(xi - 1, yi),
            .right_idx  = xi == (nx - 1) ? NULL_INDEX : to_vec_idx(xi + 1, yi),
            .bottom_idx = yi == 0        ? NULL_INDEX : to_vec_idx(xi, yi - 1),
            .top_idx    = yi == (ny - 1) ? NULL_INDEX : to_vec_idx(xi, yi + 1),
        };
        // clang-format on
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void periodic_boundary(std::underlying_type_t<Side> sides = ALL) noexcept {
    // Left side
    if ((sides & LEFT) != 0) {
      for (size_t yi = 0; yi < m_ny; ++yi) {
        const auto left_idx        = to_vec_idx(0, yi);
        const auto right_idx       = to_vec_idx(m_nx - 1, yi);
        m_cells[left_idx].left_idx = right_idx;
      }
    }

    // Right side
    if ((sides & RIGHT) != 0) {
      for (size_t yi = 0; yi < m_ny; ++yi) {
        const auto left_idx          = to_vec_idx(0, yi);
        const auto right_idx         = to_vec_idx(m_nx - 1, yi);
        m_cells[right_idx].right_idx = left_idx;
      }
    }

    // Bottom side
    if ((sides & BOTTOM) != 0) {
      for (size_t xi = 0; xi < m_nx; ++xi) {
        const auto bottom_idx          = to_vec_idx(xi, 0);
        const auto top_idx             = to_vec_idx(xi, m_ny - 1);
        m_cells[bottom_idx].bottom_idx = top_idx;
      }
    }

    // Top side
    if ((sides & TOP) != 0) {
      for (size_t xi = 0; xi < m_nx; ++xi) {
        const auto bottom_idx    = to_vec_idx(xi, 0);
        const auto top_idx       = to_vec_idx(xi, m_ny - 1);
        m_cells[top_idx].top_idx = bottom_idx;
      }
    }
  }

 private:
  // -----------------------------------------------------------------------------------------------
  constexpr void set_boundary_impl(std::underlying_type_t<Side> sides, size_t idx) noexcept {
    // Left side
    if ((sides & LEFT) != 0) {
      for (size_t yi = 0; yi < m_ny; ++yi) {
        const auto left_idx        = to_vec_idx(0, yi);
        m_cells[left_idx].left_idx = idx;
      }
    }

    // Right side
    if ((sides & RIGHT) != 0) {
      for (size_t yi = 0; yi < m_ny; ++yi) {
        const auto right_idx         = to_vec_idx(m_nx - 1, yi);
        m_cells[right_idx].right_idx = idx;
      }
    }

    // Bottom side
    if ((sides & BOTTOM) != 0) {
      for (size_t xi = 0; xi < m_nx; ++xi) {
        const auto bottom_idx          = to_vec_idx(xi, 0);
        m_cells[bottom_idx].bottom_idx = idx;
      }
    }

    // Top side
    if ((sides & TOP) != 0) {
      for (size_t xi = 0; xi < m_nx; ++xi) {
        const auto top_idx       = to_vec_idx(xi, m_ny - 1);
        m_cells[top_idx].top_idx = idx;
      }
    }
  }

 public:
  // -----------------------------------------------------------------------------------------------
  constexpr void zero_flux_boundary(std::underlying_type_t<Side> sides = ALL) noexcept {
    set_boundary_impl(sides, ZERO_FLUX_INDEX);
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void same_value_boundary(std::underlying_type_t<Side> sides = ALL) noexcept {
    set_boundary_impl(sides, SAME_VALUE_INDEX);
  }

  // -----------------------------------------------------------------------------------------------
  template <typename FUNC>
  constexpr void fill_center(FUNC f) noexcept {
    const auto get_centroid =
        []<typename PointType>(const SmallVector<PointType>& points) -> PointType {
      auto centroid = PointType::Zero();
      for (const auto& p : points) {
        centroid += p;
      }
      return centroid / static_cast<PassiveFloat>(points.size());
    };

    for (auto& cell : m_cells) {
      if (cell.is_cartesian()) {
        auto& value = cell.get_cartesian().value;
        value       = f(cell.template x_min<SIM_C>() +
                      static_cast<PassiveFloat>(0.5) * cell.template dx<SIM_C>(),
                  cell.template y_min<SIM_C>() +
                      static_cast<PassiveFloat>(0.5) * cell.template dy<SIM_C>());
      } else if (cell.is_cut()) {
        const auto left_centroid  = get_centroid(cell.template get_left_points<SIM_C>());
        const auto right_centroid = get_centroid(cell.template get_right_points<SIM_C>());

        auto& cell_value       = cell.get_cut();
        cell_value.left_value  = f(left_centroid.x, left_centroid.y);
        cell_value.right_value = f(right_centroid.x, right_centroid.y);
      } else {
        Igor::Panic("Unknown cell type with variant index {}.", cell.cell_type.index());
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <typename FUNC>
  constexpr void fill_four_point(FUNC f) noexcept {
    for (auto& cell : m_cells) {
      if (cell.is_cartesian()) {
        auto& value = cell.get_cartesian().value;
        value       = (f(cell.template x_min<SIM_C>() + cell.template dx<SIM_C>() / 4,
                   cell.template y_min<SIM_C>() + cell.template dy<SIM_C>() / 4) +
                 f(cell.template x_min<SIM_C>() + 3 * cell.template dx<SIM_C>() / 4,
                   cell.template y_min<SIM_C>() + cell.template dy<SIM_C>() / 4) +
                 f(cell.template x_min<SIM_C>() + cell.template dx<SIM_C>() / 4,
                   cell.template y_min<SIM_C>() + 3 * cell.template dy<SIM_C>() / 4) +
                 f(cell.template x_min<SIM_C>() + 3 * cell.template dx<SIM_C>() / 4,
                   cell.template y_min<SIM_C>() + 3 * cell.template dy<SIM_C>() / 4)) /
                4;
      } else if (cell.is_cut()) {
        auto& cut_value = cell.get_cut();

        // Left side
        {
          const auto corners = cell.template get_left_points<SIM_C>();
          for (size_t j = 0; j < corners.size(); ++j) {
            auto point = SimCoord<ActiveFloat>::Zero();
            for (size_t i = 0; i < corners.size(); ++i) {
              point += (1 + (i == j)) * corners[i] / static_cast<PassiveFloat>(corners.size() + 1);
            }

            cut_value.left_value += f(point.x, point.y);
          }
          cut_value.left_value /= static_cast<PassiveFloat>(corners.size());
        }

        // Right side
        {
          const auto corners = cell.template get_right_points<SIM_C>();
          for (size_t j = 0; j < corners.size(); ++j) {
            auto point = SimCoord<ActiveFloat>::Zero();
            for (size_t i = 0; i < corners.size(); ++i) {
              point += (1 + (i == j)) * corners[i] / static_cast<PassiveFloat>(corners.size() + 1);
            }

            cut_value.right_value += f(point.x, point.y);
          }
          cut_value.right_value /= static_cast<PassiveFloat>(corners.size());
        }
      } else {
        Igor::Panic("Unknown cell type with variant index {}.", cell.cell_type.index());
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <size_t N = 15UZ, typename FUNC>
  constexpr void fill_quad(FUNC f) noexcept {
    for (auto& cell : m_cells) {
      if (cell.is_cartesian()) {
        const auto polygon         = cell.template get_cartesian_polygon<SIM_C>();
        const auto integral        = quadrature<N>(f, polygon);
        cell.get_cartesian().value = integral / polygon.area();
      } else {
        // Left subcell
        {
          const auto polygon        = cell.template get_cut_left_polygon<SIM_C>();
          const auto integral       = quadrature<N>(f, polygon);
          cell.get_cut().left_value = integral / polygon.area();
        }
        // Right subcell
        {
          const auto polygon         = cell.template get_cut_right_polygon<SIM_C>();
          const auto integral        = quadrature<N>(f, polygon);
          cell.get_cut().right_value = integral / polygon.area();
        }
      }
    }
  }

  // -------------------------------------------------------------------------------------------------
  template <typename PointType>
  [[nodiscard]] constexpr auto
  find_sides(const PointType& p, const Cell<ActiveFloat, PassiveFloat>& cell) const noexcept
      -> Side {
    constexpr CoordType coord_type = PointType2CoordType<PointType>;

    const int on_x = approx_eq(cell.template x_min<coord_type>(), ad::value(p.x)) * ON_MIN +
                     approx_eq(cell.template x_min<coord_type>() + cell.template dx<coord_type>(),
                               ad::value(p.x)) *
                         ON_MAX;
    const int on_y = approx_eq(cell.template y_min<coord_type>(), ad::value(p.y)) * ON_MIN +
                     approx_eq(cell.template y_min<coord_type>() + cell.template dy<coord_type>(),
                               ad::value(p.y)) *
                         ON_MAX;
    IGOR_ASSERT(on_x != NOT_ON || on_y != NOT_ON,
                "Point {} is not on a edge in x-direction or in y-direction of cell {}.",
                p,
                cell);
    Side loc = static_cast<Side>((on_x == ON_MIN) * LEFT | (on_x == ON_MAX) * RIGHT |
                                 (on_y == ON_MIN) * BOTTOM | (on_y == ON_MAX) * TOP);
    IGOR_ASSERT(count_sides(loc) <= 2,
                "The cut can be on a maximum of two sides, but is on {}.",
                count_sides(loc));
    return loc;
  };

  // -----------------------------------------------------------------------------------------------
  enum class ClassifyError : uint8_t { OK, BOTH_CUTS_ON_SAME_SIDE, CUT_ON_MULTIPLE_SIDES };

  template <Point2D_c PointType>
  [[nodiscard]] constexpr auto classify_cut(const Cell<ActiveFloat, PassiveFloat>& cell,
                                            PointType& cut_entry_point,
                                            PointType& cut_exit_point,
                                            Side& cut_entry_loc,
                                            Side& cut_exit_loc) const noexcept -> ClassifyError {
    constexpr CoordType coord_type = PointType2CoordType<PointType>;

    IGOR_ASSERT(point_in_cell(cut_entry_point, cell),
                "Cut-entry-point ({}) {} is not in cell {}.",
                Igor::type_name(cut_entry_point),
                cut_entry_point,
                cell);
    IGOR_ASSERT(point_in_cell(cut_exit_point, cell),
                "Cut-exit-point ({}) {} is not in cell {}.",
                Igor::type_name(cut_exit_point),
                cut_exit_point,
                cell);

    auto remove_common = [](Side& side1, Side& side2) {
      const Side common = static_cast<Side>(side1 & side2);
      side1             = static_cast<Side>(side1 ^ common);
      side2             = static_cast<Side>(side2 ^ common);
    };

    auto decide_on_corner = [](Side& side) {
      if (side == static_cast<Side>(LEFT | TOP) || side == static_cast<Side>(RIGHT | TOP)) {
        side = TOP;
      } else if (side == static_cast<Side>(LEFT | BOTTOM) ||
                 side == static_cast<Side>(RIGHT | BOTTOM)) {
        side = BOTTOM;
      }
    };

    cut_entry_loc = find_sides(cut_entry_point, cell);
    cut_exit_loc  = find_sides(cut_exit_point, cell);

    if (cut_entry_loc == cut_exit_loc) { return ClassifyError::BOTH_CUTS_ON_SAME_SIDE; }

    remove_common(cut_entry_loc, cut_exit_loc);
    decide_on_corner(cut_entry_loc);
    decide_on_corner(cut_exit_loc);

    if (count_sides(cut_entry_loc) == 0 || count_sides(cut_exit_loc) == 0) {
      return ClassifyError::BOTH_CUTS_ON_SAME_SIDE;
    }

    if (count_sides(cut_entry_loc) != 1 || count_sides(cut_exit_loc) != 1) {
      Igor::Warn("Expected cut1 on cut2 to be on only one side after removing common sides, but "
                 "cut_entry_loc = {} and cut_exit_loc = {}.",
                 cut_entry_loc,
                 cut_exit_loc);
      return ClassifyError::CUT_ON_MULTIPLE_SIDES;
    }

    // Make cut relative to cell size, i.e. cut in [0, 1]
    cut_entry_point.x =
        (cut_entry_point.x - cell.template x_min<coord_type>()) / cell.template dx<coord_type>();
    cut_entry_point.y =
        (cut_entry_point.y - cell.template y_min<coord_type>()) / cell.template dy<coord_type>();
    cut_exit_point.x =
        (cut_exit_point.x - cell.template x_min<coord_type>()) / cell.template dx<coord_type>();
    cut_exit_point.y =
        (cut_exit_point.y - cell.template y_min<coord_type>()) / cell.template dy<coord_type>();

    IGOR_ASSERT(cut_entry_point.x >= -EPS<PassiveFloat>() &&
                    cut_entry_point.x - 1 <= EPS<PassiveFloat>(),
                "Cut-entry-x {} is not in [0, 1].",
                cut_entry_point.x);
    IGOR_ASSERT(cut_entry_point.y >= -EPS<PassiveFloat>() &&
                    cut_entry_point.y - 1 <= EPS<PassiveFloat>(),
                "Cut-entry-y {} is not in [0, 1].",
                cut_entry_point.y);
    IGOR_ASSERT(cut_exit_point.x >= -EPS<PassiveFloat>() &&
                    cut_exit_point.x - 1 <= EPS<PassiveFloat>(),
                "Cut-exit-x {} is not in [0, 1].",
                cut_exit_point.x);
    IGOR_ASSERT(cut_exit_point.y >= -EPS<PassiveFloat>() &&
                    cut_exit_point.y - 1 <= EPS<PassiveFloat>(),
                "Cut-exit-y {} is not in [0, 1].",
                cut_exit_point.y);

    return ClassifyError::OK;
  }

  // -----------------------------------------------------------------------------------------------
  template <typename ParamCurve>
  [[nodiscard]] constexpr auto cut_curve(ParamCurve curve,
                                         size_t num_discretization_points = 0UZ) noexcept -> bool {
    if (num_discretization_points == 0UZ) { num_discretization_points = size() / 2UZ; }

    std::vector<GridCoord<ActiveFloat>> points(num_discretization_points);
    std::generate(std::begin(points),
                  std::end(points),
                  [this, num_discretization_points, curve, i = 0UZ]() mutable {
                    const auto t = static_cast<PassiveFloat>(i++) /
                                   static_cast<PassiveFloat>(num_discretization_points - 1);
                    return to_grid_coord(curve(t));
                  });
    return cut_piecewise_linear<ExtendType::MAX>(points);
  }

  // -----------------------------------------------------------------------------------------------
  template <ExtendType extend_type, typename T>
  [[nodiscard]] auto cut_piecewise_linear(const std::vector<SimCoord<T>>& points) noexcept -> bool {
    std::vector<GridCoord<T>> points_grid_coord(points.size());
    std::transform(std::cbegin(points),
                   std::cend(points),
                   std::begin(points_grid_coord),
                   [this](const SimCoord<T>& p) { return to_grid_coord(p); });
    return cut_piecewise_linear<extend_type>(points_grid_coord);
  }

  // -----------------------------------------------------------------------------------------------
  template <ExtendType extend_type, typename T>
  [[nodiscard]] auto cut_piecewise_linear(std::vector<GridCoord<T>> points) noexcept -> bool {
    if (!m_cut_cell_idxs.empty()) {
      Igor::Warn("Grid is already cut, this method expects the grid to no be cut.");
      return false;
    }

    auto get_whole_min = [](T x1, T x2) -> int {
      const T x = std::min(x1, x2);
      if (x - std::floor(x) < EPS<T>()) { return static_cast<int>(std::floor(x)); }
      return static_cast<int>(std::ceil(x));
    };
    auto get_whole_max = [](T x1, T x2) -> int {
      const T x = std::max(x1, x2);
      if (std::ceil(x) - x < EPS<T>()) { return static_cast<int>(std::ceil(x)); }
      return static_cast<int>(std::floor(x));
    };

    // Remove points outside of grid
    std::erase_if(points, [this](const GridCoord<T>& p) { return !point_in_grid(p); });
    if (points.size() < 2) {
      Igor::Warn("Only {} points inside of the grid, require at least 2.", points.size());
      return false;
    }

    // - Extend end points -------------------------------------------------------------------------
    if constexpr (extend_type != ExtendType::NONE) {
      // - Handle first point; extend curve backwards to cut cell containing the first point -------
      {
        auto& p0       = points[0];
        const auto& p1 = points[1];

        const auto s = p1 - p0;
        const std::array<T, 4> rs{
            (std::floor(p0.x) - p0.x) / s.x,
            (std::ceil(p0.x) - p0.x) / s.x,
            (std::floor(p0.y) - p0.y) / s.y,
            (std::ceil(p0.y) - p0.y) / s.y,
        };

        T r_entry = -std::numeric_limits<T>::max();
        for (auto r : rs) {
          if constexpr (extend_type == ExtendType::MAX) {
            if (r < EPS<T>() && r > r_entry) { r_entry = r; }
          } else {
            if (std::abs(r) < std::abs(r_entry)) { r_entry = r; }
          }
        }
        IGOR_ASSERT(r_entry != -std::numeric_limits<T>::max(),
                    "Could not find the extended intersection point for p0={} and p1={}. Potential "
                    "r-values are rs={}",
                    p0,
                    p1,
                    rs);
        p0 += r_entry * s;
      }

      // - Handle last point; extend curve forwards to cut cell containing the last point ----------
      {
        const auto& p0 = points[points.size() - 2];
        auto& p1       = points[points.size() - 1];

        const auto s = p1 - p0;
        const std::array<T, 4> rs{
            (std::floor(p1.x) - p1.x) / s.x,
            (std::ceil(p1.x) - p1.x) / s.x,
            (std::floor(p1.y) - p1.y) / s.y,
            (std::ceil(p1.y) - p1.y) / s.y,
        };

        T r_exit = std::numeric_limits<T>::max();
        for (auto r : rs) {
          if constexpr (extend_type == ExtendType::MAX) {
            if (r > -EPS<T>() && r < r_exit) { r_exit = r; }
          } else {
            if (std::abs(r) < std::abs(r_exit)) { r_exit = r; }
          }
        }
        IGOR_ASSERT(r_exit != std::numeric_limits<T>::max(),
                    "Could not find the extended intersection point for p0={} and p1={}. Potential "
                    "r-values are rs={}",
                    p0,
                    p1,
                    rs);

        p1 += r_exit * s;
      }
    }
    // - Extend end points -------------------------------------------------------------------------

    std::vector<GridCoord<T>> intersect_points{};
    ptrdiff_t begin_points_idx = 0;
    for (size_t i = 0; i < points.size() - 1; ++i) {
      const auto& p0 = points[i];
      const auto& p1 = points[i + 1];

      const auto slope = p1 - p0;

      const bool search_whole_x = !approx_eq(slope.x, static_cast<T>(0));
      const bool search_whole_y = !approx_eq(slope.y, static_cast<T>(0));
      if (!(search_whole_x || search_whole_y)) {
        Igor::Warn("p0 = {} and p1 = {} are the same point.", p0, p1);
        continue;
      }

      // Intersections where x in N
      if (search_whole_x) {
        for (int x = get_whole_min(p0.x, p1.x); x <= get_whole_max(p0.x, p1.x); ++x) {
          const T r = (static_cast<T>(x) - p0.x) / slope.x;
          const T y = p0.y + r * slope.y;
          intersect_points.emplace_back(static_cast<T>(x), y);
        }
      }

      // Intersections where y in N
      if (search_whole_y) {
        for (int y = get_whole_min(p0.y, p1.y); y <= get_whole_max(p0.y, p1.y); ++y) {
          const T r = (static_cast<T>(y) - p0.y) / slope.y;
          const T x = p0.x + r * slope.x;

          // If x is a whole number, the point was already added in the previous loop, can ignore it
          if (search_whole_x && approx_eq<PassiveFloat>(ad::value(x), std::round(ad::value(x)))) {
            continue;
          }

          intersect_points.emplace_back(x, static_cast<T>(y));
        }
      }

      std::sort(std::next(std::begin(intersect_points), begin_points_idx),
                std::end(intersect_points),
                [&p0](const GridCoord<T>& lhs, const GridCoord<T>& rhs) {
                  const auto lhs_dist = (lhs - p0).norm();
                  const auto rhs_dist = (rhs - p0).norm();
                  return lhs_dist < rhs_dist;
                });
      begin_points_idx = static_cast<ptrdiff_t>(intersect_points.size());
    }

    // Remove duplicate points
    // TODO: Can we do this on the fly while filling intersect_points?
    intersect_points.erase(std::unique(intersect_points.begin(),
                                       intersect_points.end(),
                                       [](const GridCoord<T>& lhs, const GridCoord<T>& rhs) {
                                         return (lhs - rhs).norm() < EPS<T>();
                                       }),
                           intersect_points.end());

    m_cut_cell_idxs.reserve(intersect_points.size() - 1);
    for (size_t i = 0; i < intersect_points.size() - 1;) {
      auto entry_point = intersect_points[i];
      auto exit_point  = intersect_points[i + 1];
      if (approx_eq((entry_point - exit_point).norm(), static_cast<T>(0))) {
        Igor::Warn("entry_point {} and exit_point {} are the same point, ignore.",
                   entry_point,
                   exit_point);
        continue;
      }

      const auto mid_point = (entry_point + exit_point) / 2;
      const auto cell_idx  = find_cell(mid_point);
      IGOR_ASSERT(cell_idx != NULL_INDEX,
                  "Expected mid_point {} to be in grid [{}, {}]x[{}, {}], but is not.",
                  mid_point,
                  0,
                  m_nx,
                  0,
                  m_ny);

      auto& cell_to_cut = m_cells[cell_idx];
      Side entry_loc;
      Side exit_loc;
      switch (classify_cut(cell_to_cut, entry_point, exit_point, entry_loc, exit_loc)) {
        case ClassifyError::BOTH_CUTS_ON_SAME_SIDE:
          {
            const auto num_sides_entry = count_sides(find_sides(entry_point, cell_to_cut));
            const auto num_sides_exit  = count_sides(find_sides(exit_point, cell_to_cut));
            IGOR_ASSERT(!(num_sides_entry == 2 && num_sides_exit == 2),
                        "The case where both intersection points lay on a corner should be handled "
                        "correctly.");

            // Entry point lays on single side => remove entry and exit point and clean up last
            // cut-cell
            if (num_sides_entry == 1 && num_sides_exit == 1) {
              // Remove intersection points that are on the same side of the cell
              intersect_points.erase(
                  std::next(std::begin(intersect_points), static_cast<ptrdiff_t>(i)),
                  std::next(std::begin(intersect_points), static_cast<ptrdiff_t>(i + 2)));

              // Reset i to continue at the correct index
              // TODO: This assert fails very often
              IGOR_ASSERT(i >= 1,
                          "Index must be greater or equal to 1, but is {}. I don't think we can do "
                          "anything in this case.",
                          i);
              i -= 1;

              // Repair cut cells by removing last cut as it was incorrect
              const auto last_cut_idx = m_cut_cell_idxs.back();
              m_cut_cell_idxs.pop_back();
              const auto old_value = m_cells[last_cut_idx].get_cut().left_value;  // == right_value
              m_cells[last_cut_idx].cell_type = CartesianValue<ActiveFloat>{.value = old_value};
            }
            // Entry point lays on corner => remove only exit point and leave last cut-cell
            else if (num_sides_entry == 2) {
              // Remove exit point
              intersect_points.erase(
                  std::next(std::begin(intersect_points), static_cast<ptrdiff_t>(i + 1)),
                  std::next(std::begin(intersect_points), static_cast<ptrdiff_t>(i + 2)));
            } else if (num_sides_exit == 2) {
              // Remove entry point
              intersect_points.erase(
                  std::next(std::begin(intersect_points), static_cast<ptrdiff_t>(i)),
                  std::next(std::begin(intersect_points), static_cast<ptrdiff_t>(i + 1)));

              // Reset i to continue at the correct index
              IGOR_ASSERT(i >= 1,
                          "Index must be greater or equal to 1, but is {}. I don't think we can do "
                          "anything in this case.",
                          i);
              i -= 1;

              // Repair cut cells by removing last cut as it was incorrect
              const auto last_cut_idx = m_cut_cell_idxs.back();
              m_cut_cell_idxs.pop_back();
              const auto old_value = m_cells[last_cut_idx].get_cut().left_value;  // == right_value
              m_cells[last_cut_idx].cell_type = CartesianValue<ActiveFloat>{.value = old_value};
            }

            continue;
          }
        case ClassifyError::CUT_ON_MULTIPLE_SIDES:
          {
            Igor::Warn("Could not classify the cut, cut is invalid.");
            Igor::Debug("Points used to cut the grid piecewise linear: {}", points);
            Igor::Debug("entry_point = {}", entry_point);
            Igor::Debug("exit_point = {}", exit_point);
            return false;
          }
        case ClassifyError::OK: break;
      }

      ActiveFloat old_value = 0.0;
      if (cell_to_cut.is_cartesian()) { old_value = cell_to_cut.get_cartesian().value; }

      cell_to_cut.cell_type = CutValue<ActiveFloat>{
          .left_value    = old_value,
          .right_value   = old_value,
          .rel_cut_entry = {entry_point.x, entry_point.y},
          .entry_loc     = entry_loc,
          .rel_cut_exit  = {exit_point.x, exit_point.y},
          .exit_loc      = exit_loc,
      };

      m_cut_cell_idxs.push_back(cell_idx);
      ++i;
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void merge_cut_cells() noexcept {
    for (size_t cell_idx : m_cut_cell_idxs) {
      IGOR_ASSERT(is_cell(cell_idx), "Expected cell {} to be a valid cell, but is not.", cell_idx);
      auto& cell = m_cells[cell_idx];
      IGOR_ASSERT(
          cell.is_cut(), "Expected cell {} to be cut, but is cartesian: {}", cell_idx, cell);

      const auto left_value   = cell.get_cut().left_value;
      const auto left_polygon = cell.template get_cut_left_polygon<SIM_C>();

      const auto right_value   = cell.get_cut().right_value;
      const auto right_polygon = cell.template get_cut_right_polygon<SIM_C>();

      const auto cell_polygon = cell.template get_cartesian_polygon<SIM_C>();

      cell.cell_type = CartesianValue<ActiveFloat>{
          .value = (left_value * left_polygon.area() + right_value * right_polygon.area()) /
                   cell_polygon.area()};
    }
    m_cut_cell_idxs.clear();
  }

  // -----------------------------------------------------------------------------------------------
  template <typename PointType>
  [[nodiscard]] constexpr auto point_in_grid(const PointType& point) const noexcept -> bool {
    if constexpr (is_SimCoord_v<PointType>) {
      return point.x - m_x_min >= -EPS<PassiveFloat>() &&
             m_x_max - point.x >= -EPS<PassiveFloat>() &&  //
             point.y - m_y_min >= -EPS<PassiveFloat>() && m_y_max - point.y >= -EPS<PassiveFloat>();
    } else {
      return point.x >= -EPS<PassiveFloat>() &&
             static_cast<PassiveFloat>(m_nx) - point.x >= -EPS<PassiveFloat>() &&  //
             point.y >= -EPS<PassiveFloat>() &&
             static_cast<PassiveFloat>(m_ny) - point.y >= -EPS<PassiveFloat>();
    }
  }

  // -----------------------------------------------------------------------------------------------
  template <Point2D_c PointType>
  [[nodiscard]] constexpr auto find_cell(const PointType& pos) const noexcept -> size_t {
    const auto gpos = [this, &pos] {
      if constexpr (is_SimCoord_v<PointType>) {
        return to_grid_coord(pos);
      } else {
        (void)this;
        return pos;
      }
    }();

    if (gpos.x < -EPS<PassiveFloat>() ||
        gpos.x > static_cast<PassiveFloat>(m_nx) + EPS<PassiveFloat>()) {
      Igor::Warn(
          "x-component of position Sim({}) ^= Grid({}) is not in grid with nx={}", pos, gpos, m_nx);
      return NULL_INDEX;
    }
    if (gpos.y < -EPS<PassiveFloat>() ||
        gpos.y > static_cast<PassiveFloat>(m_ny) + EPS<PassiveFloat>()) {
      Igor::Warn(
          "y-component of position Sim({}) ^= Grid({}) is not in grid with ny={}", pos, gpos, m_ny);
      return NULL_INDEX;
    }

    const size_t xi = [this, &gpos] {
      if (approx_eq(ad::value(gpos.x), static_cast<PassiveFloat>(m_nx))) {
        return m_nx - 1;
      } else {
        return static_cast<size_t>(std::floor(gpos.x));
      }
    }();

    const size_t yi = [this, &gpos] {
      if (approx_eq(ad::value(gpos.y), static_cast<PassiveFloat>(m_ny))) {
        return m_ny - 1;
      } else {
        return static_cast<size_t>(std::floor(gpos.y));
      }
    }();

    IGOR_ASSERT(xi < m_nx,
                "gpos.x={} translates to xi={} which is out of bounds for nx={}.",
                gpos.x,
                xi,
                m_nx);
    IGOR_ASSERT(yi < m_ny,
                "gpos.y={} translates to yi={} which is out of bounds for ny={}.",
                gpos.y,
                yi,
                m_ny);

    return to_vec_idx(xi, yi);
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto get_shock_curve() const noexcept
      -> std::vector<SimCoord<ActiveFloat>> {
    if (m_cut_cell_idxs.empty()) { return {}; }

    std::vector<SimCoord<ActiveFloat>> curve_points(m_cut_cell_idxs.size() + 1);

    // - First point on shock curve ----------------------------------------------------------------
    {
      const auto& cell = m_cells[m_cut_cell_idxs.front()];
      curve_points[0]  = cell.template cut_entry<SIM_C>();
    }

    // - Intermediate points on shock curve --------------------------------------------------------
    for (size_t i = 0; i < m_cut_cell_idxs.size() - 1; ++i) {
      const auto& cell_prev = m_cells[m_cut_cell_idxs[i]];
      const auto& cell_next = m_cells[m_cut_cell_idxs[i + 1]];

      const auto cut_prev = cell_prev.template cut_exit<SIM_C>();
      const auto cut_next = cell_next.template cut_entry<SIM_C>();
      IGOR_ASSERT((cut_prev - cut_next).norm() < EPS<PassiveFloat>(),
                  "Expected cut_a and cut_b to be the same point but cut_a = {} and cut_b = "
                  "{}.\ncell_a = {}\ncell_b = {}",
                  cut_prev,
                  cut_next,
                  cell_prev,
                  cell_next);

      // This is done purely to mix together the derivatives, although they might be the same
      curve_points[i + 1] = (cut_prev + cut_next) / 2;
    }

    // - Last point on shock curve -----------------------------------------------------------------
    {
      const auto& cell    = m_cells[m_cut_cell_idxs.back()];
      curve_points.back() = cell.template cut_exit<SIM_C>();
    }

    return curve_points;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto get_cuts() const noexcept -> std::vector<SimCoord<ActiveFloat>> {
    if (m_cut_cell_idxs.empty()) { return {}; }

    std::vector<SimCoord<ActiveFloat>> cut_points(2 * m_cut_cell_idxs.size());

    for (size_t i = 0; i < m_cut_cell_idxs.size(); ++i) {
      const auto& cell      = m_cells[m_cut_cell_idxs[i]];
      cut_points[2 * i + 0] = cell.template cut_entry<SIM_C>();
      cut_points[2 * i + 1] = cell.template cut_exit<SIM_C>();
    }

    return cut_points;
  }

  // -----------------------------------------------------------------------------------------------
  static constexpr std::array<Side, 8UZ> NEIGHBOUR_ORDER = {
      BOTTOM,
      static_cast<Side>(BOTTOM | RIGHT),
      RIGHT,
      static_cast<Side>(TOP | RIGHT),
      TOP,
      static_cast<Side>(TOP | LEFT),
      LEFT,
      static_cast<Side>(BOTTOM | LEFT),
  };

  // Order is bottom, bottom-right, right, top-right, top, top-left, left, bottom-left
  [[nodiscard]] constexpr auto get_neighbours(size_t idx) const noexcept -> SmallVector<size_t> {
    SmallVector<size_t> neighbours(8UZ, NULL_INDEX);
    if (!is_cell(idx)) { return neighbours; }
    const auto& center_cell = m_cells[idx];
    size_t counter          = 0UZ;

    // Bottom cell
    if (is_cell(center_cell.bottom_idx)) { neighbours[counter] = center_cell.bottom_idx; }
    counter += 1;

    // Bottom right diagonal cell
    if (is_cell(center_cell.right_idx) && is_cell(center_cell.bottom_idx)) {
      IGOR_ASSERT(
          m_cells[center_cell.right_idx].bottom_idx == m_cells[center_cell.bottom_idx].right_idx,
          "Diagonal cell is different depending on path, grid was incorrectly constructed.");
      if (is_cell(m_cells[center_cell.right_idx].bottom_idx)) {
        neighbours[counter] = m_cells[center_cell.right_idx].bottom_idx;
      }
    }
    counter += 1;

    // Right cell
    if (is_cell(center_cell.right_idx)) { neighbours[counter] = center_cell.right_idx; }
    counter += 1;

    // Top right diagonal cell
    if (is_cell(center_cell.right_idx) && is_cell(center_cell.top_idx)) {
      IGOR_ASSERT(
          m_cells[center_cell.right_idx].top_idx == m_cells[center_cell.top_idx].right_idx,
          "Diagonal cell is different depending on path, grid was incorrectly constructed.");
      if (is_cell(m_cells[center_cell.right_idx].top_idx)) {
        neighbours[counter] = m_cells[center_cell.right_idx].top_idx;
      }
    }
    counter += 1;

    // Top cell
    if (is_cell(center_cell.top_idx)) { neighbours[counter] = center_cell.top_idx; }
    counter += 1;

    // Top left diagonal cell
    if (is_cell(center_cell.left_idx) && is_cell(center_cell.top_idx)) {
      IGOR_ASSERT(
          m_cells[center_cell.left_idx].top_idx == m_cells[center_cell.top_idx].left_idx,
          "Diagonal cell is different depending on path, grid was incorrectly constructed.");
      if (is_cell(m_cells[center_cell.left_idx].top_idx)) {
        neighbours[counter] = m_cells[center_cell.left_idx].top_idx;
      }
    }
    counter += 1;

    // Left cell
    if (is_cell(center_cell.left_idx)) { neighbours[counter] = center_cell.left_idx; }
    counter += 1;

    // Bottom left diagonal cell
    if (is_cell(center_cell.left_idx) && is_cell(center_cell.bottom_idx)) {
      IGOR_ASSERT(
          m_cells[center_cell.left_idx].bottom_idx == m_cells[center_cell.bottom_idx].left_idx,
          "Diagonal cell is different depending on path, grid was incorrectly constructed.");
      if (is_cell(m_cells[center_cell.left_idx].bottom_idx)) {
        neighbours[counter] = m_cells[center_cell.left_idx].bottom_idx;
      }
    }
    counter += 1;

    IGOR_ASSERT(counter == 8UZ, "Expected array with exactly 8 entries.");
    return neighbours;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto eval(const SimCoord<PassiveFloat>& point) const noexcept
      -> ActiveFloat {
    return eval(to_grid_coord(point));
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto eval(const GridCoord<PassiveFloat>& point) const noexcept
      -> ActiveFloat {
    const auto cell_idx = find_cell(point);
    if (cell_idx == NULL_INDEX) {
      Igor::Panic(
          "Point {} is not in grid [{}, {}]x[{}, {}].", point, m_x_min, m_x_max, m_y_min, m_y_max);
    }

    const auto& cell = m_cells[cell_idx];
    if (cell.is_cartesian()) {
      return cell.get_cartesian().value;
    } else if (cell.is_cut()) {
      const auto point_in_left  = cell.template get_cut_left_polygon<GRID_C>().contains(point);
      const auto point_in_right = cell.template get_cut_right_polygon<GRID_C>().contains(point);

      if (point_in_left) {
        return cell.get_cut().left_value;
      } else if (point_in_right) {
        return cell.get_cut().right_value;
      } else {
#ifdef ZAP_WARN_GRID_EVAL
        if (point_in_left && point_in_right) {
          Igor::Warn("Point {} is in both subcells, take the average.", point);
        } else if (!(point_in_left || point_in_right)) {
          Igor::Warn("Point {} is in neither subcell, probably a rounding error, take the average.",
                     point);
        }
#endif  // ZAP_WARN_GRID_EVAL
        return (cell.get_cut().left_value + cell.get_cut().right_value) / 2;
      }
    } else {
      Igor::Panic("Unknown cell type with variant index {}", cell.cell_type.index());
      std::unreachable();
    }
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto mass() const noexcept -> ActiveFloat {
    ActiveFloat mass = 0.0;
    for (const auto& cell : m_cells) {
      IGOR_ASSERT(cell.is_cartesian() || cell.is_cut(), "Invalid cell type.");
      if (cell.is_cartesian()) {
        mass += cell.get_cartesian().value * cell.template dx<SIM_C>() * cell.template dy<SIM_C>();
      } else {
        const auto left_area  = cell.template get_cut_left_polygon<SIM_C>().area();
        const auto right_area = cell.template get_cut_right_polygon<SIM_C>().area();
        mass += cell.get_cut().left_value * left_area + cell.get_cut().right_value * right_area;
      }
    }
    return mass;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto nx() const noexcept -> size_t { return m_nx; }
  [[nodiscard]] constexpr auto ny() const noexcept -> size_t { return m_ny; }
  [[nodiscard]] constexpr auto size() const noexcept -> size_t { return m_cells.size(); }
  [[nodiscard]] constexpr auto min_delta() const noexcept -> PassiveFloat {
    return std::min(m_dx, m_dy);
  }
  [[nodiscard]] constexpr auto x_min() const noexcept -> PassiveFloat { return m_x_min; }
  [[nodiscard]] constexpr auto x_max() const noexcept -> PassiveFloat { return m_x_max; }
  [[nodiscard]] constexpr auto y_min() const noexcept -> PassiveFloat { return m_y_min; }
  [[nodiscard]] constexpr auto y_max() const noexcept -> PassiveFloat { return m_y_max; }
  [[nodiscard]] constexpr auto scale_x() const noexcept -> PassiveFloat { return 1 / m_dx; }
  [[nodiscard]] constexpr auto scale_y() const noexcept -> PassiveFloat { return 1 / m_dy; }
  [[nodiscard]] constexpr auto cells() noexcept -> std::vector<Cell<ActiveFloat, PassiveFloat>>& {
    return m_cells;
  }
  [[nodiscard]] constexpr auto cells() const noexcept
      -> const std::vector<Cell<ActiveFloat, PassiveFloat>>& {
    return m_cells;
  }
  [[nodiscard]] constexpr auto cut_cell_idxs() const noexcept -> const std::vector<size_t>& {
    return m_cut_cell_idxs;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto operator[](size_t idx) const noexcept
      -> const Cell<ActiveFloat, PassiveFloat>& {
    IGOR_ASSERT(idx < m_nx * m_ny, "Index {} is out of bounds with max {}.", idx, m_nx * m_ny);
    return m_cells[idx];
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto operator[](size_t idx) noexcept -> Cell<ActiveFloat, PassiveFloat>& {
    IGOR_ASSERT(idx < m_nx * m_ny, "Index {} is out of bounds with max {}.", idx, m_nx * m_ny);
    return m_cells[idx];
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto begin() noexcept { return m_cells.begin(); }
  [[nodiscard]] constexpr auto begin() const noexcept { return m_cells.begin(); }
  [[nodiscard]] constexpr auto cbegin() const noexcept { return m_cells.cbegin(); }
  [[nodiscard]] constexpr auto end() noexcept { return m_cells.end(); }
  [[nodiscard]] constexpr auto end() const noexcept { return m_cells.end(); }
  [[nodiscard]] constexpr auto cend() const noexcept { return m_cells.cend(); }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto is_cell(size_t idx) const noexcept -> bool {
    return idx < m_cells.size();
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto on_boundary(size_t idx) const noexcept -> Side {
    const auto [xi, yi] = to_mat_idx(idx);
    // clang-format off
    return static_cast<Side>(
        LEFT   * (xi == 0)        |
        RIGHT  * (xi == m_nx - 1) |
        BOTTOM * (yi == 0)        |
        TOP    * (yi == m_ny - 1)
    );
    // clang-format on
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto translate(Geometry::Polygon<GridCoord<ActiveFloat>> polygon,
                                         Side side) const noexcept
      -> Geometry::Polygon<GridCoord<ActiveFloat>> {
    IGOR_ASSERT(m_nx > 1 && m_ny > 1,
                "This function assumes that both nx and ny are greater than one, but they are {} "
                "and {} respectively",
                m_nx,
                m_ny);

    if (side == 0) { return polygon; }

    PassiveFloat dx = 0;
    PassiveFloat dy = 0;
    if ((side & LEFT) > 0) { dx -= static_cast<PassiveFloat>(m_nx); }
    if ((side & RIGHT) > 0) { dx += static_cast<PassiveFloat>(m_nx); }
    if ((side & BOTTOM) > 0) { dy -= static_cast<PassiveFloat>(m_ny); }
    if ((side & TOP) > 0) { dy += static_cast<PassiveFloat>(m_ny); }

    return Geometry::translate(std::move(polygon), dx, dy);
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void dump_cells(std::ostream& out) const noexcept {
    for (size_t i = 0; i < m_cells.size(); ++i) {
      const auto& cell = m_cells[i];
      out << i << ": " << cell;
      if (i + 1 == m_cells.size()) {
        out << '\n';
      } else {
        out << ",\n";
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  friend struct Cell<ActiveFloat, PassiveFloat>;
};

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_GRID_HPP_
