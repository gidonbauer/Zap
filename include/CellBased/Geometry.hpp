#ifndef ZAP_CELL_BASED_GEOMETRY_HPP_
#define ZAP_CELL_BASED_GEOMETRY_HPP_

#include <algorithm>

#include <Eigen/Dense>

#include <AD/ad.hpp>

#include "CellBased/Definitions.hpp"
#include "CellBased/SmallVector.hpp"

namespace Zap::CellBased::Geometry {

template <Point2D_c PointType>
class Polygon {
  static_assert(
      std::is_same_v<decltype(std::declval<PointType>().x), decltype(std::declval<PointType>().y)>,
      "Expect x- and y-component of PointType to have the same type.");
  using ActiveFloat  = decltype(std::declval<PointType>().x);
  using PassiveFloat = std::conditional_t<ad::mode<ActiveFloat>::is_ad_type,
                                          typename ad::mode<ActiveFloat>::passive_t,
                                          ActiveFloat>;

  SmallVector<PointType> m_points{};

 public:
  // -----------------------------------------------------------------------------------------------
  constexpr Polygon() noexcept = default;

  constexpr Polygon(SmallVector<PointType> points) noexcept
      : m_points(std::move(points)) {
    remove_duplicate_points();
    sort_points_counter_clockwise();
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void sort_points_counter_clockwise() noexcept {
    PointType center = PointType::Zero();
    for (const auto& p : m_points) {
      center += p;
    }
    center /= static_cast<PassiveFloat>(m_points.size());

    // Sort points in counter clockwise order
    std::sort(std::begin(m_points),
              std::end(m_points),
              [&center](const PointType& p1, const PointType& p2) {
                const auto a1 = std::atan2(p1.x - center.x, p1.y - center.y);
                const auto a2 = std::atan2(p2.x - center.x, p2.y - center.y);
                return a1 > a2;
              });
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void remove_duplicate_points() noexcept {
    m_points.erase(std::unique(std::begin(m_points),
                               std::end(m_points),
                               [](const PointType& p1, const PointType& p2) {
                                 return (p1 - p2).norm() < 50 * EPS<PassiveFloat>;
                               }),
                   std::end(m_points));
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto area() const noexcept -> ActiveFloat {
    ActiveFloat double_area{0};
    for (size_t i = 0; i < m_points.size(); ++i) {
      const auto& p1 = m_points[i];
      const auto& p2 = m_points[(i + 1) % m_points.size()];
      double_area += p1.x * p2.y - p1.y * p2.x;
    }
    return double_area / 2;
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void add_point(PointType p) noexcept {
    if (std::find_if(std::cbegin(m_points), std::cend(m_points), [&p](const PointType& e) {
          return (p - e).norm() <= 1e-8;
        }) == std::cend(m_points)) {
      m_points.push_back(std::move(p));
      sort_points_counter_clockwise();
    }
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void clear() noexcept { return m_points.clear(); }
  [[nodiscard]] constexpr auto empty() const noexcept -> bool { return m_points.empty(); }
  [[nodiscard]] constexpr auto size() const noexcept -> size_t { return m_points.size(); }
  [[nodiscard]] constexpr auto operator[](size_t idx) noexcept -> PointType& {
    assert(idx < size());
    return m_points[idx];
  }
  [[nodiscard]] constexpr auto operator[](size_t idx) const noexcept -> const PointType& {
    assert(idx < size());
    return m_points[idx];
  }
  [[nodiscard]] constexpr auto points() noexcept -> SmallVector<PointType>& { return m_points; }
  [[nodiscard]] constexpr auto points() const noexcept -> const SmallVector<PointType>& {
    return m_points;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto point_in_polygon(const PointType& p) const noexcept -> bool {
    if (size() < 2) { return false; }

    for (size_t i = 0; i < size(); ++i) {
      const auto& p1 = m_points[i];
      const auto& p2 = m_points[(i + 1) % size()];
      const auto n   = PointType{p1.y - p2.y, -(p1.x - p2.x)};
      const auto v   = n.x * (p.x - p1.x) + n.y * (p.y - p1.y);
      if (v < 0) { return false; }
    }
    return true;
  }
};

// - Calculate intersection polygon using Sutherland-Hodgman algorithm -----------------------------
template <Point2D_c PointType>
[[nodiscard]] constexpr auto intersection(const Polygon<PointType>& polygon1,
                                          const Polygon<PointType>& polygon2) noexcept
    -> Polygon<PointType> {
  static_assert(
      std::is_same_v<decltype(std::declval<PointType>().x), decltype(std::declval<PointType>().y)>,
      "Expect x- and y-component of PointType to have the same type.");
  using Float = decltype(std::declval<PointType>().x);

  // Utility function to check if point is inside a polygon edge
  const auto point_on_line =
      [](const PointType& p, const PointType& a, const PointType& b) -> bool {
    // Check if point p is to the left of line segment ab
    return (b.x - a.x) * (p.y - a.y) >= (b.y - a.y) * (p.x - a.x) - EPS<Float>;
  };

  // Utility function to compute intersection point of line segment ab with cd
  const auto line_intersect = [](const PointType& a,
                                 const PointType& b,
                                 const PointType& c,
                                 const PointType& d) -> PointType {
    Float a1 = b.y - a.y;
    Float b1 = a.x - b.x;
    Float c1 = a1 * a.x + b1 * a.y;

    Float a2 = d.y - c.y;
    Float b2 = c.x - d.x;
    Float c2 = a2 * c.x + b2 * c.y;

    Float determinant = a1 * b2 - a2 * b1;
    IGOR_ASSERT(std::abs(determinant) >= EPS<Float>,
                "Cannot find intersection of lines {}->{} and {}->{}: determinant is {}.",
                a,
                b,
                c,
                d,
                determinant);

    return {(b2 * c1 - b1 * c2) / determinant, (a1 * c2 - a2 * c1) / determinant};
  };

  if (polygon1.empty() || polygon2.empty()) { return Polygon<PointType>{}; }

  Polygon output_list = polygon1;

  for (size_t i = 0; i < polygon2.size(); i++) {
    Polygon input_list = output_list;
    output_list.clear();

    const auto A = polygon2[i];
    const auto B = polygon2[(i + 1) % polygon2.size()];

    for (size_t j = 0; j < input_list.size(); j++) {
      const auto P = input_list[j];
      const auto Q = input_list[(j + 1) % input_list.size()];

      if (point_on_line(Q, A, B)) {
        if (!point_on_line(P, A, B)) { output_list.add_point(line_intersect(A, B, P, Q)); }
        output_list.add_point(Q);
      } else if (point_on_line(P, A, B)) {
        output_list.add_point(line_intersect(A, B, P, Q));
      }
    }
  }

  return output_list;
}

}  // namespace Zap::CellBased::Geometry

#endif  // ZAP_CELL_BASED_GEOMETRY_HPP_
