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
          return (p - e).norm() <= EPS<PassiveFloat>;
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
  [[nodiscard]] constexpr auto contains(const PointType& p) const noexcept -> bool {
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

// - Utility function to compute intersection point of line segment ab with cd ---------------------
template <Point2D_c PointType>
[[nodiscard]] constexpr auto
line_intersect(const PointType& A, const PointType& B, const PointType& C, const PointType& D)
    -> std::optional<PointType> {
  static_assert(
      std::is_same_v<decltype(std::declval<PointType>().x), decltype(std::declval<PointType>().y)>,
      "Expect x- and y-component of PointType to have the same type.");
  using Float = decltype(std::declval<PointType>().x);

  const Float det = (B.x - A.x) * (C.y - D.y) - (B.y - A.y) * (C.x - D.x);
  if (std::abs(det) < EPS<Float>) { return std::nullopt; }

  const Float r = ((C.y - D.y) * (C.x - A.x) + (D.x - C.x) * (C.y - A.y)) / det;
  const Float s = ((A.y - B.y) * (C.x - A.x) + (B.x - A.x) * (C.y - A.y)) / det;

  if (!(0 - EPS<Float> <= r && r <= 1 + EPS<Float>) ||
      !(0 - EPS<Float> <= s && s <= 1 + EPS<Float>)) {
    return std::nullopt;
  }

  return std::optional<PointType>({A.x + r * (B.x - A.x), A.y + r * (B.y - A.y)});
}

// - Calculate intersection of two convex polygons -------------------------------------------------
// Taken from
// https://codereview.stackexchange.com/questions/243743/the-intersection-of-two-polygons-in-c
template <Point2D_c PointType>
[[nodiscard]] constexpr auto operator&(const Polygon<PointType>& polygon1,
                                       const Polygon<PointType>& polygon2) noexcept
    -> Polygon<PointType> {
  // TODO: Consider looking into the Sutherland-Hodgman algorithm
  if (polygon1.empty() || polygon2.empty()) { return Polygon<PointType>{}; }

  Polygon<PointType> intersect;
  for (size_t i = 0; i < polygon1.size(); ++i) {
    // Add points of p1 contained in p2.
    if (polygon2.contains(polygon1[i])) { intersect.add_point(polygon1[i]); }
    for (size_t j = 0; j < polygon2.size(); ++j) {
      // Add points of intersection.
      auto cross = line_intersect(polygon1[i],
                                  polygon1[(i + 1) % polygon1.size()],
                                  polygon2[j],
                                  polygon2[(j + 1) % polygon2.size()]);
      if (cross.has_value()) { intersect.add_point(std::move(*cross)); }
    }
  }

  // Add points of p2 contained in p1.
  for (size_t i = 0; i < polygon2.size(); ++i) {
    if (polygon1.contains(polygon2[i])) { intersect.add_point(polygon2[i]); }
  }

  return intersect;
}

// -------------------------------------------------------------------------------------------------
template <Point2D_c PointType, typename PassiveFloat>
[[nodiscard]] constexpr auto translate(Polygon<PointType> polygon,
                                       PassiveFloat dx,
                                       PassiveFloat dy) noexcept -> Polygon<PointType> {
  static_assert(
      std::is_same_v<decltype(std::declval<PointType>().x), decltype(std::declval<PointType>().y)>,
      "Expect x- and y-component of PointType to have the same type.");
  using PT_ActiveFloat  = decltype(std::declval<PointType>().x);
  using PT_PassiveFloat = std::conditional_t<ad::mode<PT_ActiveFloat>::is_ad_type,
                                             typename ad::mode<PT_ActiveFloat>::passive_t,
                                             PT_ActiveFloat>;
  static_assert(std::is_same_v<PassiveFloat, PT_PassiveFloat>,
                "Expect PassiveFloat and PT_PassiveFloat (from PointType) to be the same type.");

  for (auto& p : polygon.points()) {
    p.x += dx;
    p.y += dy;
  }
  return polygon;
}

}  // namespace Zap::CellBased::Geometry

#endif  // ZAP_CELL_BASED_GEOMETRY_HPP_
