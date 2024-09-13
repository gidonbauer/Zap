#ifndef ZAP_CELL_BASED_GEOMETRY_HPP_
#define ZAP_CELL_BASED_GEOMETRY_HPP_

#include <algorithm>

#include <Eigen/Dense>

#include "CellBased/Definitions.hpp"
#include "CellBased/SmallVector.hpp"
#include "Igor.hpp"

namespace Zap::CellBased::Geometry {

template <typename Float>
class Polygon {
  using Point = Eigen::Vector<Float, POINT_SIZE>;
  SmallVector<Point> m_points{};

 public:
  // -----------------------------------------------------------------------------------------------
  constexpr Polygon() noexcept = default;

  constexpr Polygon(SmallVector<Point> points) noexcept
      : m_points(std::move(points)) {
    remove_duplicate_points();
    sort_points_counter_clockwise();
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void sort_points_counter_clockwise() noexcept {
    Point center = Point::Zero();
    for (const auto& p : m_points) {
      center += p;
    }
    center /= static_cast<Float>(m_points.size());

    // Sort points in counter clockwise order
    std::sort(
        std::begin(m_points), std::end(m_points), [&center](const Point& p1, const Point& p2) {
          const auto a1 = std::atan2(p1(X) - center(X), p1(Y) - center(Y));
          const auto a2 = std::atan2(p2(X) - center(X), p2(Y) - center(Y));
          return a1 > a2;
        });
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void remove_duplicate_points() noexcept {
    const auto first_duplicate =
        std::unique(std::begin(m_points), std::end(m_points), [](const Point& p1, const Point& p2) {
          return (p1 - p2).norm() < 1e-8;
        });
    m_points.erase(first_duplicate, std::end(m_points));
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto area() const noexcept -> Float {
    Float double_area{0};
    for (size_t i = 0; i < m_points.size(); ++i) {
      const auto& p1 = m_points[i];
      const auto& p2 = m_points[(i + 1) % m_points.size()];
      double_area += p1(X) * p2(Y) - p1(Y) * p2(X);
    }
    return double_area / 2;
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void add_point(Point p) noexcept {
    if (std::find_if(std::cbegin(m_points), std::cend(m_points), [&p](const Point& e) {
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
  [[nodiscard]] constexpr auto operator[](size_t idx) noexcept -> Point& {
    assert(idx < size());
    return m_points[idx];
  }
  [[nodiscard]] constexpr auto operator[](size_t idx) const noexcept -> const Point& {
    assert(idx < size());
    return m_points[idx];
  }
  [[nodiscard]] constexpr auto points() noexcept -> SmallVector<Point>& { return m_points; }
  [[nodiscard]] constexpr auto points() const noexcept -> const SmallVector<Point>& {
    return m_points;
  }
};

// - Calculate intersection polygon using Sutherland-Hodgman algorithm -----------------------------
template <std::floating_point Float>
[[nodiscard]] constexpr auto
intersection(const Polygon<Float>& polygon1,
             const Polygon<Float>& polygon2) noexcept -> Polygon<Float> {
  // Utility function to check if point is inside a polygon edge
  const auto point_on_line = [](const Eigen::Vector<Float, POINT_SIZE>& p,
                                const Eigen::Vector<Float, POINT_SIZE>& a,
                                const Eigen::Vector<Float, POINT_SIZE>& b) -> bool {
    // Check if point p is to the left of line segment ab
    return (b(X) - a(X)) * (p(Y) - a(Y)) >= (b(Y) - a(Y)) * (p(X) - a(X)) - EPS<Float>;
  };

  // Utility function to compute intersection point of line segment ab with cd
  const auto line_intersect =
      [](const Eigen::Vector<Float, POINT_SIZE>& a,
         const Eigen::Vector<Float, POINT_SIZE>& b,
         const Eigen::Vector<Float, POINT_SIZE>& c,
         const Eigen::Vector<Float, POINT_SIZE>& d) -> Eigen::Vector<Float, POINT_SIZE> {
    Float a1 = b(Y) - a(Y);
    Float b1 = a(X) - b(X);
    Float c1 = a1 * a(X) + b1 * a(Y);

    Float a2 = d(Y) - c(Y);
    Float b2 = c(X) - d(X);
    Float c2 = a2 * c(X) + b2 * c(Y);

    Float determinant = a1 * b2 - a2 * b1;
    assert(std::abs(determinant) >= EPS<Float>);

    return {(b2 * c1 - b1 * c2) / determinant, (a1 * c2 - a2 * c1) / determinant};
  };

  if (polygon1.empty() || polygon2.empty()) { return Polygon<Float>{}; }

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
        if (!point_on_line(P, A, B)) {
          output_list.add_point(line_intersect(A, B, P, Q));  //
        }
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
