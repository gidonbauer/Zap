#ifndef ZAP_CELL_BASED_RECONSTRUCT_SHOCK_HPP_
#define ZAP_CELL_BASED_RECONSTRUCT_SHOCK_HPP_

#include "CellBased/Definitions.hpp"

namespace Zap::CellBased {

// TODO: Think about spline reconstruction.

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, typename PassiveFloat>
[[nodiscard]] constexpr auto
piecewise_linear(PassiveFloat t, const std::vector<SimCoord<ActiveFloat>>& points) noexcept
    -> SimCoord<ActiveFloat> {
  assert(t >= 0 && t <= 1);
  const auto segment = t * static_cast<PassiveFloat>(points.size() - 1);
  auto first_idx     = static_cast<size_t>(std::floor(segment));
  auto last_idx      = static_cast<size_t>(std::floor(segment + 1));
  if (last_idx == points.size()) {
    first_idx -= 1;
    last_idx -= 1;
  }

  assert(first_idx < points.size());
  assert(last_idx < points.size());
  const auto t_prime = segment - static_cast<PassiveFloat>(first_idx);
  assert(t_prime >= 0 && t_prime <= 1);

  return points[first_idx] + t_prime * (points[last_idx] - points[first_idx]);
}

template <typename ActiveFloat, typename PassiveFloat>
[[nodiscard]] constexpr auto
piecewise_linear(const std::vector<PassiveFloat>& ts,
                 const std::vector<SimCoord<ActiveFloat>>& points) noexcept
    -> std::vector<SimCoord<ActiveFloat>> {
  std::vector<SimCoord<ActiveFloat>> curve_points(ts.size());
  std::transform(std::cbegin(ts),
                 std::cend(ts),
                 std::begin(curve_points),
                 [&points](PassiveFloat t) { return piecewise_linear(t, points); });
  return curve_points;
}

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat, typename PassiveFloat>
[[nodiscard]] constexpr auto smoothstep(PassiveFloat t,
                                        const std::vector<SimCoord<ActiveFloat>>& points) noexcept
    -> SimCoord<ActiveFloat> {
  assert(t >= 0 && t <= 1);
  const auto segment = t * static_cast<PassiveFloat>(points.size() - 1);
  auto first_idx     = static_cast<size_t>(std::floor(segment));
  auto last_idx      = static_cast<size_t>(std::floor(segment + 1));
  if (last_idx == points.size()) {
    first_idx -= 1;
    last_idx -= 1;
  }

  assert(first_idx < points.size());
  assert(last_idx < points.size());
  auto t_prime = segment - static_cast<PassiveFloat>(first_idx);
  t_prime      = 3 * t_prime * t_prime - 2 * t_prime * t_prime * t_prime;

  return points[first_idx] + t_prime * (points[last_idx] - points[first_idx]);
}

template <typename ActiveFloat, typename PassiveFloat>
[[nodiscard]] constexpr auto smoothstep(const std::vector<PassiveFloat>& ts,
                                        const std::vector<SimCoord<ActiveFloat>>& points) noexcept
    -> std::vector<SimCoord<ActiveFloat>> {
  std::vector<SimCoord<ActiveFloat>> curve_points(ts.size());
  std::transform(std::cbegin(ts),
                 std::cend(ts),
                 std::begin(curve_points),
                 [&points](PassiveFloat t) { return smoothstep(t, points); });
  return curve_points;
}

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_RECONSTRUCT_SHOCK_HPP_
