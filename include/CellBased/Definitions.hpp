#ifndef ZAP_CELL_BASED_DEFINITIONS_HPP_
#define ZAP_CELL_BASED_DEFINITIONS_HPP_

#include <concepts>
#include <cstddef>

#include <Eigen/Core>

namespace Zap::CellBased {

enum : size_t { X, Y, POINT_SIZE };

template <typename T>
using Point = Eigen::Vector<T, POINT_SIZE>;

enum Side : int {
  BOTTOM = 0b0001,
  RIGHT  = 0b0010,
  TOP    = 0b0100,
  LEFT   = 0b1000,
  ALL    = LEFT | RIGHT | BOTTOM | TOP,
};

template <std::floating_point Float>
inline constexpr Float EPS;
template <>
inline constexpr float EPS<float> = 1e-6f;
template <>
inline constexpr double EPS<double> = 1e-8;

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_DEFINITIONS_HPP_
