#ifndef ZAP_MATRIX_HPP_
#define ZAP_MATRIX_HPP_

#include <Eigen/Core>

namespace Zap::MatBased {

template <typename T>
using Vector = Eigen::Vector<T, Eigen::Dynamic>;
template <typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
constexpr void swap(Matrix<T>& a, Matrix<T>& b) noexcept {
  a.swap(b);
}

}  // namespace Zap::MatBased

#endif  // ZAP_MATRIX_HPP_
