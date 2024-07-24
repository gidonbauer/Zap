#ifndef ZAP_CELL_BASED_EIGEN_DECOMP_HPP_
#define ZAP_CELL_BASED_EIGEN_DECOMP_HPP_

#include <type_traits>

#include <Eigen/Dense>

namespace Zap::CellBased {

// -------------------------------------------------------------------------------------------------
namespace SingleEq {
static constexpr auto eig_vals_x = []<typename Float>(const Float& u) constexpr noexcept {
  static_assert(std::is_floating_point_v<Float>);
  return u;
};
static constexpr auto eig_vecs_x = []<typename Float>(const Float& /*u*/) constexpr noexcept {
  static_assert(std::is_floating_point_v<Float>);
  return Float{1};
};

static constexpr auto eig_vals_y = eig_vals_x;
static constexpr auto eig_vecs_y = eig_vecs_x;
}  // namespace SingleEq

// -------------------------------------------------------------------------------------------------
namespace DefaultSystem {
static constexpr auto eig_vals_x =
    []<typename Float, int DIM>(
        const Eigen::Vector<Float, DIM>& u) noexcept -> Eigen::Matrix<Float, DIM, DIM> {
  static_assert(std::is_floating_point_v<Float>);
  static_assert(DIM == 2);
  return Eigen::DiagonalMatrix<Float, DIM, DIM>{u(0), u(0)};
};
static constexpr auto eig_vecs_x =
    []<typename Float, int DIM>(
        const Eigen::Vector<Float, DIM>& /*u*/) noexcept -> Eigen::Matrix<Float, DIM, DIM> {
  static_assert(std::is_floating_point_v<Float>);
  static_assert(DIM == 2);
  return Eigen::Matrix<Float, DIM, DIM>::Identity();
};

static constexpr auto eig_vals_y =
    []<typename Float, int DIM>(
        const Eigen::Vector<Float, DIM>& u) noexcept -> Eigen::Matrix<Float, DIM, DIM> {
  static_assert(std::is_floating_point_v<Float>);
  static_assert(DIM == 2);
  return Eigen::DiagonalMatrix<Float, DIM, DIM>{u(1), u(1)};
};
static constexpr auto eig_vecs_y =
    []<typename Float, int DIM>(
        const Eigen::Vector<Float, DIM>& /*u*/) noexcept -> Eigen::Matrix<Float, DIM, DIM> {
  static_assert(std::is_floating_point_v<Float>);
  static_assert(DIM == 2);
  return Eigen::Matrix<Float, DIM, DIM>::Identity();
};
}  // namespace DefaultSystem

// -------------------------------------------------------------------------------------------------
namespace ExtendedSystem {
static constexpr auto eig_vals_x =
    []<typename Float, int DIM>(
        const Eigen::Vector<Float, DIM>& u) noexcept -> Eigen::Matrix<Float, DIM, DIM> {
  static_assert(std::is_floating_point_v<Float>);
  static_assert(DIM == 2);
  return Eigen::DiagonalMatrix<Float, DIM, DIM>{2 * u(0), u(0)};
};

static constexpr auto eig_vecs_x =
    []<typename Float, int DIM>(
        const Eigen::Vector<Float, DIM>& u) noexcept -> Eigen::Matrix<Float, DIM, DIM> {
  static_assert(std::is_floating_point_v<Float>);
  static_assert(DIM == 2);
  Eigen::Matrix<Float, DIM, DIM> V{};

  // First eigenvector
  if (std::abs(u(1)) >= 1e-6) {
    V(0, 0) = u(0) / u(1);
    V(1, 0) = 1;
  } else {
    V(0, 0) = 1;
    V(1, 0) = 0;
  }
  // V(0, 0) /= std::sqrt(V(0, 0) * V(0, 0) + V(1, 0) * V(1, 0));
  // V(1, 0) /= std::sqrt(V(0, 0) * V(0, 0) + V(1, 0) * V(1, 0));

  // Second eigenvector
  V(0, 1) = 0;
  V(1, 1) = 1;

  return V;
};

static constexpr auto eig_vals_y =
    []<typename Float, int DIM>(
        const Eigen::Vector<Float, DIM>& u) noexcept -> Eigen::Matrix<Float, DIM, DIM> {
  static_assert(std::is_floating_point_v<Float>);
  static_assert(DIM == 2);
  return Eigen::DiagonalMatrix<Float, DIM, DIM>{2 * u(1), u(1)};
};

static constexpr auto eig_vecs_y =
    []<typename Float, int DIM>(
        const Eigen::Vector<Float, DIM>& u) noexcept -> Eigen::Matrix<Float, DIM, DIM> {
  static_assert(std::is_floating_point_v<Float>);
  static_assert(DIM == 2);
  Eigen::Matrix<Float, DIM, DIM> V{};

  // First eigenvector
  if (std::abs(u(1)) >= 1e-6) {
    V(0, 0) = u(0) / u(1);
    V(1, 0) = 1;
  } else {
    V(0, 0) = 0;
    V(1, 0) = 1;
  }
  // V(0, 0) /= std::sqrt(V(0, 0) * V(0, 0) + V(1, 0) * V(1, 0));
  // V(1, 0) /= std::sqrt(V(0, 0) * V(0, 0) + V(1, 0) * V(1, 0));

  // Second eigenvector
  V(0, 1) = 1;
  V(1, 1) = 0;

  return V;
};
}  // namespace ExtendedSystem

// namespace ExtendedSystem {
// static constexpr auto eig_vals_x =
//     []<typename Float, int DIM>(
//         const Eigen::Vector<Float, DIM>& u) noexcept -> Eigen::Matrix<Float, DIM, DIM> {
//   static_assert(std::is_floating_point_v<Float>);
//   static_assert(DIM == 2);
//   return Eigen::DiagonalMatrix<Float, DIM, DIM>{
//       (3 * u(0) + std::abs(u(0))) / 2,
//       (3 * u(0) - std::abs(u(0))) / 2,
//   };
// };
//
// static constexpr auto eig_vecs_x =
//     []<typename Float, int DIM>(
//         const Eigen::Vector<Float, DIM>& u) noexcept -> Eigen::Matrix<Float, DIM, DIM> {
//   static_assert(std::is_floating_point_v<Float>);
//   static_assert(DIM == 2);
//   Eigen::Matrix<Float, DIM, DIM> V{};
//
//   // First eigenvector
//   if (u(0) >= 0 && std::abs(u(1)) >= 1e-6) {
//     V(0, 0) = u(0) / u(1);
//     V(1, 0) = 1;
//   } else if (u(0) >= 0) {
//     V(0, 0) = 1;
//     V(1, 0) = 0;
//   } else {
//     V(0, 0) = 0;
//     V(1, 0) = 1;
//   }
//   V(0, 0) /= std::sqrt(V(0, 0) * V(0, 0) + V(1, 0) * V(1, 0));
//   V(1, 0) /= std::sqrt(V(0, 0) * V(0, 0) + V(1, 0) * V(1, 0));
//
//   // Second eigenvector
//   if (u(0) < 0 && std::abs(u(1)) >= 1e-6) {
//     V(0, 1) = u(0) / u(1);
//     V(1, 1) = 1;
//   } else if (u(0) < 0) {
//     V(0, 1) = 1;
//     V(1, 1) = 0;
//   } else {
//     V(0, 1) = 0;
//     V(1, 1) = 1;
//   }
//   V(0, 1) /= std::sqrt(V(0, 1) * V(0, 1) + V(1, 1) * V(1, 1));
//   V(1, 1) /= std::sqrt(V(0, 1) * V(0, 1) + V(1, 1) * V(1, 1));
//
//   return V;
// };
//
// static constexpr auto eig_vals_y =
//     []<typename Float, int DIM>(
//         const Eigen::Vector<Float, DIM>& u) noexcept -> Eigen::Matrix<Float, DIM, DIM> {
//   static_assert(std::is_floating_point_v<Float>);
//   static_assert(DIM == 2);
//   return Eigen::DiagonalMatrix<Float, DIM, DIM>{
//       (3 * u(1) + std::abs(u(1))) / 2,
//       (3 * u(1) - std::abs(u(1))) / 2,
//   };
// };
//
// static constexpr auto eig_vecs_y =
//     []<typename Float, int DIM>(
//         const Eigen::Vector<Float, DIM>& u) noexcept -> Eigen::Matrix<Float, DIM, DIM> {
//   static_assert(std::is_floating_point_v<Float>);
//   static_assert(DIM == 2);
//   Eigen::Matrix<Float, DIM, DIM> V{};
//
//   // First eigenvector
//   if (u(1) > 0) {
//     V(0, 0) = u(0) / u(1);
//     V(1, 0) = 1;
//   } else {
//     V(0, 0) = 1;
//     V(1, 0) = 0;
//   }
//   V(0, 0) /= std::sqrt(V(0, 0) * V(0, 0) + V(1, 0) * V(1, 0));
//   V(1, 0) /= std::sqrt(V(0, 0) * V(0, 0) + V(1, 0) * V(1, 0));
//
//   // Second eigenvector
//   if (u(1) >= 0) {
//     V(0, 1) = 1;
//     V(1, 1) = 0;
//   } else {
//     V(0, 1) = u(0) / u(1);
//     V(1, 1) = 1;
//   }
//   V(0, 1) /= std::sqrt(V(0, 1) * V(0, 1) + V(1, 1) * V(1, 1));
//   V(1, 1) /= std::sqrt(V(0, 1) * V(0, 1) + V(1, 1) * V(1, 1));
//
//   return V;
// };
// }  // namespace ExtendedSystem

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_EIGEN_DECOMP_HPP_
