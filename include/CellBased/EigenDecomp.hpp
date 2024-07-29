#ifndef ZAP_CELL_BASED_EIGEN_DECOMP_HPP_
#define ZAP_CELL_BASED_EIGEN_DECOMP_HPP_

#include <concepts>

#include <Eigen/Dense>

namespace Zap::CellBased {

// -------------------------------------------------------------------------------------------------
namespace SingleEq {

struct A {
  static constexpr size_t DIM = 1;

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto mat(Float u) const noexcept -> Float {
    return u;
  }

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto eig_vals(Float u) const noexcept -> Float {
    return u;
  }

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto eig_vecs(Float /*u*/) const noexcept -> Float {
    return 1;
  }
};

using B = A;

}  // namespace SingleEq

// -------------------------------------------------------------------------------------------------
namespace DefaultSystem {

struct A {
  static constexpr size_t DIM = 2;

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto
  mat(const Eigen::Vector<Float, DIM>& u) const noexcept -> Eigen::Matrix<Float, DIM, DIM> {
    return Eigen::DiagonalMatrix<Float, DIM, DIM>{u(0), u(0)};
  }

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto
  eig_vals(const Eigen::Vector<Float, DIM>& u) const noexcept -> Eigen::Matrix<Float, DIM, DIM> {
    return Eigen::DiagonalMatrix<Float, DIM, DIM>{u(0), u(0)};
  }

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto eig_vecs(const Eigen::Vector<Float, DIM>& /*u*/) const noexcept
      -> Eigen::Matrix<Float, DIM, DIM> {
    return Eigen::Matrix<Float, DIM, DIM>::Identity();
  }
};

struct B {
  static constexpr size_t DIM = 2;

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto
  mat(const Eigen::Vector<Float, DIM>& u) const noexcept -> Eigen::Matrix<Float, DIM, DIM> {
    return Eigen::DiagonalMatrix<Float, DIM, DIM>{u(1), u(1)};
  }

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto
  eig_vals(const Eigen::Vector<Float, DIM>& u) const noexcept -> Eigen::Matrix<Float, DIM, DIM> {
    return Eigen::DiagonalMatrix<Float, DIM, DIM>{u(1), u(1)};
  }

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto eig_vecs(const Eigen::Vector<Float, DIM>& /*u*/) const noexcept
      -> Eigen::Matrix<Float, DIM, DIM> {
    return Eigen::Matrix<Float, DIM, DIM>::Identity();
  }
};

}  // namespace DefaultSystem

// -------------------------------------------------------------------------------------------------
namespace ExtendedSystem {

class A {
 private:
  static constexpr double m_eps = 1e-6;

 public:
  static constexpr size_t DIM = 2;

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto
  mat(const Eigen::Vector<Float, DIM>& u) const noexcept -> Eigen::Matrix<Float, DIM, DIM> {
    return Eigen::Matrix<Float, DIM, DIM>{2 * u(0), 0, u(1), u(0)};
  }

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto
  eig_vals(const Eigen::Vector<Float, DIM>& u) const noexcept -> Eigen::Matrix<Float, DIM, DIM> {
    return Eigen::DiagonalMatrix<Float, DIM, DIM>{2 * u(0), u(0)};
  }

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto
  eig_vecs(const Eigen::Vector<Float, DIM>& u) const noexcept -> Eigen::Matrix<Float, DIM, DIM> {
    Eigen::Matrix<Float, DIM, DIM> V{};

    // First eigenvector
    if (std::abs(u(1)) >= m_eps) {
      V(0, 0) = u(0) / u(1);
      V(1, 0) = 1;
    } else {
      V(0, 0) = 1;
      V(1, 0) = 0;
    }

    // Second eigenvector
    V(0, 1) = 0;
    V(1, 1) = 1;

    return V;
  }
};

struct B {
 private:
  static constexpr double m_eps = 1e-6;

 public:
  static constexpr size_t DIM = 2;

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto
  mat(const Eigen::Vector<Float, DIM>& u) const noexcept -> Eigen::Matrix<Float, DIM, DIM> {
    return Eigen::Matrix<Float, DIM, DIM>{u(1), u(0), 0, 2 * u(1)};
  }

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto
  eig_vals(const Eigen::Vector<Float, DIM>& u) const noexcept -> Eigen::Matrix<Float, DIM, DIM> {
    return Eigen::DiagonalMatrix<Float, DIM, DIM>{2 * u(1), u(1)};
  }

  template <std::floating_point Float>
  [[nodiscard]] constexpr auto
  eig_vecs(const Eigen::Vector<Float, DIM>& u) const noexcept -> Eigen::Matrix<Float, DIM, DIM> {
    Eigen::Matrix<Float, DIM, DIM> V{};

    // First eigenvector
    if (std::abs(u(1)) >= m_eps) {
      V(0, 0) = u(0) / u(1);
      V(1, 0) = 1;
    } else {
      V(0, 0) = 0;
      V(1, 0) = 1;
    }

    // Second eigenvector
    V(0, 1) = 1;
    V(1, 1) = 0;

    return V;
  }
};

}  // namespace ExtendedSystem

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_EIGEN_DECOMP_HPP_
