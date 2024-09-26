#ifndef ZAP_SCALAR_BOUNDARY_CONDITIONS_HPP_
#define ZAP_SCALAR_BOUNDARY_CONDITIONS_HPP_

#include "MatBased/Matrix.hpp"

namespace Zap::MatBased {

// - Copy the value next to the boundary -----------------------------------------------------------
template <typename Float, typename NumericalFluxX, typename NumericalFluxY>
constexpr void copy_boundary_conditions(Matrix<Float>& u_next,
                                        const Matrix<Float>& u_curr,
                                        Float /*dt*/,
                                        Float /*dx*/,
                                        Float /*dy*/,
                                        NumericalFluxX /*numerical_flux_x*/,
                                        NumericalFluxY /*numerical_flux_y*/) noexcept {
  assert(u_curr.rows() == u_next.rows());
  assert(u_curr.cols() == u_next.cols());

  for (int row = 0; row < u_curr.rows(); ++row) {
    u_next(row, 0)                 = u_next(row, 1);
    u_next(row, u_next.cols() - 1) = u_next(row, u_next.cols() - 2);
  }
  for (int col = 0; col < u_curr.cols(); ++col) {
    u_next(0, col)                 = u_next(1, col);
    u_next(u_next.rows() - 1, col) = u_next(u_next.rows() - 2, col);
  }
}

// - Zero flux over the boundaries -----------------------------------------------------------------
template <typename Float, typename NumericalFluxX, typename NumericalFluxY>
constexpr void zero_flux_boundary(Matrix<Float>& u_next,
                                  const Matrix<Float>& u_curr,
                                  Float dt,
                                  Float dx,
                                  Float dy,
                                  NumericalFluxX numerical_flux_x,
                                  NumericalFluxY numerical_flux_y) noexcept {
  assert(u_curr.rows() == u_next.rows());
  assert(u_curr.cols() == u_next.cols());

  // Update bottom left corner
  {
    const auto F_x_plus = numerical_flux_x(u_curr(0, 0), u_curr(0, 1));
    const auto F_y_plus = numerical_flux_y(u_curr(0, 0), u_curr(1, 0));

    u_next(0, 0) = u_curr(0, 0) - (dt / dx) * (F_x_plus /*- F_x_minus*/) -
                   (dt / dy) * (F_y_plus /*- F_y_minus*/);
  }

  // Update bottom right corner
  {
    const auto F_x_minus =
        numerical_flux_x(u_curr(0, u_curr.cols() - 2), u_curr(0, u_curr.cols() - 1));
    const auto F_y_plus =
        numerical_flux_y(u_curr(0, u_curr.cols() - 1), u_curr(1, u_curr.cols() - 1));

    u_next(0, u_curr.cols() - 1) = u_curr(0, u_curr.cols() - 1) -
                                   (dt / dx) * (/*F_x_plus -*/ F_x_minus) -
                                   (dt / dy) * (F_y_plus /*- F_y_minus*/);
  }

  // Update top left corner
  {
    const auto F_x_plus =
        numerical_flux_x(u_curr(u_curr.rows() - 1, 0), u_curr(u_curr.rows() - 1, 1));
    const auto F_y_minus =
        numerical_flux_y(u_curr(u_curr.rows() - 2, 0), u_curr(u_curr.rows() - 1, 0));

    u_next(u_curr.rows() - 1, 0) = u_curr(u_curr.rows() - 1, 0) -
                                   (dt / dx) * (F_x_plus /*- F_x_minus*/) -
                                   (dt / dy) * (/*F_y_plus -*/ F_y_minus);
  }

  // Update top right corner
  {
    const auto F_x_minus = numerical_flux_x(u_curr(u_curr.rows() - 1, u_curr.cols() - 2),
                                            u_curr(u_curr.rows() - 1, u_curr.cols() - 1));
    const auto F_y_minus = numerical_flux_y(u_curr(u_curr.rows() - 2, u_curr.cols() - 1),
                                            u_curr(u_curr.rows() - 1, u_curr.cols() - 1));

    u_next(u_curr.rows() - 1, u_curr.cols() - 1) = u_curr(u_curr.rows() - 1, u_curr.cols() - 1) -
                                                   (dt / dx) * (/*F_x_plus -*/ F_x_minus) -
                                                   (dt / dy) * (/*F_y_plus -*/ F_y_minus);
  }

  // Update left and right side
  for (int row = 1; row < u_curr.rows() - 1; ++row) {
    // Left
    {
      const auto F_x_plus  = numerical_flux_x(u_curr(row, 0), u_curr(row, 1));
      const auto F_y_minus = numerical_flux_y(u_curr(row - 1, 0), u_curr(row, 0));
      const auto F_y_plus  = numerical_flux_y(u_curr(row, 0), u_curr(row + 1, 0));

      u_next(row, 0) = u_curr(row, 0) - (dt / dx) * (F_x_plus /*- F_x_minus*/) -
                       (dt / dy) * (F_y_plus - F_y_minus);
    }

    // Right
    {
      const auto F_x_minus =
          numerical_flux_x(u_curr(row, u_next.cols() - 2), u_curr(row, u_next.cols() - 1));
      const auto F_y_minus =
          numerical_flux_y(u_curr(row - 1, u_next.cols() - 1), u_curr(row, u_next.cols() - 1));
      const auto F_y_plus =
          numerical_flux_y(u_curr(row, u_next.cols() - 1), u_curr(row + 1, u_next.cols() - 1));

      u_next(row, u_next.cols() - 1) = u_curr(row, u_next.cols() - 1) -
                                       (dt / dx) * (/*F_x_plus -*/ F_x_minus) -
                                       (dt / dy) * (F_y_plus - F_y_minus);
    }
  }

  // Update top and bottom side
  for (int col = 1; col < u_curr.cols() - 1; ++col) {
    // Bottom
    {
      const auto F_x_minus = numerical_flux_x(u_curr(0, col - 1), u_curr(0, col));
      const auto F_x_plus  = numerical_flux_x(u_curr(0, col), u_curr(0, col + 1));
      const auto F_y_plus  = numerical_flux_y(u_curr(0, col), u_curr(0 + 1, col));

      u_next(0, col) = u_curr(0, col) - (dt / dx) * (F_x_plus - F_x_minus) -
                       (dt / dy) * (F_y_plus /*- F_y_minus*/);
    }

    // Top
    {
      const auto F_x_minus =
          numerical_flux_x(u_curr(u_next.rows() - 1, col - 1), u_curr(u_next.rows() - 1, col));
      const auto F_x_plus =
          numerical_flux_x(u_curr(u_next.rows() - 1, col), u_curr(u_next.rows() - 1, col + 1));
      const auto F_y_minus =
          numerical_flux_y(u_curr(u_next.rows() - 2, col), u_curr(u_next.rows() - 1, col));

      u_next(u_next.rows() - 1, col) = u_curr(u_next.rows() - 1, col) -
                                       (dt / dx) * (F_x_plus - F_x_minus) -
                                       (dt / dy) * (/*F_y_plus -*/ F_y_minus);
    }
  }
}

// - Rotated periodic boundaries left -> 90° * bottom and bottom -> -90° * left --------------------
template <typename Float, typename NumericalFluxX, typename NumericalFluxY>
constexpr void rotated_periodic_boundary(Matrix<Float>& u_next,
                                         const Matrix<Float>& u_curr,
                                         Float dt,
                                         Float dx,
                                         Float dy,
                                         NumericalFluxX numerical_flux_x,
                                         NumericalFluxY numerical_flux_y) noexcept {
  assert(u_curr.rows() == u_next.rows());
  assert(u_curr.cols() == u_next.cols());
  assert(u_curr.rows() == u_curr.cols());

  // Update bottom left corner
  {
    const auto F_x_minus = numerical_flux_x(u_curr(0, 0), u_curr(0, 0));
    const auto F_y_minus = numerical_flux_y(u_curr(0, 0), u_curr(0, 0));
    const auto F_x_plus  = numerical_flux_x(u_curr(0, 0), u_curr(0, 1));
    const auto F_y_plus  = numerical_flux_y(u_curr(0, 0), u_curr(1, 0));

    u_next(0, 0) =
        u_curr(0, 0) - (dt / dx) * (F_x_plus - F_x_minus) - (dt / dy) * (F_y_plus - F_y_minus);
  }

  // Update bottom right corner
  {
    const auto F_x_minus =
        numerical_flux_x(u_curr(0, u_curr.cols() - 2), u_curr(0, u_curr.cols() - 1));
    const auto F_y_minus =
        numerical_flux_y(u_curr(u_curr.cols() - 1, 0), u_curr(0, u_curr.cols() - 1));
    const auto F_y_plus =
        numerical_flux_y(u_curr(0, u_curr.cols() - 1), u_curr(1, u_curr.cols() - 1));

    // zero flux over right boundary!
    u_next(0, u_curr.cols() - 1) = u_curr(0, u_curr.cols() - 1) -
                                   (dt / dx) * (/*F_x_plus -*/ F_x_minus) -
                                   (dt / dy) * (F_y_plus - F_y_minus);
  }

  // Update top left corner
  {
    const auto F_x_minus =
        numerical_flux_x(u_curr(0, u_curr.rows() - 1), u_curr(u_curr.rows() - 1, 0));
    const auto F_x_plus =
        numerical_flux_x(u_curr(u_curr.rows() - 1, 0), u_curr(u_curr.rows() - 1, 1));
    const auto F_y_minus =
        numerical_flux_y(u_curr(u_curr.rows() - 2, 0), u_curr(u_curr.rows() - 1, 0));

    u_next(u_curr.rows() - 1, 0) = u_curr(u_curr.rows() - 1, 0) -
                                   (dt / dx) * (F_x_plus - F_x_minus) -
                                   (dt / dy) * (/* F_y_plus -*/ F_y_minus);
  }

  // Update top right corner
  {
    const auto F_x_minus = numerical_flux_x(u_curr(u_curr.rows() - 1, u_curr.cols() - 2),
                                            u_curr(u_curr.rows() - 1, u_curr.cols() - 1));
    const auto F_y_minus = numerical_flux_y(u_curr(u_curr.rows() - 2, u_curr.cols() - 1),
                                            u_curr(u_curr.rows() - 1, u_curr.cols() - 1));

    u_next(u_curr.rows() - 1, u_curr.cols() - 1) = u_curr(u_curr.rows() - 1, u_curr.cols() - 1) -
                                                   (dt / dx) * (/*F_x_plus -*/ F_x_minus) -
                                                   (dt / dy) * (/*F_y_plus -*/ F_y_minus);
  }

  // Update left and right side
  for (int row = 1; row < u_curr.rows() - 1; ++row) {
    // Left
    {
      const auto F_x_minus = numerical_flux_x(u_curr(0, row), u_curr(row, 0));
      const auto F_x_plus  = numerical_flux_x(u_curr(row, 0), u_curr(row, 1));
      const auto F_y_minus = numerical_flux_y(u_curr(row - 1, 0), u_curr(row, 0));
      const auto F_y_plus  = numerical_flux_y(u_curr(row, 0), u_curr(row + 1, 0));

      u_next(row, 0) =
          u_curr(row, 0) - (dt / dx) * (F_x_plus - F_x_minus) - (dt / dy) * (F_y_plus - F_y_minus);
    }
    // Right
    {
      const auto F_x_minus =
          numerical_flux_x(u_curr(row, u_next.cols() - 2), u_curr(row, u_next.cols() - 1));
      const auto F_y_minus =
          numerical_flux_y(u_curr(row - 1, u_next.cols() - 1), u_curr(row, u_next.cols() - 1));
      const auto F_y_plus =
          numerical_flux_y(u_curr(row, u_next.cols() - 1), u_curr(row + 1, u_next.cols() - 1));

      u_next(row, u_next.cols() - 1) = u_curr(row, u_next.cols() - 1) -
                                       (dt / dx) * (/*F_x_plus -*/ F_x_minus) -
                                       (dt / dy) * (F_y_plus - F_y_minus);
    }
  }

  // Update top and bottom side
  for (int col = 1; col < u_curr.cols() - 1; ++col) {
    // Top
    {
      const auto F_x_minus =
          numerical_flux_x(u_curr(u_next.rows() - 1, col - 1), u_curr(u_next.rows() - 1, col));
      const auto F_x_plus =
          numerical_flux_x(u_curr(u_next.rows() - 1, col), u_curr(u_next.rows() - 1, col + 1));
      const auto F_y_minus =
          numerical_flux_y(u_curr(u_next.rows() - 2, col), u_curr(u_next.rows() - 1, col));

      u_next(u_next.rows() - 1, col) = u_curr(u_next.rows() - 1, col) -
                                       (dt / dx) * (F_x_plus - F_x_minus) -
                                       (dt / dy) * (/* F_y_plus - */ F_y_minus);
    }
    // Bottom
    {
      const auto F_x_minus = numerical_flux_x(u_curr(0, col - 1), u_curr(0, col));
      const auto F_x_plus  = numerical_flux_x(u_curr(0, col), u_curr(0, col + 1));
      const auto F_y_plus  = numerical_flux_y(u_curr(0, col), u_curr(0 + 1, col));
      const auto F_y_minus = numerical_flux_y(u_curr(col, 0), u_curr(0, col));

      u_next(0, col) =
          u_curr(0, col) - (dt / dx) * (F_x_plus - F_x_minus) - (dt / dy) * (F_y_plus - F_y_minus);
    }
  }
}

// - Regular periodic boundaries -------------------------------------------------------------------
template <typename Float, typename NumericalFluxX, typename NumericalFluxY>
void periodic_boundary(Matrix<Float>& u_next,
                       const Matrix<Float>& u_curr,
                       Float dt,
                       Float dx,
                       Float dy,
                       NumericalFluxX numerical_flux_x,
                       NumericalFluxY numerical_flux_y) noexcept {
  assert(u_curr.rows() == u_next.rows());
  assert(u_curr.cols() == u_next.cols());

  constexpr auto wrapping_idx = [](Eigen::Index idx,
                                   Eigen::Index offset,
                                   Eigen::Index size) constexpr noexcept -> Eigen::Index {
    return (idx + size + offset) % size;
  };

#pragma omp parallel for
  for (Eigen::Index row = 0; row < u_curr.rows(); ++row) {
    // Update left side
    {
      const auto F_x_minus = numerical_flux_x(u_curr(row, u_curr.cols() - 1), u_curr(row, 0));
      const auto F_x_plus  = numerical_flux_x(u_curr(row, 0), u_curr(row, 1));
      const auto F_y_minus =
          numerical_flux_y(u_curr(wrapping_idx(row, -1, u_curr.rows()), 0), u_curr(row, 0));
      const auto F_y_plus =
          numerical_flux_y(u_curr(row, 0), u_curr(wrapping_idx(row, 1, u_curr.rows()), 0));

      u_next(row, 0) =
          u_curr(row, 0) - (dt / dx) * (F_x_plus - F_x_minus) - (dt / dy) * (F_y_plus - F_y_minus);
    }

    // Update right side
    {
      const auto F_x_minus =
          numerical_flux_x(u_curr(row, u_next.cols() - 2), u_curr(row, u_next.cols() - 1));
      const auto F_x_plus =
          numerical_flux_x(u_curr(row, u_curr.cols() - 1), u_curr(0, u_curr.cols() - 1));
      const auto F_y_minus =
          numerical_flux_y(u_curr(wrapping_idx(row, -1, u_curr.rows()), u_next.cols() - 1),
                           u_curr(row, u_next.cols() - 1));
      const auto F_y_plus =
          numerical_flux_y(u_curr(row, u_next.cols() - 1),
                           u_curr(wrapping_idx(row, 1, u_curr.rows()), u_next.cols() - 1));

      u_next(row, u_next.cols() - 1) = u_curr(row, u_next.cols() - 1) -
                                       (dt / dx) * (F_x_plus - F_x_minus) -
                                       (dt / dy) * (F_y_plus - F_y_minus);
    }
  }

#pragma omp parallel for
  for (Eigen::Index col = 1; col < u_curr.cols() - 1; ++col) {
    // Update bottom side
    {
      const auto F_x_minus = numerical_flux_x(u_curr(0, col - 1), u_curr(0, col));
      const auto F_x_plus  = numerical_flux_x(u_curr(0, col), u_curr(0, col + 1));
      const auto F_y_minus = numerical_flux_y(u_curr(u_curr.rows() - 1, col), u_curr(0, col));
      const auto F_y_plus  = numerical_flux_y(u_curr(0, col), u_curr(1, col));

      u_next(0, col) =
          u_curr(0, col) - (dt / dx) * (F_x_plus - F_x_minus) - (dt / dy) * (F_y_plus - F_y_minus);
    }

    // Update top side
    {
      const auto F_x_minus =
          numerical_flux_x(u_curr(u_next.rows() - 1, col - 1), u_curr(u_next.rows() - 1, col));
      const auto F_x_plus =
          numerical_flux_x(u_curr(u_next.rows() - 1, col), u_curr(u_next.rows() - 1, col + 1));
      const auto F_y_minus =
          numerical_flux_y(u_curr(u_next.rows() - 2, col), u_curr(u_next.rows() - 1, col));
      const auto F_y_plus = numerical_flux_y(u_curr(u_curr.rows() - 1, col), u_curr(0, col));

      u_next(u_next.rows() - 1, col) = u_curr(u_next.rows() - 1, col) -
                                       (dt / dx) * (F_x_plus - F_x_minus) -
                                       (dt / dy) * (F_y_plus - F_y_minus);
    }
  }
}

// - Boundary assumes that the value outside is the same as the cell value -------------------------
template <typename Float, typename NumericalFluxX, typename NumericalFluxY>
void equal_value_boundary(Matrix<Float>& u_next,
                          const Matrix<Float>& u_curr,
                          Float dt,
                          Float dx,
                          Float dy,
                          NumericalFluxX numerical_flux_x,
                          NumericalFluxY numerical_flux_y) noexcept {
  assert(u_curr.rows() == u_next.rows());
  assert(u_curr.cols() == u_next.cols());

  // Update bottom left corner
  {
    const auto F_x_minus = numerical_flux_x(u_curr(0, 0), u_curr(0, 0));
    const auto F_x_plus  = numerical_flux_x(u_curr(0, 0), u_curr(0, 1));
    const auto F_y_minus = numerical_flux_y(u_curr(0, 0), u_curr(0, 0));
    const auto F_y_plus  = numerical_flux_y(u_curr(0, 0), u_curr(1, 0));

    u_next(0, 0) =
        u_curr(0, 0) - (dt / dx) * (F_x_plus - F_x_minus) - (dt / dy) * (F_y_plus - F_y_minus);
  }

  // Update bottom right corner
  {
    const auto F_x_minus =
        numerical_flux_x(u_curr(0, u_curr.cols() - 2), u_curr(0, u_curr.cols() - 1));
    const auto F_x_plus =
        numerical_flux_x(u_curr(0, u_curr.cols() - 1), u_curr(0, u_curr.cols() - 1));
    const auto F_y_minus =
        numerical_flux_y(u_curr(0, u_curr.cols() - 1), u_curr(1, u_curr.cols() - 1));
    const auto F_y_plus =
        numerical_flux_y(u_curr(0, u_curr.cols() - 1), u_curr(1, u_curr.cols() - 1));

    u_next(0, u_curr.cols() - 1) = u_curr(0, u_curr.cols() - 1) -
                                   (dt / dx) * (F_x_plus - F_x_minus) -
                                   (dt / dy) * (F_y_plus - F_y_minus);
  }

  // Update top left corner
  {
    const auto F_x_minus =
        numerical_flux_x(u_curr(u_curr.rows() - 1, 0), u_curr(u_curr.rows() - 1, 0));
    const auto F_x_plus =
        numerical_flux_x(u_curr(u_curr.rows() - 1, 0), u_curr(u_curr.rows() - 1, 1));
    const auto F_y_minus =
        numerical_flux_y(u_curr(u_curr.rows() - 2, 0), u_curr(u_curr.rows() - 1, 0));
    const auto F_y_plus =
        numerical_flux_y(u_curr(u_curr.rows() - 1, 0), u_curr(u_curr.rows() - 1, 0));

    u_next(u_curr.rows() - 1, 0) = u_curr(u_curr.rows() - 1, 0) -
                                   (dt / dx) * (F_x_plus - F_x_minus) -
                                   (dt / dy) * (F_y_plus - F_y_minus);
  }

  // Update top right corner
  {
    const auto F_x_minus = numerical_flux_x(u_curr(u_curr.rows() - 1, u_curr.cols() - 2),
                                            u_curr(u_curr.rows() - 1, u_curr.cols() - 1));
    const auto F_x_plus  = numerical_flux_x(u_curr(u_curr.rows() - 1, u_curr.cols() - 1),
                                           u_curr(u_curr.rows() - 1, u_curr.cols() - 1));
    const auto F_y_minus = numerical_flux_y(u_curr(u_curr.rows() - 2, u_curr.cols() - 1),
                                            u_curr(u_curr.rows() - 1, u_curr.cols() - 1));
    const auto F_y_plus  = numerical_flux_y(u_curr(u_curr.rows() - 1, u_curr.cols() - 1),
                                           u_curr(u_curr.rows() - 1, u_curr.cols() - 1));

    u_next(u_curr.rows() - 1, u_curr.cols() - 1) = u_curr(u_curr.rows() - 1, u_curr.cols() - 1) -
                                                   (dt / dx) * (F_x_plus - F_x_minus) -
                                                   (dt / dy) * (F_y_plus - F_y_minus);
  }

  // Update left and right side
  for (int row = 1; row < u_curr.rows() - 1; ++row) {
    // Left
    {
      const auto F_x_minus = numerical_flux_x(u_curr(row, 0), u_curr(row, 0));
      const auto F_x_plus  = numerical_flux_x(u_curr(row, 0), u_curr(row, 1));
      const auto F_y_minus = numerical_flux_y(u_curr(row - 1, 0), u_curr(row, 0));
      const auto F_y_plus  = numerical_flux_y(u_curr(row, 0), u_curr(row + 1, 0));

      u_next(row, 0) =
          u_curr(row, 0) - (dt / dx) * (F_x_plus - F_x_minus) - (dt / dy) * (F_y_plus - F_y_minus);
    }

    // Right
    {
      const auto F_x_minus =
          numerical_flux_x(u_curr(row, u_next.cols() - 2), u_curr(row, u_next.cols() - 1));
      const auto F_x_plus =
          numerical_flux_x(u_curr(row, u_next.cols() - 1), u_curr(row, u_next.cols() - 1));
      const auto F_y_minus =
          numerical_flux_y(u_curr(row - 1, u_next.cols() - 1), u_curr(row, u_next.cols() - 1));
      const auto F_y_plus =
          numerical_flux_y(u_curr(row, u_next.cols() - 1), u_curr(row + 1, u_next.cols() - 1));

      u_next(row, u_next.cols() - 1) = u_curr(row, u_next.cols() - 1) -
                                       (dt / dx) * (F_x_plus - F_x_minus) -
                                       (dt / dy) * (F_y_plus - F_y_minus);
    }
  }

  // Update top and bottom side
  for (int col = 1; col < u_curr.cols() - 1; ++col) {
    // Bottom
    {
      const auto F_x_minus = numerical_flux_x(u_curr(0, col - 1), u_curr(0, col));
      const auto F_x_plus  = numerical_flux_x(u_curr(0, col), u_curr(0, col + 1));
      const auto F_y_minus = numerical_flux_y(u_curr(0, col), u_curr(0, col));
      const auto F_y_plus  = numerical_flux_y(u_curr(0, col), u_curr(1, col));

      u_next(0, col) =
          u_curr(0, col) - (dt / dx) * (F_x_plus - F_x_minus) - (dt / dy) * (F_y_plus - F_y_minus);
    }

    // Top
    {
      const auto F_x_minus =
          numerical_flux_x(u_curr(u_next.rows() - 1, col - 1), u_curr(u_next.rows() - 1, col));
      const auto F_x_plus =
          numerical_flux_x(u_curr(u_next.rows() - 1, col), u_curr(u_next.rows() - 1, col + 1));
      const auto F_y_minus =
          numerical_flux_y(u_curr(u_next.rows() - 2, col), u_curr(u_next.rows() - 1, col));
      const auto F_y_plus =
          numerical_flux_y(u_curr(u_next.rows() - 1, col), u_curr(u_next.rows() - 1, col));

      u_next(u_next.rows() - 1, col) = u_curr(u_next.rows() - 1, col) -
                                       (dt / dx) * (F_x_plus - F_x_minus) -
                                       (dt / dy) * (F_y_plus - F_y_minus);
    }
  }
}

}  // namespace Zap::MatBased

#endif  // ZAP_SCALAR_BOUNDARY_CONDITIONS_HPP_
