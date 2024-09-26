#ifndef ZAP_SCALAR_SOLVER_HPP_
#define ZAP_SCALAR_SOLVER_HPP_

#include <optional>

#include <Eigen/Dense>

#include "Igor.hpp"
#include "MatBased/Matrix.hpp"

namespace Zap::MatBased {

template <typename Float>
[[nodiscard]] constexpr auto abs_max(const Matrix<Float>& mat) noexcept -> Float {
  assert(mat.rows() > 0 && mat.cols() > 0);
  auto abs_max = std::abs(mat(0, 0));
  for (int col = 0; col < mat.cols(); ++col) {
    for (int row = 0; row < mat.rows(); ++row) {
      abs_max = std::max(abs_max, std::abs(mat(row, col)));
    }
  }
  return abs_max;
}

// - Solve 2D Burgers equation ---------------------------------------------------------------------
//  ^
//  | (m, 0)
//  |   .
//  y (2,0) (2,1)
//  | (1,0) (1,1) (1,2)
//  | (0,0) (0,1) (0,2) ... (0,n)
//  .------ x ---->
template <typename Float, typename BoundaryCondition, typename UWriter, typename TWriter>
[[nodiscard]] auto solve_2d_burgers(const Vector<Float>& x,      // x-grid
                                    const Vector<Float>& y,      // y-grid
                                    const Matrix<Float>& u0,     // initial value
                                    Float tend,                  // Final time
                                    BoundaryCondition boundary,  // Type of boundaries
                                    UWriter& u_writer,           // Save intermediate u in file
                                    TWriter& t_writer            // Save intermediate t in file
                                    ) -> std::optional<Matrix<Float>> {
  assert(y.size() == u0.rows());
  assert(x.size() == u0.cols());
  assert(x.size() >= 2);
  assert(y.size() >= 2);

  // -----------------------------------------------------------------------------------------------

  constexpr auto f = [](const auto& u) constexpr noexcept {
    return static_cast<Float>(0.5) * u * u;
  };

  constexpr auto g = [](const auto& u) constexpr noexcept {
    return static_cast<Float>(0.5) * u * u;
  };

  // From LeVeque: Numerical Methods for Conservation Laws 2nd edition (13.24)
  constexpr auto godunov_flux_x = [f]<typename T>(const T& u_left,
                                                  const T& u_right) constexpr noexcept -> T {
    constexpr auto zero = static_cast<T>(0);
    if (u_left <= u_right) {
      // min u in [u_left, u_right] f(u) = 0.5 * u^2
      if (u_left <= zero && u_right >= zero) { return f(zero); }
      return f(std::min(std::abs(u_left), std::abs(u_right)));
    }
    // max u in [u_right, u_left] f(u) = 0.5 * u^2
    return f(std::max(std::abs(u_left), std::abs(u_right)));
  };

  // From LeVeque: Numerical Methods for Conservation Laws 2nd edition (13.24)
  constexpr auto godunov_flux_y = [g]<typename T>(const T& u_left,
                                                  const T& u_right) constexpr noexcept -> T {
    constexpr auto zero = static_cast<T>(0);
    if (u_left <= u_right) {
      // min u in [u_left, u_right] f(u) = 0.5 * u^2
      if (u_left <= zero && u_right >= zero) { return g(zero); }
      return g(std::min(std::abs(u_left), std::abs(u_right)));
    }
    // max u in [u_right, u_left] f(u) = 0.5 * u^2
    return g(std::max(std::abs(u_left), std::abs(u_right)));
  };

  // -----------------------------------------------------------------------------------------------

  Matrix<Float> u_curr = u0;
  Matrix<Float> u_next = Matrix<Float>::Zero(u0.rows(), u0.cols());

  if (!u_writer.write_data(u_curr) || !t_writer.write_data(static_cast<Float>(0))) {
    return std::nullopt;
  }

  // Assume uniform grid
  const auto dx = x[1] - x[0];
  const auto dy = y[1] - y[0];

  for (Float t = 0.0; t < tend;) {
    const Float CFL_factor = abs_max(u_curr);
    if (std::isnan(CFL_factor) || std::isinf(CFL_factor)) {
      Igor::Warn("CFL_factor is invalid at time t={}: CFL_factor = {}", t, CFL_factor);
      return std::nullopt;
    }
    const Float dt = std::min(0.5 * std::min(dx, dy) / CFL_factor, tend - t);

    // Solve for interior points
#pragma omp parallel for
    for (Eigen::Index row = 1; row < u0.rows() - 1; ++row) {
      for (Eigen::Index col = 1; col < u0.cols() - 1; ++col) {
        // clang-format off
        const auto F_x_minus = godunov_flux_x(u_curr(row,     col - 1), u_curr(row,     col    ));
        const auto F_x_plus  = godunov_flux_x(u_curr(row,     col    ), u_curr(row,     col + 1));
        const auto F_y_minus = godunov_flux_y(u_curr(row - 1, col    ), u_curr(row,     col    ));
        const auto F_y_plus  = godunov_flux_y(u_curr(row,     col    ), u_curr(row + 1, col    ));

        u_next(row, col) = u_curr(row, col) -
                           (dt / dx) * (F_x_plus - F_x_minus) -
                           (dt / dy) * (F_y_plus - F_y_minus);
        // clang-format on
      }
    }

    // Solve for boundary cells
    boundary(u_next, u_curr, dt, dx, dy, godunov_flux_x, godunov_flux_y);

    // Update time
    t += dt;

    Zap::MatBased::swap(u_curr, u_next);

    // - Save intermediate results to file to not exceed memory ------------------------------------
    if (!u_writer.write_data(u_curr) || !t_writer.write_data(t)) { return std::nullopt; }
  }

  return u_curr;
}

}  // namespace Zap::MatBased

#endif  // ZAP_SCALAR_SOLVER_HPP_
