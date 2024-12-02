#ifndef ZAP_SHOCK_AD_1D_SOLVER_HPP_
#define ZAP_SHOCK_AD_1D_SOLVER_HPP_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>

#include <AD/ad.hpp>

namespace Zap::ShockAD_1D {

// -------------------------------------------------------------------------------------------------
// Evaluate piecewise constant reconstruction of solution u
template <typename AT, typename PT>
[[nodiscard]] constexpr auto
piecewise_constant(const std::vector<PT>& x,  // x grid
                   const std::vector<AT>& u,  // solution at grid points at time t
                   const PT& x_eval           // x position at which u is evaluated
                   ) noexcept -> PT {
  assert(x.size() == u.size());
  assert(x.size() >= 2UL);

  // x_eval is below lowest grid point => constant extrapolation
  if (x_eval <= x.front()) { return ad::value(u.front()); }

  // x_eval is above highest grid point => constant extrapolation
  if (x_eval >= x.back()) { return ad::value(u.back()); }

  const auto dx = x[1] - x[0];

  // Find first element in x that is above x_eval
  const auto first_over = std::upper_bound(std::cbegin(x), std::cend(x), x_eval);
  assert(first_over != std::cbegin(x) && first_over != std::cend(x));

  // x_eval is within the volume of the higher value,
  // i.e. x_eval in [x_over - 1/2 * dx, x_over + 1/2 * dx]
  if (*first_over - x_eval <= 0.5 * dx) {
    return ad::value(u[static_cast<size_t>(std::distance(std::cbegin(x), first_over))]);
  }

  // x_eval is within the volume of the lower value,
  // i.e. x_eval in [x_under - 1/2 * dx, x_under + 1/2 * dx]
  const auto last_under = std::prev(first_over);
  return ad::value(u[static_cast<size_t>(std::distance(std::cbegin(x), last_under))]);
}

// -------------------------------------------------------------------------------------------------
// Evaluate piecewise linear reconstruction of solution u
template <typename AT, typename PT>
[[nodiscard]] constexpr auto
piecewise_linear(const std::vector<PT>& x,  // x grid
                 const std::vector<AT>& u,  // solution at grid points at time t
                 const AT& x_eval           // x position at which u is evaluated
                 ) noexcept -> AT {
  assert(x.size() == u.size());
  assert(x.size() >= 2UL);

  // x_eval is below lowest grid point => constant extrapolation
  if (x_eval <= x.front()) { return u.front(); }

  // x_eval is above highest grid point => constant extrapolation
  if (x_eval >= x.back()) { return u.back(); }

  // Find first element in x that is above x_eval and last element under x_eval
  const auto first_over     = std::upper_bound(std::cbegin(x), std::cend(x), x_eval);
  const auto first_over_idx = std::distance(std::cbegin(x), first_over);
  assert(first_over != std::cbegin(x) && first_over != std::cend(x));

  const auto last_under     = std::prev(first_over);
  const auto last_under_idx = std::distance(std::cbegin(x), last_under);

  // Linearly interpolate between the values
  return u[static_cast<size_t>(last_under_idx)] +
         (x_eval - *last_under) / (*first_over - *last_under) *
             (u[static_cast<size_t>(first_over_idx)] - u[static_cast<size_t>(last_under_idx)]);
}

// -------------------------------------------------------------------------------------------------
// Solve Burgers equation using a Lax Friedrichs scheme
//   -> Return t, u, x1
template <typename AT, typename PT>
[[nodiscard]] constexpr auto
solve_1d_burgers(const std::vector<PT>& x,   // x grid
                 const std::vector<AT>& u0,  // initial value at t=0
                 const PT& x1_0,             // location of first (and only) shock at time t=0
                 const PT& tend,             // Final time
                 const PT& C,                // Un-smoothing parameters for derivative calculation
                 const PT& alpha             // Un-smoothing parameters for derivative calculation
                 ) -> std::tuple<std::vector<PT>, std::vector<std::vector<AT>>, std::vector<AT>> {
  assert(x.size() == u0.size() && "x and u0 must have same size");
  assert(x.size() >= 2UL && "x and u0 must have at least two entries");

  // Flux function for Burgers equation
  const auto f = [](const auto& u) constexpr { return 0.5 * u * u; };

  // Step size in x direction, assume equidistant x grid
  const auto dx = x[1] - x[0];

  // Factor for step size in t direction from Courant-Friedrichs-Lewy condition:
  //   C*dt <= dx with C = max_i |f'(u_i)|
  //   Assume u_i <= 2.0 => f'(u_i) = u_i <= 2.0 => C = 2.0 is acceptable (holds for eps <= 1)
  constexpr PT CFL_factor = 2.0;

  // Number of steps in x direction
  const auto x_steps = x.size();

  // Approximate number of steps in t direction
  const auto t_steps = static_cast<size_t>(std::floor(CFL_factor * tend / dx)) + 1UL;

  // Initialize solution of Burgers equation u
  //   -> Matrix where each row is the solution at a time step
  std::vector<std::vector<AT>> u;
  u.reserve(t_steps);
  u.emplace_back(x_steps);
  std::copy(std::cbegin(u0), std::cend(u0), std::begin(u[0]));

  // Initialize solution of shock location x1
  std::vector<AT> x1;
  x1.reserve(t_steps);
  x1.push_back(x1_0);

  // Initialize t grid
  std::vector<PT> t;
  t.reserve(t_steps);
  t.push_back(0.0);

  // Lax Friedrichs schme:
  for (size_t idx = 1UL; t.back() < tend; ++idx) {
    u.emplace_back(x_steps);
    assert(u.size() == idx + 1UL);

    const auto& u_prev = u[idx - 1UL];
    auto& u_curr       = u[idx];

    // Calculate step size in t direction via CFL condition
    const auto dt = std::min(dx / CFL_factor, tend - t.back());

    // Initialize boundary values
    u_curr[0]           = u_prev[0];
    u_curr[x_steps - 1] = u_prev[x_steps - 1];

    // Calculate u for next time step by the Lax Friedrichs scheme
    for (size_t i = 1; i < x_steps - 1; ++i) {
      u_curr[i] = 0.5 * (u_prev[i + 1] + u_prev[i - 1]) -
                  0.5 * dt / dx * (f(u_prev[i + 1]) - f(u_prev[i - 1]));
    }

    // Calculate shock loaction for the next time step
    //   -> Use piecewise constant reconstruction of u
    //   -> this calculation should not be differentiated
    // Note: We know that for Burgers equation
    //   (f(u_plus) - f(u_minus)) / (u_plus - u_minus) = 1/2 * (u_plus + u_minus)
    const PT x1_next_value = [&]() -> PT {
      const PT x_plus  = ad::value(x1.back()) + dx;  // TODO: Is this correct?
      const PT x_minus = ad::value(x1.back()) - dx;  // TODO: Is this correct?
      const PT u_plus  = piecewise_constant(x, u_prev, x_plus);
      const PT u_minus = piecewise_constant(x, u_prev, x_minus);
      return ad::value(x1.back()) + dt * static_cast<PT>(0.5) * (u_plus + u_minus);
    }();

    // Calculate smoothed shock location for next time step
    //   -> Use piecewise linear reconstruction of u
    //   -> Only use the result for derivative calculation -> Overwrite value later
    AT x1_next = [&]() -> AT {
      const AT x_plus  = x1.back() + C * std::pow(dx, alpha);  // TODO: Is this correct?
      const AT x_minus = x1.back() - C * std::pow(dx, alpha);  // TODO: Is this correct?
      const AT u_plus  = piecewise_linear(x, u_prev, x_plus);
      const AT u_minus = piecewise_linear(x, u_prev, x_minus);
      return x1.back() + dt * static_cast<PT>(0.5) * (u_plus + u_minus);
    }();
    ad::value(x1_next) = x1_next_value;

    x1.push_back(x1_next);

    // Update t grid
    t.push_back(t.back() + dt);
  }

  return {t, u, x1};
}

}  // namespace Zap::ShockAD_1D

#endif  // ZAP_SHOCK_AD_1D_SOLVER_HPP_
