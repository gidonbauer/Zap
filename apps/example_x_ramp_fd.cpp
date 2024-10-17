#include "CellBased/EigenDecomp.hpp"
#include "CellBased/Solver.hpp"
#include "IO/NoopWriter.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Timer.hpp"

using namespace Zap::CellBased;
using Zap::IO::NoopWriter;

// - Setup -----------------------------------------------------------------------------------------
using Float           = double;
constexpr size_t DIM  = 1;
constexpr Float X_MIN = 0.0;
constexpr Float X_MAX = 2.0;
constexpr Float Y_MIN = 0.0;
constexpr Float Y_MAX = 2.0;

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto chi(Float x, Float x_min, Float x_max) noexcept -> Float {
  return static_cast<Float>(x >= x_min && x <= x_max);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto x_shock(Float t, Float eps) noexcept -> Float {
  return std::sqrt(1 + (1 + eps) * t);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto xi_shock(Float t) noexcept -> Float {
  return t / (2 * std::sqrt(1 + t));
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto u_eps(Float x, Float t, Float eps) noexcept -> Float {
  return (1 + eps) * x / (1 + (1 + eps) * t) * chi(x, Float{0}, x_shock(t, eps));
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto v(Float x, Float t) noexcept -> Float {
  return x / ((1 + t) * (1 + t)) * chi(x, 0, x_shock(t, Float{0}));
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto run(Float eps, size_t nx, size_t ny, Float tend) noexcept
    -> UniformGrid<Float, Float, DIM> {
  Igor::ScopeTimer timer(std::format("Solver (nx={}, ny={}, tend={}, eps={})", nx, ny, tend, eps));

  UniformGrid<Float, Float, DIM> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.same_value_boundary();

  {
    std::vector<SimCoord<Float>> points{
        {x_shock(0, eps), 0.0},
        {x_shock(0, eps), (Y_MAX - Y_MIN) + Y_MIN},
    };
    if (!grid.cut_piecewise_linear(points)) {
      Igor::Panic("Could not cut the grid along the initial shock.");
    }
  }

  auto u0 = [eps](Float x, Float /*y*/) { return u_eps(x, 0, eps); };
  grid.fill_four_point(u0);

  NoopWriter grid_writer;
  NoopWriter t_writer;

  auto solver    = make_solver<ExtendType::MAX>(SingleEq::A{}, SingleEq::B{});
  const auto res = solver.solve(grid, tend, grid_writer, t_writer, 0.25);
  IGOR_ASSERT(res.has_value(), "Expected the solver to succeed, but did not.");
  return *res;
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto v_fd(const UniformGrid<Float, Float, DIM>& sol_p_eps,
                        const UniformGrid<Float, Float, DIM>& sol_m_eps,
                        Float eps,
                        SimCoord<Float> p) -> Float {
  return ((sol_p_eps.eval(p) - sol_m_eps.eval(p)) / (2 * eps))(0);
}

// -------------------------------------------------------------------------------------------------
template <typename FUNC>
[[nodiscard]] auto L1_norm(FUNC f, size_t n) noexcept -> Float {
  // Use simpsons rule
  if (n % 2 == 1) { n += 1; }
  Float L1       = 0;
  const Float dx = (X_MAX - X_MIN) / static_cast<Float>(n);
  for (size_t i = 1; i <= n / 2; ++i) {
    L1 += std::abs(f(static_cast<Float>(2 * i - 2) * dx + X_MIN)) +
          4 * std::abs(f(static_cast<Float>(2 * i - 1) * dx + X_MIN)) +
          std::abs(f(static_cast<Float>(2 * i) * dx + X_MIN));
  }
  return L1 * dx / 3;
}

// -------------------------------------------------------------------------------------------------
auto main() -> int {
  const size_t nx   = 5000;
  const size_t ny   = 5;
  const Float tend  = 1.0;
  const size_t n_L1 = 100'000;

  const Float eps = 1e-8;

  const auto sol       = run(0, nx, ny, tend);
  const auto sol_p_eps = run(eps, nx, ny, tend);
  const auto sol_m_eps = run(-eps, nx, ny, tend);

  const auto shock_curve       = sol.get_shock_curve();
  const auto shock_curve_p_eps = sol_p_eps.get_shock_curve();
  const auto shock_curve_m_eps = sol_m_eps.get_shock_curve();

  const auto avg_shock_x = std::transform_reduce(std::cbegin(shock_curve),
                                                 std::cend(shock_curve),
                                                 Float{0},
                                                 std::plus<>{},
                                                 [](const auto& p) { return p.x; }) /
                           static_cast<Float>(shock_curve.size());
  const auto avg_shock_x_p_eps = std::transform_reduce(std::cbegin(shock_curve_p_eps),
                                                       std::cend(shock_curve_p_eps),
                                                       Float{0},
                                                       std::plus<>{},
                                                       [](const auto& p) { return p.x; }) /
                                 static_cast<Float>(shock_curve_p_eps.size());
  const auto avg_shock_x_m_eps = std::transform_reduce(std::cbegin(shock_curve_m_eps),
                                                       std::cend(shock_curve_m_eps),
                                                       Float{0},
                                                       std::plus<>{},
                                                       [](const auto& p) { return p.x; }) /
                                 static_cast<Float>(shock_curve_m_eps.size());

  const auto xi_shock_fd = (avg_shock_x_p_eps - avg_shock_x_m_eps) / (2 * eps);

  const auto L1_err_u = L1_norm(
      [&sol, tend](Float x) { return sol.eval(SimCoord<Float>{x, 1.0})(0) - u_eps(x, tend, 0.0); },
      n_L1);
  const auto L1_u = L1_norm([&sol](Float x) { return sol.eval(SimCoord<Float>{x, 1.0})(0); }, n_L1);
  const auto L1_u_ex = L1_norm([tend](Float x) { return u_eps(x, tend, 0.0); }, n_L1);

  const auto L1_err_v = L1_norm(
      [&sol_m_eps, &sol_p_eps, eps, tend](Float x) {
        return v_fd(sol_p_eps, sol_m_eps, eps, SimCoord<Float>{x, 1.0}) - v(x, tend);
      },
      n_L1);
  const auto L1_v_fd =
      L1_norm([&sol_m_eps, &sol_p_eps, eps](
                  Float x) { return v_fd(sol_p_eps, sol_m_eps, eps, SimCoord<Float>{x, 1.0}); },
              n_L1);
  const auto L1_v_ex = L1_norm([tend](Float x) { return v(x, tend); }, n_L1);

  std::cout << '\n';
  Igor::Info("tend = {}", tend);
  Igor::Info("eps  = {}", eps);
  Igor::Info("dx   = {}", (X_MAX - X_MIN) / static_cast<Float>(nx));

  std::cout << '\n';
  Igor::Info("||u - u_ex||_L1                  = {:.6f}", L1_err_u);
  Igor::Info("||u - u_ex||_L1 / ||u_ex||_L1    = {:.6f}", L1_err_u / L1_u_ex);
  Igor::Info("||u||_L1                         = {:.6f}", L1_u);
  Igor::Info("||u_ex||_L1                      = {:.6f}", L1_u_ex);

  std::cout << '\n';
  Igor::Info("||v_fd - v_ex||_L1               = {:.6f}", L1_err_v);
  Igor::Info("||v_fd - v_ex||_L1 / ||v_ex||_L1 = {:.6f}", L1_err_v / L1_v_ex);
  Igor::Info("||v_fd||_L1                      = {:.6f}", L1_v_fd);
  Igor::Info("||v_ex||_L1                      = {:.6f}", L1_v_ex);

  std::cout << '\n';
  Igor::Info("|x_fd - x_ex|                    = {:.6f}", std::abs(avg_shock_x - x_shock(tend, 0)));
  Igor::Info("|x_fd - x_ex| / x_ex             = {:.6f}",
             std::abs(avg_shock_x - x_shock(tend, 0)) / x_shock(tend, 0));
  Igor::Info("x_fd                             = {:.6f}", avg_shock_x);
  Igor::Info("x_ex                             = {:.6f}", x_shock(tend, 0));

  std::cout << '\n';
  Igor::Info("|xi_fd - xi_ex|                  = {:.6f}", std::abs(xi_shock_fd - xi_shock(tend)));
  Igor::Info("|xi_fd - xi_ex| / xi_ex          = {:.6f}",
             std::abs(xi_shock_fd - xi_shock(tend)) / xi_shock(tend));
  Igor::Info("xi_fd                            = {:.6f}", xi_shock_fd);
  Igor::Info("xi_ex                            = {:.6f}", xi_shock(tend));
}
