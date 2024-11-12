#include <numbers>

// #define ZAP_2ND_ORDER_CORRECTION

#include "CellBased/Grid.hpp"
#include "CellBased/Solver.hpp"

#include "MatBased/BoundaryConditions.hpp"
#include "MatBased/Solver.hpp"

#include "IO/NoopWriter.hpp"
#include "IO/ReadClawpackClassic.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Macros.hpp"
#include "Igor/Timer.hpp"

#include "Common.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "compare_clawpack/"

using namespace Zap::IO;
using namespace Zap::CellBased;

// - Setup -----------------------------------------------------------------------------------------
using Float           = double;
constexpr Float X_MIN = 0.0;
constexpr Float X_MAX = 5.0;
constexpr Float Y_MIN = 0.0;
constexpr Float Y_MAX = 5.0;
constexpr Float R     = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;
constexpr Float TEND  = 1.0;

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto sqr(auto x) noexcept { return x * x; }

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto u0(Float x, Float y) noexcept -> Float {
  return (sqr(x - X_MIN) + sqr(y - Y_MIN)) *
         static_cast<Float>((sqr(x - X_MIN) + sqr(y - Y_MIN)) <= sqr(R));
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto init_shock(Float t) noexcept -> SimCoord<Float> {
  return {
      R * std::cos(std::numbers::pi_v<Float> / 2 * t),
      R * std::sin(std::numbers::pi_v<Float> / 2 * t),
  };
};

// -------------------------------------------------------------------------------------------------
auto run_cell_based(size_t nx, size_t ny, Float tend, Float CFL_safety_factor) noexcept
    -> std::optional<Zap::CellBased::UniformGrid<Float, Float>> {
  Igor::ScopeTimer timer(fmt::format("Cut-cell solver with nx={}, ny={}, tend={}", nx, ny, tend));

  Zap::CellBased::UniformGrid<Float, Float> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);

  grid.periodic_boundary();
  if (!grid.cut_curve(init_shock)) { return std::nullopt; }
  grid.fill_four_point(u0);

  Zap::CellBased::Solver<Zap::CellBased::ExtendType::NEAREST> solver;

  NoopWriter u_writer;
  NoopWriter t_writer;
  return solver.solve(grid, tend, u_writer, t_writer, CFL_safety_factor);
}

// -------------------------------------------------------------------------------------------------
auto run_mat_based(size_t nx, size_t ny, Float tend, Float CFL_safety_factor) noexcept
    -> std::optional<std::tuple<Zap::MatBased::Vector<Float>,
                                Zap::MatBased::Vector<Float>,
                                Zap::MatBased::Matrix<Float>>> {
  const auto dx = (X_MAX - X_MIN) / static_cast<Float>(nx);
  const auto dy = (Y_MAX - Y_MIN) / static_cast<Float>(ny);

  Zap::MatBased::Vector<Float> x(nx);
  std::generate(std::begin(x), std::end(x), [i = 0, dx]() mutable {
    return X_MIN + dx * (static_cast<Float>(i++) + static_cast<Float>(0.5));
  });
  Zap::MatBased::Vector<Float> y(static_cast<int>(ny));
  std::generate(std::begin(y), std::end(y), [i = 0, dy]() mutable {
    return Y_MIN + dy * (static_cast<Float>(i++) + static_cast<Float>(0.5));
  });

  Zap::MatBased::Matrix<Float> u0_grid(static_cast<int>(ny), static_cast<int>(nx));
  for (int yi = 0; yi < static_cast<int>(ny); ++yi) {
    for (int xi = 0; xi < static_cast<int>(nx); ++xi) {
      // - Fill four point -------------
      u0_grid(yi, xi) = (u0(x(xi) - dx / 4, y(yi) - dy / 4) + u0(x(xi) + dx / 4, y(yi) - dy / 4) +
                         u0(x(xi) - dx / 4, y(yi) + dy / 4) + u0(x(xi) + dx / 4, y(yi) + dy / 4)) /
                        4;
    }
  }

  auto boundary = [](Zap::MatBased::Matrix<Float>& u_next,
                     const Zap::MatBased::Matrix<Float>& u_curr,
                     Float dt,
                     Float local_dx,
                     Float local_dy,
                     auto numerical_flux_x,
                     auto numerical_flux_y) {
    return Zap::MatBased::periodic_boundary(
        u_next, u_curr, dt, local_dx, local_dy, numerical_flux_x, numerical_flux_y);
  };

  Igor::ScopeTimer timer(fmt::format("Godunov solver with nx={}, ny={}, tend={}", nx, ny, tend));

  Zap::IO::NoopWriter u_writer;
  Zap::IO::NoopWriter t_writer;
  const auto res = Zap::MatBased::solve_2d_burgers(
      x, y, u0_grid, tend, boundary, u_writer, t_writer, CFL_safety_factor);

  if (!res.has_value()) { return std::nullopt; }
  return std::make_tuple(x, y, *res);
}

// -------------------------------------------------------------------------------------------------
constexpr auto eval_mat_at(Float x,
                           Float y,
                           const Zap::MatBased::Vector<Float>& xs,
                           const Zap::MatBased::Vector<Float>& ys,
                           const Zap::MatBased::Matrix<Float>& u) noexcept -> Float {
  assert(ys.size() == u.rows());
  assert(xs.size() == u.cols());
  assert(xs.size() >= 2);
  assert(ys.size() >= 2);

  if (std::isnan(x) || std::isnan(y)) { return std::numeric_limits<Float>::quiet_NaN(); }

  const auto dx    = xs(1) - xs(0);
  const auto x_idx = static_cast<Eigen::Index>(std::clamp(
      std::round((x - xs(0)) / dx), static_cast<Float>(0), static_cast<Float>(xs.rows() - 1)));

  const auto dy    = ys(1) - ys(0);
  const auto y_idx = static_cast<Eigen::Index>(std::clamp(
      std::round((y - ys(0)) / dy), static_cast<Float>(0), static_cast<Float>(ys.rows() - 1)));

  return u(y_idx, x_idx);
}
// -------------------------------------------------------------------------------------------------
void compare(const ClawpackSolution<Float>& refernce_sol,
             size_t nx,
             size_t ny,
             Float tend,
             Float CFL_safety_factor,
             std::ostream& out) noexcept {
  out << nx << 'x' << ny << ":\n";
  {
    const auto sol = run_cell_based(nx, ny, tend, CFL_safety_factor);
    if (!sol.has_value()) {
      Igor::Warn("Solver with nx={}, ny={}, tend={}", nx, ny, tend);
      return;
    }

    const auto abs_diff = [&](Float x, Float y) {
      return std::abs(refernce_sol(x, y) - sol->eval(SimCoord<Float>{x, y}));
    };

    const auto L1_error = simpsons_rule_2d(abs_diff, X_MIN, X_MAX, Y_MIN, Y_MAX, 2000, 2000);
    // Igor::Info("Cell based: nx={}, ny={}, tend={} => L1={}", nx, ny, tend, L1_error);
    out << "L1_error_cell = " << L1_error << '\n';
  }
  {
    const auto sol = run_mat_based(nx, ny, tend, CFL_safety_factor);
    if (!sol.has_value()) {
      Igor::Warn("Solver with nx={}, ny={}, tend={}", nx, ny, tend);
      return;
    }

    const auto abs_diff = [&](Float x, Float y) {
      return std::abs(refernce_sol(x, y) -
                      eval_mat_at(x, y, std::get<0>(*sol), std::get<1>(*sol), std::get<2>(*sol)));
    };

    const auto L1_error = simpsons_rule_2d(abs_diff, X_MIN, X_MAX, Y_MIN, Y_MAX, 2000, 2000);
    // Igor::Info("Mat based: nx={}, ny={}, tend={} => L1={}", nx, ny, tend, L1_error);
    out << "L1_error_mat = " << L1_error << '\n';
  }
  out << "------------------------------------------------------------\n" << std::flush;
  std::cout << '\n';
}

// -------------------------------------------------------------------------------------------------
auto main() -> int {
  const auto refernce_sol = [] {
    try {
      // Solution for nx = 1000, ny = 1000, tend = 1.0
      const std::string q_filename = "../apps/clawpack_solution/fort.q0050";
      const std::string t_filename = "../apps/clawpack_solution/fort.t0050";

      return ClawpackSolution<Float>(q_filename, t_filename);
    } catch (const std::exception& e) {
      Igor::Panic("Could not read clawpack solution: {}", e.what());
      std::unreachable();
    }
  }();

  IGOR_ASSERT(approx_eq(X_MIN, refernce_sol.x_min()) && approx_eq(X_MAX, refernce_sol.x_max()) &&
                  approx_eq(Y_MIN, refernce_sol.y_min()) && approx_eq(Y_MAX, refernce_sol.y_max()),
              "Setup does not match with reference solution.");

  constexpr std::array ns = {
      5UZ,
      11UZ,
      // 15UZ,
      // 21UZ,
      // 31UZ,
      41UZ,
      51UZ,
      75UZ,
      101UZ,
      151UZ,
      201UZ,
      251UZ,
      301UZ,
      351UZ,
      401UZ,
      // 451UZ,
      // 501UZ,
  };

  constexpr auto result_file = "clawpack_compare_res.txt";
  std::ofstream out(result_file);
  if (!out) { Igor::Warn("Could not open file `{}`: {}", result_file, std::strerror(errno)); }

  for (size_t n : ns) {
    compare(refernce_sol, n, n, TEND, 0.2, out);
  }
  Igor::Info("Saved results to `{}`.", result_file);
}
