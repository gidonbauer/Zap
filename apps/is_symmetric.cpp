#include <numbers>

// #define ZAP_SERIAL
// #define ZAP_SINGLE_ITERATION
#define ZAP_NO_CHUNK_INFO
#include "CellBased/Solver.hpp"
#include "IO/NoopWriter.hpp"

#include "Igor/Logging.hpp"

using namespace Zap::CellBased;
using namespace Zap::IO;

// - Setup -----------------------------------------------------------------------------------------
using ActiveFloat            = double;
using PassiveFloat           = double;
constexpr PassiveFloat X_MIN = 0.0;
constexpr PassiveFloat X_MAX = 5.0;
constexpr PassiveFloat Y_MIN = 0.0;
constexpr PassiveFloat Y_MAX = 5.0;

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto sqr(auto x) noexcept { return x * x; }

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto u0(PassiveFloat x, PassiveFloat y) noexcept -> ActiveFloat {
  return (sqr(x - X_MIN) + sqr(y - Y_MIN)) *
         static_cast<ActiveFloat>((sqr(x - X_MIN) + sqr(y - Y_MIN)) <=
                                  sqr((X_MIN + X_MAX + Y_MIN + Y_MAX) / 4));
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto init_shock(PassiveFloat t) noexcept -> Zap::CellBased::SimCoord<ActiveFloat> {
  const auto r = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;
  return {
      r * std::cos(std::numbers::pi_v<PassiveFloat> / 2 * t),
      r * std::sin(std::numbers::pi_v<PassiveFloat> / 2 * t),
  };
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto
num_non_symmetric(const UniformGrid<ActiveFloat, PassiveFloat>& numerical_solution,
                  size_t N         = 2000,
                  PassiveFloat tol = EPS<PassiveFloat>) noexcept -> size_t {
  size_t n = 0UZ;
  for (size_t yi = 0; yi < N; ++yi) {
    for (size_t xi = 0; xi < N; ++xi) {
      const PassiveFloat x = X_MIN + (X_MAX - X_MIN) * static_cast<PassiveFloat>(xi) /
                                         static_cast<PassiveFloat>(N - 1);
      const PassiveFloat y = Y_MIN + (Y_MAX - Y_MIN) * static_cast<PassiveFloat>(yi) /
                                         static_cast<PassiveFloat>(N - 1);

      const auto a = numerical_solution.eval(SimCoord(x, y));
      const auto b = numerical_solution.eval(SimCoord(y, x));
      n += static_cast<size_t>(!approx_eq(a, b, tol));
      // if (!approx_eq(a, b, 1e-2)) {
      //   Igor::Debug("Numerical solution is not symmetric.");
      //   Igor::Debug("f({}, {}) = {}", x, y, a);
      //   Igor::Debug("f({}, {}) = {}", y, x, b);
      //   return false;
      // }
    }
  }
  return n;
}

[[nodiscard]] constexpr auto run(size_t nx, size_t ny, PassiveFloat tend) noexcept
    -> std::optional<UniformGrid<ActiveFloat, PassiveFloat>> {
  UniformGrid<ActiveFloat, PassiveFloat> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.periodic_boundary();

  if (!grid.cut_curve(init_shock)) {
    Igor::Warn("Could not cut initital shock.");
    return std::nullopt;
  }

  grid.fill_four_point(u0);

  if (num_non_symmetric(grid) != 0UZ) {
    Igor::Warn("Initital condition for nx={} and ny={} is non-symmetric.", nx, ny);
    return std::nullopt;
  }

  Solver<ExtendType::NEAREST> solver;
  NoopWriter grid_writer;
  NoopWriter time_writer;
  return solver.solve(grid, tend, grid_writer, time_writer);
}

// -------------------------------------------------------------------------------------------------
auto main() -> int {
  const PassiveFloat tend = 1.0;
  const PassiveFloat tol  = 1e-4;

  Igor::Info("Number of non-symmetric points:");
  for (size_t n = 11; n <= 401; n += n > 100 ? 50 : 10) {
    const auto res = run(n, n, tend);
    if (!res.has_value()) {
      Igor::Warn("Solver for nx={}, ny={}, tend={} failed.", n, n, tend);
    } else {
      const auto N   = 2000UZ;
      const auto nns = num_non_symmetric(*res, N, tol);
      Igor::Info("{:>3}x{:<3} => {:>7} ({:>5.2f}%) ({:>5.2f}%)",
                 n,
                 n,
                 nns,
                 static_cast<PassiveFloat>(nns * 100) / sqr(static_cast<PassiveFloat>(N)),
                 static_cast<PassiveFloat>(res->cut_cell_idxs().size() * 100) /
                     static_cast<PassiveFloat>(res->size()));
    }
  }
}
