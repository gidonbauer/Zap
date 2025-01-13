#include <benchmark/benchmark.h>

#include <cmath>
#include <numbers>

#define ZAP_NO_CHUNK_INFO
// #define ZAP_STATIC_CUT
// #define ZAP_TANGENTIAL_CORRECTION

#include "CellBased/Grid.hpp"
#include "CellBased/Solver.hpp"

#include "IO/NoopWriter.hpp"

using namespace Zap::IO;
using namespace Zap::CellBased;

// - Setup -----------------------------------------------------------------------------------------
using ActiveFloat                        = double;
using PassiveFloat                       = double;
constexpr PassiveFloat X_MIN             = 0.0;
constexpr PassiveFloat X_MAX             = 5.0;
constexpr PassiveFloat Y_MIN             = 0.0;
constexpr PassiveFloat Y_MAX             = 5.0;
constexpr PassiveFloat TEND              = 1.0;
constexpr PassiveFloat CFL_SAFETY_FACTOR = 0.2;

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
static void run_cell_solver(benchmark::State& state) noexcept {
  IGOR_ASSERT(
      state.range(0) > 0, "Expected range value to be greater than zero, is {}", state.range(0));
  const auto n = static_cast<size_t>(state.range(0));
  for (auto _ : state) {
    state.PauseTiming();
    Zap::CellBased::UniformGrid<ActiveFloat, PassiveFloat> grid(X_MIN, X_MAX, n, Y_MIN, Y_MAX, n);
    // grid.same_value_boundary();
    grid.periodic_boundary();
    if (!grid.cut_curve(init_shock)) {
      Igor::Panic("Could not cut grid for nx={} and ny={}", n, n);
    }
    grid.fill_four_point(u0);
    state.ResumeTiming();

    NoopWriter u_writer;
    NoopWriter t_writer;
    const auto sol = solve_2d_burgers<Zap::CellBased::ExtendType::NEAREST>(
        grid, TEND, u_writer, t_writer, CFL_SAFETY_FACTOR);
    if (!sol.has_value()) { Igor::Panic("Solver failed for nx={} and ny={}", n, n); }
  }
}

BENCHMARK(run_cell_solver)
    ->Unit(benchmark::kSecond)
    ->Iterations(10)
    ->DenseRange(11, 101, 10)
    ->DenseRange(151, 401, 50);

BENCHMARK_MAIN();
