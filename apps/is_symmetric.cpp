#include <numbers>

#include <boost/multiprecision/cpp_bin_float.hpp>
namespace mp = boost::multiprecision;

#include <fmt/format.h>
#include <fmt/ostream.h>

#include "Igor/Logging.hpp"
#include "Igor/Timer.hpp"

// - Setup -----------------------------------------------------------------------------------------
using PassiveFloat = mp::number<mp::cpp_bin_float<128>, mp::expression_template_option::et_on>;
using ActiveFloat  = PassiveFloat;
const PassiveFloat X_MIN = 0.0;
const PassiveFloat X_MAX = 5.0;
const PassiveFloat Y_MIN = 0.0;
const PassiveFloat Y_MAX = 5.0;
const PassiveFloat R     = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;

// - Multiprecision adapter ------------------------------------------------------------------------
template <>
struct fmt::formatter<PassiveFloat> : fmt::ostream_formatter {};

namespace std {

using mp::abs;
using mp::atan2;
using mp::ceil;
using mp::cos;
using mp::floor;
using mp::isinf;
using mp::isnan;
using mp::round;
using mp::sin;
using mp::sqrt;

}  // namespace std
// - Multiprecision adapter ------------------------------------------------------------------------

// #define ZAP_SERIAL
// #define ZAP_SINGLE_ITERATION
#define ZAP_NO_TANGENTIAL_CORRECTION
#define ZAP_NO_CHUNK_INFO
#include "CellBased/Solver.hpp"
#include "IO/NoopWriter.hpp"

#include "Igor/Logging.hpp"

using namespace Zap::CellBased;
using namespace Zap::IO;

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto sqr(auto x) noexcept { return x * x; }

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto u0(PassiveFloat x, PassiveFloat y) noexcept -> ActiveFloat {
  return (sqr(x - X_MIN) + sqr(y - Y_MIN)) *
         static_cast<ActiveFloat>((sqr(x - X_MIN) + sqr(y - Y_MIN)) <= sqr(R));
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto init_shock(PassiveFloat t) noexcept -> Zap::CellBased::SimCoord<ActiveFloat> {
  return {
      R * std::cos(std::numbers::pi / 2 * t),
      R * std::sin(std::numbers::pi / 2 * t),
  };
}

// -------------------------------------------------------------------------------------------------
struct SymmetryResults {
  size_t num_non_symmetric;
  size_t num_points_checked;
  PassiveFloat non_symmetric_error;
  PassiveFloat entire_mass;
};

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto check_symmetry(const UniformGrid<ActiveFloat, PassiveFloat>& numerical_solution,
                                  size_t N         = 2000,
                                  PassiveFloat tol = EPS<PassiveFloat>()) noexcept
    -> SymmetryResults {
  // Igor::ScopeTimer timer("Symmetry check");

  size_t num_non_symmetric         = 0;
  size_t num_points_checked        = 0;
  PassiveFloat non_symmetric_error = 0.0;

#pragma omp parallel for schedule(dynamic) reduction(+ : num_non_symmetric)                        \
    reduction(+ : num_points_checked) reduction(+ : non_symmetric_error)
  for (size_t yi = 0; yi < N; ++yi) {
    for (size_t xi = yi + 1; xi < N; ++xi) {
      const PassiveFloat x = X_MIN + (X_MAX - X_MIN) * static_cast<PassiveFloat>(xi) /
                                         static_cast<PassiveFloat>(N - 1);
      const PassiveFloat y = Y_MIN + (Y_MAX - Y_MIN) * static_cast<PassiveFloat>(yi) /
                                         static_cast<PassiveFloat>(N - 1);

      const auto a = numerical_solution.eval(SimCoord(x, y));
      const auto b = numerical_solution.eval(SimCoord(y, x));
      if (!approx_eq(a, b, tol)) {
        num_non_symmetric += 1;
        non_symmetric_error += std::abs(a - b);
      }
      num_points_checked += 1;
    }
  }
  non_symmetric_error /= static_cast<PassiveFloat>(num_points_checked);

  return SymmetryResults{
      .num_non_symmetric   = num_non_symmetric,
      .num_points_checked  = num_points_checked,
      .non_symmetric_error = non_symmetric_error,
      .entire_mass         = numerical_solution.mass(),
  };
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto run(size_t nx, size_t ny, PassiveFloat tend, PassiveFloat cfl) noexcept
    -> std::optional<UniformGrid<ActiveFloat, PassiveFloat>> {
  UniformGrid<ActiveFloat, PassiveFloat> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.periodic_boundary();

  if (!grid.cut_curve(init_shock)) {
    Igor::Warn("Could not cut initital shock.");
    return std::nullopt;
  }

  grid.fill_four_point(u0);

  // if (check_symmetry(grid).num_non_symmetric != 0UZ) {
  //   Igor::Warn("Initital condition for nx={} and ny={} is non-symmetric.", nx, ny);
  //   return std::nullopt;
  // }

  Solver<ExtendType::NEAREST> solver;
  NoopWriter grid_writer;
  NoopWriter time_writer;

  // Igor::ScopeTimer timer("Solver");
  return solver.solve(grid, tend, grid_writer, time_writer, cfl);
}

// -------------------------------------------------------------------------------------------------
auto main() -> int {
  const PassiveFloat tend = 1.0;
  const PassiveFloat tol  = 1e-4;
  const size_t N          = 500UZ;

  // Igor::Info("| Gridsize | #non-symmetric | %non-symmetric |  %cut  | abs. non-symmetric error |
  // "
  //            "rel. non-symmetric error |");
  // Igor::Info("-------------------------------------------------------------------------------------"
  //            "------------------------");
  Igor::Info("Gridsize,#non-symmetric,%non-symmetric,%cut,abs. non-symmetric error,"
             "rel. non-symmetric error");
  for (size_t n = 11; n <= 401; n += n > 150 ? 50 : 10) {
    const PassiveFloat cfl = n <= 101 ? 0.5 : 0.2;
    const auto res         = run(n, n, tend, cfl);
    if (!res.has_value()) {
      Igor::Warn("Solver for nx={}, ny={}, tend={} failed.", n, n, tend);
    } else {
      const auto sym_res = check_symmetry(*res, N, tol);
      // Igor::Info("| {:>3}x{:<3}  |"
      //            "    {:>8}    |"
      //            "     {:>5.2f}%     |"
      //            " {:>5.2f}% |"
      //            "       {:>12.8f}       |"
      //            "       {:>12.8f}       |",
      //            n,
      //            n,
      //            sym_res.num_non_symmetric,
      //            static_cast<PassiveFloat>(sym_res.num_non_symmetric * 100) /
      //                static_cast<PassiveFloat>(sym_res.num_points_checked),
      //            static_cast<PassiveFloat>(res->cut_cell_idxs().size() * 100) /
      //                static_cast<PassiveFloat>(res->size()),
      //            sym_res.non_symmetric_error,
      //            sym_res.non_symmetric_error / sym_res.entire_mass);
      Igor::Info(
          "{}x{},{},{}%,{}%,{},{}",
          n,
          n,
          sym_res.num_non_symmetric,
          static_cast<PassiveFloat>(static_cast<PassiveFloat>(sym_res.num_non_symmetric * 100) /
                                    static_cast<PassiveFloat>(sym_res.num_points_checked)),
          static_cast<PassiveFloat>(static_cast<PassiveFloat>(res->cut_cell_idxs().size() * 100) /
                                    static_cast<PassiveFloat>(res->size())),
          sym_res.non_symmetric_error,
          static_cast<PassiveFloat>(sym_res.non_symmetric_error / sym_res.entire_mass));
    }
    std::cout << std::flush;
  }
}
