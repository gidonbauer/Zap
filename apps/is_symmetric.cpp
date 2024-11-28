#include <numbers>

#ifdef USE_MP
#include <boost/multiprecision/cpp_bin_float.hpp>
namespace mp = boost::multiprecision;
#endif  // USE_MP

#include <fmt/format.h>
#include <fmt/ostream.h>

#include "Igor/Logging.hpp"
#include "Igor/Timer.hpp"

#ifdef USE_MP

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

#else

// - Setup -----------------------------------------------------------------------------------------
using PassiveFloat           = double;
using ActiveFloat            = PassiveFloat;
constexpr PassiveFloat X_MIN = 0.0;
constexpr PassiveFloat X_MAX = 5.0;
constexpr PassiveFloat Y_MIN = 0.0;
constexpr PassiveFloat Y_MAX = 5.0;
constexpr PassiveFloat R     = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;

#endif  // USE_MP

// #define ZAP_SERIAL
// #define ZAP_SINGLE_ITERATION
// #define ZAP_STATIC_CUT
// #define ZAP_NO_TANGENTIAL_CORRECTION
#define ZAP_NO_CHUNK_INFO
#include "CellBased/Solver.hpp"
#include "IO/NoopWriter.hpp"

#include "Igor/Logging.hpp"

#include "Common.hpp"

using namespace Zap::CellBased;
using namespace Zap::IO;

// -------------------------------------------------------------------------------------------------
struct Args {
  PassiveFloat tend = 1.0;
  PassiveFloat tol  = 1e-6;
  size_t N          = 2000UZ;
  bool plus_one     = true;
};

// - Print usage -----------------------------------------------------------------------------------
template <Igor::detail::Level level>
void usage(std::string_view prog, std::ostream& out) noexcept {
  Args args{};

  out << Igor::detail::level_repr(level) << "Usage: " << prog
      << " [--tend tend] [--tol tol] [--N N] [--no-plus-one]\n";
  out << "\t--tend            Final time for simulation, default is " << args.tend << '\n';
  out << "\t--tol             Tolerance to classiyfy a value as equal, default is " << args.tol
      << '\n';
  out << "\t--N               Number of points in x- and y-direction to check for symmetry, "
         "default is "
      << args.N << '\n';
  out << "\t--no-plus-one     Do not make number of cells odd, default is " << std::boolalpha
      << !args.plus_one << '\n';
}

// - Parse command line arguments ------------------------------------------------------------------
[[nodiscard]] auto parse_args(int& argc, char**& argv) noexcept -> std::optional<Args> {
  Args args{};

  const auto prog = next_arg(argc, argv);
  IGOR_ASSERT(prog.has_value(), "Could not get program name from command line arguments.");

  for (auto arg = next_arg(argc, argv); arg.has_value(); arg = next_arg(argc, argv)) {
    using namespace std::string_view_literals;
    if (arg == "--tend"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --tend.\n";
        return std::nullopt;
      }
      args.tend = parse_double(*arg);
    } else if (arg == "--tol"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --tol.\n";
        return std::nullopt;
      }
      args.tol = parse_double(*arg);
    } else if (arg == "--N"sv || arg == "-N") {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --N.\n";
        return std::nullopt;
      }
      args.N = parse_size_t(*arg);
    } else if (arg == "--no-plus-one") {
      args.plus_one = false;
    } else if (arg == "-h"sv || arg == "--help"sv) {
      usage<Igor::detail::Level::INFO>(*prog, std::cout);
      std::exit(0);
    } else {
      usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
      std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC) << "Unknown argument `"
                << *arg << "`.\n";
      return std::nullopt;
    }
  }
  return args;
}

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

  NoopWriter grid_writer;
  NoopWriter time_writer;

  // Igor::ScopeTimer timer("Solver");
  return solve_2d_burgers<ExtendType::NEAREST>(grid, tend, grid_writer, time_writer, cfl);
}

// -------------------------------------------------------------------------------------------------
auto main(int argc, char** argv) -> int {
  const auto args = parse_args(argc, argv);
  if (!args.has_value()) { return 1; }

  Igor::Info("tend     = {}", args->tend);
  Igor::Info("tol      = {}", args->tol);
  Igor::Info("N        = {}", args->N);
  Igor::Info("plus-one = {}", args->plus_one);

#ifndef USE_MP
  Igor::Info("| Gridsize | #non-symmetric | %non-symmetric |  %cut  | abs. non-symmetric error | "
             "rel. non-symmetric error |");
  Igor::Info("-------------------------------------------------------------------------------------"
             "------------------------");
#else
  Igor::Info("Gridsize,#non-symmetric,%non-symmetric,%cut,abs. non-symmetric error,"
             "rel. non-symmetric error");
#endif  // USE_MP
  for (size_t n = 10 + static_cast<size_t>(args->plus_one);
       n <= 400 + static_cast<size_t>(args->plus_one);
       n += n > 150 ? 50 : 10) {
    const PassiveFloat cfl = n < 50 ? 0.5 : 0.2;
    const auto res         = run(n, n, args->tend, cfl);
    if (!res.has_value()) {
      Igor::Warn("Solver for nx={}, ny={}, tend={} failed.", n, n, args->tend);
    } else {
      const auto sym_res = check_symmetry(*res, args->N, args->tol);
#ifndef USE_MP
      Igor::Info("| {:>3}x{:<3}  |"
                 "    {:>8}    |"
                 "     {:>5.2f}%     |"
                 " {:>5.2f}% |"
                 "       {:>12.8f}       |"
                 "       {:>12.8f}       |",
                 n,
                 n,
                 sym_res.num_non_symmetric,
                 static_cast<PassiveFloat>(sym_res.num_non_symmetric * 100) /
                     static_cast<PassiveFloat>(sym_res.num_points_checked),
                 static_cast<PassiveFloat>(res->cut_cell_idxs().size() * 100) /
                     static_cast<PassiveFloat>(res->size()),
                 sym_res.non_symmetric_error,
                 sym_res.non_symmetric_error / sym_res.entire_mass);
#else
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
#endif  // USE_MP
    }
    std::cout << std::flush;
  }
}
