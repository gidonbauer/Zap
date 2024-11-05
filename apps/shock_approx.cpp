#include <filesystem>
#include <numbers>

#include <AD/ad.hpp>

#define ZAP_TANGENTIAL_CORRECTION
// #define ZAP_STATIC_CUT

// #define X_RAMP

#include "CellBased/EigenDecomp.hpp"
#include "CellBased/ReconstructShock.hpp"
#include "CellBased/Solver.hpp"

#include "IO/ToNpy.hpp"
#ifdef SAVE_SOLUTION
#include "IO/IncCellWriter.hpp"
#include "IO/IncMatrixWriter.hpp"
#else
#include "IO/NoopWriter.hpp"
#endif  // SAVE_SOLUTION

#include "Igor/Logging.hpp"
#include "Igor/Macros.hpp"
#include "Igor/Timer.hpp"

#include "Common.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "shock_approx/"

using namespace Zap::CellBased;
using namespace std::string_view_literals;

// - Setup -----------------------------------------------------------------------------------------
using PassiveFloat           = double;
using ActiveFloat            = ad::gt1s<PassiveFloat>::type;
constexpr size_t DIM         = 1;
constexpr PassiveFloat X_MIN = 0.0;
constexpr PassiveFloat X_MAX = 2.0;
constexpr PassiveFloat Y_MIN = 0.0;
constexpr PassiveFloat Y_MAX = 2.0;

std::string_view g_executable_name;  // NOLINT

// - Available command line options ----------------------------------------------------------------
struct Args {
  size_t nx                             = 25;
  size_t ny                             = 25;
  PassiveFloat tend                     = 1.0;
  PassiveFloat CFL_safety_factor        = 0.25;
  std::string_view shock_reconstruction = "LINEAR"sv;
  size_t ns                             = 20;
};

// - Print usage -----------------------------------------------------------------------------------
template <Igor::detail::Level level>
void usage(std::string_view prog, std::ostream& out) noexcept {
  Args args{};

  out << Igor::detail::level_repr(level) << "Usage: " << prog
      << " [--nx nx] [--ny ny] [--tend tend] [--CFL CFL-safety-factor] [--method "
         "method] "
         "[--ns ns]\n";
  out << "\t--nx              Number of cells in x-direction, default is " << args.nx << '\n';
  out << "\t--ny              Number of cells in y-direction, default is " << args.ny << '\n';
  out << "\t--tend            Final time for simulation, default is " << args.tend << '\n';
  out << "\t--CFL             Safety factor for CFL condition, must be in (0, 1), default is "
      << args.CFL_safety_factor << '\n';
  out << "\t--method          Method used to reconstruct the shock, available options are 'ALL', "
         "'LINEAR', 'CUBIC', and 'SMOOTHSTEP', default is "
      << args.shock_reconstruction << '\n';
  out << "\t--ns              Number of points at which the reconstruction of the shock is "
         "evaluated, default is "
      << args.ns << '\n';
}

// - Parse command line arguments ------------------------------------------------------------------
[[nodiscard]] auto parse_args(int argc, char** argv) noexcept -> std::optional<Args> {
  Args args{};

  const auto prog = next_arg(argc, argv);
  IGOR_ASSERT(prog.has_value(), "Could not get program name from command line arguments.");

  for (auto arg = next_arg(argc, argv); arg.has_value(); arg = next_arg(argc, argv)) {
    using namespace std::string_view_literals;
    if (arg == "--nx"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --nx.\n";
        return std::nullopt;
      }
      args.nx = parse_size_t(*arg);
    } else if (arg == "--ny"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --ny.\n";
        return std::nullopt;
      }
      args.ny = parse_size_t(*arg);
    } else if (arg == "--ns"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --ns.\n";
        return std::nullopt;
      }
      args.ns = parse_size_t(*arg);
    } else if (arg == "--tend"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --tend.\n";
        return std::nullopt;
      }
      args.tend = parse_double(*arg);
    } else if (arg == "--CFL"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --CFL.\n";
        return std::nullopt;
      }
      args.CFL_safety_factor = parse_double(*arg);
    } else if (arg == "--method"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide method for option --method.\n";
        return std::nullopt;
      }
      args.shock_reconstruction = *arg;
    } else if (arg == "-h"sv || arg == "--help"sv) {
      usage<Igor::detail::Level::INFO>(*prog, std::cout);
      return std::nullopt;
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
[[nodiscard]] auto run(size_t nx,
                       size_t ny,
                       PassiveFloat tend,
                       PassiveFloat eps_value,
                       PassiveFloat CFL_safety_factor) noexcept
    -> std::optional<UniformGrid<ActiveFloat, PassiveFloat, DIM>> {
  Igor::ScopeTimer timer(std::format("Solver with nx={}, ny={}, eps={}", nx, ny, eps_value));

  ActiveFloat eps     = eps_value;
  ad::derivative(eps) = 1;

#ifndef X_RAMP
  constexpr PassiveFloat r = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;
  auto u0                  = [=](ActiveFloat x, ActiveFloat y) -> ActiveFloat {
    static_assert(DIM == 1);
    return (1 + eps) * (std::pow(x - X_MIN, 2) + std::pow(y - Y_MIN, 2)) *
           static_cast<ActiveFloat>((std::pow(x - X_MIN, 2) + std::pow(y - Y_MIN, 2)) <=
                                    std::pow(r, 2));
  };

  auto init_shock = [=]<typename T>(T t) -> SimCoord<T> {
    // assert(t >= 0 && t <= 1);
    return {
        r * std::cos(std::numbers::pi_v<PassiveFloat> / 2 * t),
        r * std::sin(std::numbers::pi_v<PassiveFloat> / 2 * t),
    };
  };
#else
  constexpr PassiveFloat r = (X_MIN + X_MAX) / 2;
  auto u0                  = [&](ActiveFloat x, ActiveFloat /*y*/) -> ActiveFloat {
    static_assert(DIM == 1);
    return (1 + eps) * x * static_cast<PassiveFloat>(x < r);
  };
#endif  // X_RAMP

  UniformGrid<ActiveFloat, PassiveFloat, DIM> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.periodic_boundary();

#ifndef X_RAMP
  if (!grid.cut_curve(init_shock)) { return std::nullopt; }
#else
  if (std::vector<SimCoord<PassiveFloat>> points = {{r, Y_MIN}, {r, Y_MAX}};
      !grid.cut_piecewise_linear<ExtendType::MAX>(points)) {
    return std::nullopt;
  }
#endif  // X_RAMP
  grid.fill_four_point(u0);

#ifdef SAVE_SOLUTION
  const std::string u_file = OUTPUT_DIR "u_" + std::to_string(nx) + "x" + std::to_string(ny) + "_" +
                             std::to_string(eps_value) + ".grid";
  Zap::IO::IncCellWriter<ActiveFloat, PassiveFloat, DIM> grid_writer{u_file, grid};

  const std::string t_file = OUTPUT_DIR "t_" + std::to_string(nx) + "x" + std::to_string(ny) + "_" +
                             std::to_string(eps_value) + ".mat";
  Zap::IO::IncMatrixWriter<PassiveFloat, 1, 1, 0> t_writer(t_file, 1, 1, 0);
#else
  Zap::IO::NoopWriter grid_writer;
  Zap::IO::NoopWriter t_writer;
#endif  // SAVE_SOLUTION

  auto solver = make_solver<ExtendType::NEAREST>(SingleEq::A{}, SingleEq::B{});
  return solver.solve(grid, tend, grid_writer, t_writer, CFL_safety_factor);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto reconstruct_shock(const UniformGrid<ActiveFloat, PassiveFloat, DIM>& res,
                                     size_t ns,
                                     std::string_view method) noexcept {
  if (method == "ALL"sv) { return res.get_shock_curve(); }

  std::vector<PassiveFloat> ts(ns);
  std::generate(std::begin(ts), std::end(ts), [i = 0UZ, ns]() mutable {
    return static_cast<PassiveFloat>(i++) / static_cast<PassiveFloat>(ns - 1);
  });
  if (method == "LINEAR"sv) { return piecewise_linear(ts, res.get_shock_curve()); }
  if (method == "CUBIC"sv) { return natural_cubic_spline(ts, res.get_shock_curve()); }
  if (method == "SMOOTHSTEP"sv) { return smoothstep(ts, res.get_shock_curve()); }

  usage<Igor::detail::Level::PANIC>(g_executable_name, std::cerr);
  Igor::Panic("Invaild method `{}` to reconstruction the shock", method);
  std::unreachable();
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto compare(size_t nx,
                           size_t ny,
                           PassiveFloat tend,
                           PassiveFloat eps_value,
                           PassiveFloat CFL_safety_factor,
                           size_t ns,
                           std::string_view method,
                           const std::vector<SimCoord<ActiveFloat>>& eps_zero_shock) noexcept
    -> bool {
  const auto res = run(nx, ny, tend, eps_value, CFL_safety_factor);
  if (!res.has_value()) {
    Igor::Warn("Solver for eps={} failed.", eps_value);
    return false;
  }

  const auto eps_shock = reconstruct_shock(*res, ns, method);
  if (!Zap::IO::points_to_npy(OUTPUT_DIR "shock_" + std::to_string(eps_value) + ".npy",
                              eps_shock)) {
    return false;
  }

  std::vector<SimCoord<PassiveFloat>> eps_approx_shock(eps_shock.size());
  for (size_t i = 0; i < eps_approx_shock.size(); ++i) {
    eps_approx_shock[i] = {
        ad::value(eps_zero_shock[i].x) + eps_value * ad::derivative(eps_zero_shock[i].x),
        ad::value(eps_zero_shock[i].y) + eps_value * ad::derivative(eps_zero_shock[i].y),
    };
  }
  IGOR_ASSERT(eps_shock.size() == eps_approx_shock.size(),
              "eps_shock and eps_approx_shock must have the same size but the sizes are {} and {}",
              eps_shock.size(),
              eps_approx_shock.size());
  if (!Zap::IO::points_to_npy(OUTPUT_DIR "shock_approx_" + std::to_string(eps_value) + ".npy",
                              eps_approx_shock)) {
    return false;
  }

  IGOR_ASSERT(
      eps_shock.size() % 2 == 0, "ns has to be an even number for simpsons rule, ns is {}", ns);
  ActiveFloat L1_error = 0;
  const auto delta     = 1 / static_cast<PassiveFloat>(eps_shock.size());
  for (size_t i = 1; i <= eps_shock.size() / 2; ++i) {
    L1_error += (eps_shock[2 * i - 2] - eps_approx_shock[2 * i - 2]).norm() +
                4 * (eps_shock[2 * i - 1] - eps_approx_shock[2 * i - 1]).norm() +
                (eps_shock[2 * i - 0] - eps_approx_shock[2 * i - 0]).norm();
  }
  L1_error *= delta / 3;

  Igor::Info("{} => {}", eps_value, ad::value(L1_error));

  return true;
}

// -------------------------------------------------------------------------------------------------
auto main(int argc, char** argv) -> int {
  g_executable_name = *argv;
  const auto args   = parse_args(argc, argv);
  if (!args.has_value()) { return 1; }
  assert(args->nx < std::numeric_limits<int>::max());
  assert(args->ny < std::numeric_limits<int>::max());

  {
    std::error_code ec;
    std::filesystem::create_directories(OUTPUT_DIR, ec);
    if (ec) {
      Igor::Warn("Could not create directory `{}`: {}", OUTPUT_DIR, ec.message());
      return 1;
    }
  }

  if (args->tend <= 0.0) {
    Igor::Warn("tend must be larger than 0, but is {}.", args->tend);
    return 1;
  }

  Igor::Info("nx                    = {}", args->nx);
  Igor::Info("ny                    = {}", args->ny);
  Igor::Info("tend                  = {}", args->tend);
  Igor::Info("CFL                   = {}", args->CFL_safety_factor);

#ifdef RUN_COMPARE
  Igor::Info("Reconsturction method = {}", args->shock_reconstruction);
  Igor::Info("ns                    = {}", args->ns);

  const auto res = run(args->nx, args->ny, args->tend, PassiveFloat{0}, args->CFL_safety_factor);
  if (!res.has_value()) {
    Igor::Warn("Solver for eps=0 failed.");
    return 1;
  }
  const auto eps_zero_shock = reconstruct_shock(*res, args->ns, args->shock_reconstruction);
  if (!Zap::IO::points_to_npy(OUTPUT_DIR "shock_0.npy", eps_zero_shock)) { return 1; }
  std::cout << "------------------------------------------------------------\n";

  constexpr std::array eps_values = {2e-1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5};
  bool all_good                   = true;
  for (auto eps_value : eps_values) {
    std::cout << '\n';
    const bool success = compare(args->nx,
                                 args->ny,
                                 args->tend,
                                 eps_value,
                                 args->CFL_safety_factor,
                                 args->ns,
                                 args->shock_reconstruction,
                                 eps_zero_shock);

    all_good = success && all_good;
  }

  return all_good ? 0 : 1;
#else
  constexpr std::array eps_values = {
      0.0,
      1e-5,
      1e-4,
      1e-3,
      2e-3,
      5e-3,
      1e-2,
      2e-2,
      5e-2,
      1e-1,
      2e-1,
      5e-1,
  };
  bool failed = false;
  for (auto eps_value : eps_values) {
    std::cout << '\n';
    const auto res = run(args->nx, args->ny, args->tend, eps_value, args->CFL_safety_factor);
    if (!res.has_value()) {
      Igor::Warn("Solver for eps={} failed.", eps_value);
      failed = true;
    } else {
      Igor::Info("Solver for eps={} finished successfully.", eps_value);
    }

    const auto cut_points = res->get_shock_curve();
    std::vector<SimCoord<PassiveFloat>> x_shock(cut_points.size());
    std::transform(std::cbegin(cut_points),
                   std::cend(cut_points),
                   std::begin(x_shock),
                   [](const SimCoord<ActiveFloat>& p) {
                     return SimCoord<PassiveFloat>{ad::value(p.x), ad::value(p.y)};
                   });

    std::vector<SimCoord<PassiveFloat>> xi_shock(cut_points.size());
    std::transform(std::cbegin(cut_points),
                   std::cend(cut_points),
                   std::begin(xi_shock),
                   [](const SimCoord<ActiveFloat>& p) {
                     return SimCoord<PassiveFloat>{ad::derivative(p.x), ad::derivative(p.y)};
                   });

#ifndef X_RAMP
    const std::string x_filename = OUTPUT_DIR "x_" + std::to_string(eps_value) + ".npy";
    if (!Zap::IO::points_to_npy(x_filename, x_shock)) {
      failed = true;
    } else {
      Igor::Info("Saved cut points to `{}`.", x_filename);
    }

    const std::string xi_filename = OUTPUT_DIR "xi_" + std::to_string(eps_value) + ".npy";
    if (!Zap::IO::points_to_npy(xi_filename, xi_shock)) {
      failed = true;
    } else {
      Igor::Info("Saved derivative of cut points to `{}`.", xi_filename);
    }
#else

    const auto avg_x_shock =
        std::transform_reduce(std::cbegin(x_shock),
                              std::cend(x_shock),
                              PassiveFloat{0},
                              std::plus<>{},
                              [](const SimCoord<PassiveFloat>& p) { return p.x; }) /
        static_cast<PassiveFloat>(x_shock.size());

    const auto avg_xi_shock =
        std::transform_reduce(std::cbegin(xi_shock),
                              std::cend(xi_shock),
                              PassiveFloat{0},
                              std::plus<>{},
                              [](const SimCoord<PassiveFloat>& p) { return p.x; }) /
        static_cast<PassiveFloat>(xi_shock.size());

    Igor::Info("eps = {}", eps_value);
    Igor::Info("avg_x_shock = {}", avg_x_shock);
    Igor::Info("avg_xi_shock = {}", avg_xi_shock);

#endif  // X_RAMP
  }

  return failed ? 1 : 0;

#endif  // RUN_COMPARE
}
