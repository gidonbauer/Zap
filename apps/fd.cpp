#include <filesystem>
#include <numbers>

#include <AD/ad.hpp>

#define ZAP_NO_CHUNK_INFO
#include "CellBased/Solver.hpp"

#include "IO/NoopWriter.hpp"
#include "IO/ToNpy.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Macros.hpp"
#include "Igor/Timer.hpp"

#include "Common.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "fd/"

using namespace Zap::CellBased;
using namespace Zap::IO;
using namespace std::string_view_literals;

// - Setup -----------------------------------------------------------------------------------------
using PassiveFloat           = double;
using ActiveFloat            = ad::gt1s<PassiveFloat>::type;
constexpr PassiveFloat X_MIN = 0.0;
constexpr PassiveFloat X_MAX = 5.0;
constexpr PassiveFloat Y_MIN = 0.0;
constexpr PassiveFloat Y_MAX = 5.0;

// - Available command line options ----------------------------------------------------------------
struct Args {
  size_t nx                      = 25;
  size_t ny                      = 25;
  PassiveFloat tend              = 1.0;
  PassiveFloat CFL_safety_factor = 0.25;
  bool first_order               = false;
};

// - Print usage -----------------------------------------------------------------------------------
template <Igor::detail::Level level>
void usage(std::string_view prog, std::ostream& out) noexcept {
  Args args{};

  out << Igor::detail::level_repr(level) << "Usage: " << prog
      << " [--nx nx] [--ny ny] [--tend tend] [--CFL CFL-safety-factor] [--first-order]\n";
  out << "\t--nx              Number of cells in x-direction, default is " << args.nx << '\n';
  out << "\t--ny              Number of cells in y-direction, default is " << args.ny << '\n';
  out << "\t--tend            Final time of simulation, default is " << args.tend << '\n';
  out << "\t--CFL             Safety factor for CFL condition, must be in (0, 1), default is "
      << args.CFL_safety_factor << '\n';
  out << "\t--first-order     Use first order finite differences, default is " << std::boolalpha
      << args.first_order << '\n';
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
    } else if (arg == "--tend"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --tend.\n";
        return std::nullopt;
      }
      args.tend = parse_double(*arg);
    } else if (arg == "--CFL"sv || arg == "--cfl"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --CFL.\n";
        return std::nullopt;
      }
      args.CFL_safety_factor = parse_double(*arg);
    } else if (arg == "--first-order"sv) {
      args.first_order = true;
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
    -> std::optional<UniformGrid<ActiveFloat, PassiveFloat>> {
  // Igor::ScopeTimer timer(std::format("Solver with nx={}, ny={}, eps={}", nx, ny, eps_value));

  ActiveFloat eps     = eps_value;
  ad::derivative(eps) = 1;

  constexpr PassiveFloat r = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;
  auto u0                  = [&](ActiveFloat x, ActiveFloat y) -> ActiveFloat {
    return (1 + eps) * (std::pow(x - X_MIN, 2) + std::pow(y - Y_MIN, 2)) *
           static_cast<ActiveFloat>((std::pow(x - X_MIN, 2) + std::pow(y - Y_MIN, 2)) <=
                                    std::pow(r, 2));
  };

  UniformGrid<ActiveFloat, PassiveFloat> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.periodic_boundary();

  auto init_shock = [=]<typename T>(T t) -> SimCoord<T> {
    IGOR_ASSERT(t >= 0 && t <= 1, "t must be in [0, 1] but is {}", t);
    return {
        r * std::cos(std::numbers::pi_v<PassiveFloat> / 2 * t),
        r * std::sin(std::numbers::pi_v<PassiveFloat> / 2 * t),
    };
  };

  if (!grid.cut_curve(init_shock)) { return std::nullopt; }
  grid.fill_quad(u0);

  NoopWriter grid_writer;
  NoopWriter t_writer;
  return solve_2d_burgers<ExtendType::NEAREST>(
      grid, tend, grid_writer, t_writer, CFL_safety_factor);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto run_passive(size_t nx,
                               size_t ny,
                               PassiveFloat tend,
                               PassiveFloat eps_value,
                               PassiveFloat CFL_safety_factor) noexcept
    -> std::optional<UniformGrid<PassiveFloat, PassiveFloat>> {
  // Igor::ScopeTimer timer(std::format("Solver with nx={}, ny={}, eps={}", nx, ny, eps_value));

  constexpr PassiveFloat r = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;
  auto u0                  = [=](PassiveFloat x, PassiveFloat y) -> PassiveFloat {
    return (1 + eps_value) * (std::pow(x - X_MIN, 2) + std::pow(y - Y_MIN, 2)) *
           static_cast<PassiveFloat>((std::pow(x - X_MIN, 2) + std::pow(y - Y_MIN, 2)) <=
                                     std::pow(r, 2));
  };

  UniformGrid<PassiveFloat, PassiveFloat> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.periodic_boundary();

  auto init_shock = [=]<typename T>(T t) -> SimCoord<T> {
    IGOR_ASSERT(t >= 0 && t <= 1, "t must be in [0, 1] but is {}", t);
    return {
        r * std::cos(std::numbers::pi_v<PassiveFloat> / 2 * t),
        r * std::sin(std::numbers::pi_v<PassiveFloat> / 2 * t),
    };
  };

  if (!grid.cut_curve(init_shock)) { return std::nullopt; }
  grid.fill_quad(u0);

  NoopWriter grid_writer;
  NoopWriter t_writer;
  return solve_2d_burgers<ExtendType::NEAREST>(
      grid, tend, grid_writer, t_writer, CFL_safety_factor);
}

// -------------------------------------------------------------------------------------------------
auto main(int argc, char** argv) -> int {
  const auto args = parse_args(argc, argv);
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

  const PassiveFloat eps = 1e-8;

  Igor::Info("nx   = {}", args->nx);
  Igor::Info("ny   = {}", args->ny);
  Igor::Info("tend = {}", args->tend);
  Igor::Info("CFL  = {}", args->CFL_safety_factor);
  Igor::Info("first-order FD = {}", args->first_order);
  std::cout << "------------------------------------------------------------\n";

  const auto res_ad = run(args->nx, args->ny, args->tend, PassiveFloat{0}, args->CFL_safety_factor);
  if (!res_ad.has_value()) {
    Igor::Warn("Solver for eps=0 failed.");
    return 1;
  }

  const auto shock_curve_ad = res_ad->get_shock_curve();
  std::vector<SimCoord<PassiveFloat>> shock_ad(shock_curve_ad.size());
  std::vector<SimCoord<PassiveFloat>> der_shock_ad(shock_curve_ad.size());
  for (size_t i = 0; i < shock_curve_ad.size(); ++i) {
    shock_ad[i].x     = ad::value(shock_curve_ad[i].x);
    shock_ad[i].y     = ad::value(shock_curve_ad[i].y);
    der_shock_ad[i].x = ad::derivative(shock_curve_ad[i].x);
    der_shock_ad[i].y = ad::derivative(shock_curve_ad[i].y);
  }

  const auto res_plus_eps =
      run_passive(args->nx, args->ny, args->tend, eps, args->CFL_safety_factor);
  if (!res_plus_eps.has_value()) {
    Igor::Warn("Solver for eps={} failed.", eps);
    return 1;
  }
  const auto res_minus_eps =
      run_passive(args->nx, args->ny, args->tend, -eps, args->CFL_safety_factor);
  if (!res_minus_eps.has_value()) {
    Igor::Warn("Solver for eps={} failed.", -eps);
    return 1;
  }

  const auto shock_plus_eps  = res_plus_eps->get_shock_curve();
  const auto shock_minus_eps = res_minus_eps->get_shock_curve();

  if (shock_plus_eps.size() != shock_minus_eps.size() || shock_plus_eps.size() != shock_ad.size()) {
    Igor::Warn("Incompatible shock curve sizes: eps=0 => {}, eps={} => {}, eps={} => {}",
               shock_ad.size(),
               eps,
               shock_plus_eps.size(),
               -eps,
               shock_minus_eps.size());
  }
  const auto ns = shock_ad.size();

  std::vector<SimCoord<PassiveFloat>> der_shock_fd(ns);
  for (size_t i = 0; i < ns; ++i) {
    if (args->first_order) {
      der_shock_fd[i].x = (shock_plus_eps[i].x - shock_ad[i].x) / eps;
      der_shock_fd[i].y = (shock_plus_eps[i].y - shock_ad[i].y) / eps;
    } else {
      der_shock_fd[i].x = (shock_plus_eps[i].x - shock_minus_eps[i].x) / (2 * eps);
      der_shock_fd[i].y = (shock_plus_eps[i].y - shock_minus_eps[i].y) / (2 * eps);
    }
  }

  // Save shock curve for eps=0
  {
    constexpr auto filename = OUTPUT_DIR "shock_curve.npy";
    if (!points_to_npy(filename, shock_ad)) {
      Igor::Warn("Could not save points to `{}`", filename);
      return 1;
    }
    Igor::Info("Saved shock curve to `{}`.", filename);
  }

  // Save shock curve for eps=eps
  {
    constexpr auto filename = OUTPUT_DIR "shock_curve_peps.npy";
    if (!points_to_npy(filename, shock_plus_eps)) {
      Igor::Warn("Could not save points to `{}`", filename);
      return 1;
    }
    Igor::Info("Saved shock curve for eps={} to `{}`.", eps, filename);
  }

  // Save shock curve for eps=-eps
  {
    constexpr auto filename = OUTPUT_DIR "shock_curve_meps.npy";
    if (!points_to_npy(filename, shock_minus_eps)) {
      Igor::Warn("Could not save points to `{}`", filename);
      return 1;
    }
    Igor::Info("Saved shock curve for eps={} to `{}`.", -eps, filename);
  }

  // Save derivative calculated via AD
  {
    constexpr auto filename = OUTPUT_DIR "der_ad.npy";
    if (!points_to_npy(filename, der_shock_ad)) {
      Igor::Warn("Could not save points to `{}`", filename);
      return 1;
    }
    Igor::Info("Saved derivative calculated via AD to `{}`.", filename);
  }

  // Save derivative calculated via FD
  {
    constexpr auto filename = OUTPUT_DIR "der_fd.npy";
    if (!points_to_npy(filename, der_shock_fd)) {
      Igor::Warn("Could not save points to `{}`", filename);
      return 1;
    }
    Igor::Info("Saved derivative calculated via FD to `{}`.", filename);
  }
}
