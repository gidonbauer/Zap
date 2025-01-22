#include <filesystem>
#include <numbers>

#include <AD/ad.hpp>

#define ZAP_NO_CHUNK_INFO
#include "CellBased/PredictShock.hpp"
#include "CellBased/Solver.hpp"

#include "IO/NoopWriter.hpp"
#include "IO/ToNpy.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Macros.hpp"
#include "Igor/Timer.hpp"

#include "Common.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "predict_shock/"

using namespace Zap::CellBased;
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
};

// - Print usage -----------------------------------------------------------------------------------
template <Igor::detail::Level level>
void usage(std::string_view prog, std::ostream& out) noexcept {
  Args args{};

  out << Igor::detail::level_repr(level) << "Usage: " << prog
      << " [--nx nx] [--ny ny] [--tend tend] [--CFL CFL-safety-factor]\n";
  out << "\t--nx              Number of cells in x-direction, default is " << args.nx << '\n';
  out << "\t--ny              Number of cells in y-direction, default is " << args.ny << '\n';
  out << "\t--tend            Final time for simulation, default is " << args.tend << '\n';
  out << "\t--CFL             Safety factor for CFL condition, must be in (0, 1), default is "
      << args.CFL_safety_factor << '\n';
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
    } else if (arg == "--CFL"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --CFL.\n";
        return std::nullopt;
      }
      args.CFL_safety_factor = parse_double(*arg);
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

  Zap::IO::NoopWriter grid_writer;
  Zap::IO::NoopWriter t_writer;
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

  if (args->tend <= 0.0) {
    Igor::Warn("tend must be larger than 0, but is {}.", args->tend);
    return 1;
  }

  Igor::Info("nx   = {}", args->nx);
  Igor::Info("ny   = {}", args->ny);
  Igor::Info("tend = {}", args->tend);
  Igor::Info("CFL  = {}", args->CFL_safety_factor);
  std::cout << "------------------------------------------------------------\n";

  const auto res_eps_zero =
      run(args->nx, args->ny, args->tend, PassiveFloat{0}, args->CFL_safety_factor);
  if (!res_eps_zero.has_value()) {
    Igor::Warn("Solver for eps=0 failed.");
    return 1;
  }
  const auto true_shock_zero = [&] {
    const auto active_shock = res_eps_zero->get_shock_curve();
    std::vector<SimCoord<PassiveFloat>> shock(active_shock.size());
    std::transform(std::cbegin(active_shock),
                   std::cend(active_shock),
                   std::begin(shock),
                   [](const SimCoord<ActiveFloat>& p) {
                     return SimCoord<PassiveFloat>{ad::value(p.x), ad::value(p.y)};
                   });
    return shock;
  }();

  {
    const auto filename = fmt::format("{}shock_{:.6f}.npy", OUTPUT_DIR, 0.0);
    if (!Zap::IO::points_to_npy(filename, res_eps_zero->get_shock_curve())) {
      Igor::Warn("Could not save shock to `{}`.", filename);
      return 1;
    }
  }

  const auto results_filename = fmt::format("{}shock_pred_results.txt", OUTPUT_DIR);
  std::ofstream res_out(results_filename);
  if (!res_out) {
    Igor::Warn("Could not open file `{}` for writing: {}", results_filename, std::strerror(errno));
    return 1;
  }

  // clang-format off
  constexpr std::array epss = {
      1e-5, 2e-5, 5e-5,
      1e-4, 2e-4, 5e-4,
      1e-3, 2e-3, 5e-3,
      1e-2, 2e-2, 5e-2,
      1e-1, 2e-1, 5e-1,
  };
  // clang-format on
  for (auto eps : epss) {
    const auto res_eps = run(args->nx, args->ny, args->tend, eps, args->CFL_safety_factor);
    if (!res_eps_zero.has_value()) {
      Igor::Warn("Solver for eps={} failed.", eps);
      return 1;
    }

    const auto pred_shock_eps = predict_shock_front(res_eps_zero->get_shock_curve(), eps);
    {
      const auto filename = fmt::format("{}pred_shock_{:.6f}.npy", OUTPUT_DIR, eps);
      if (!Zap::IO::points_to_npy(filename, pred_shock_eps)) {
        Igor::Warn("Could not save shock to `{}`.", filename);
        return 1;
      }
    }

    const auto true_shock_eps = [&] {
      const auto active_shock = res_eps->get_shock_curve();
      std::vector<SimCoord<PassiveFloat>> shock(active_shock.size());
      std::transform(std::cbegin(active_shock),
                     std::cend(active_shock),
                     std::begin(shock),
                     [](const SimCoord<ActiveFloat>& p) {
                       return SimCoord<PassiveFloat>{ad::value(p.x), ad::value(p.y)};
                     });
      return shock;
    }();
    {
      const auto filename = fmt::format("{}shock_{:.6f}.npy", OUTPUT_DIR, eps);
      if (!Zap::IO::points_to_npy(filename, true_shock_eps)) {
        Igor::Warn("Could not save shock to `{}`.", filename);
        return 1;
      }
    }

    const auto dif_no_pred = compare_curves<64UZ, false>(true_shock_eps, true_shock_zero, 10);
    const auto dif         = compare_curves<64UZ, false>(true_shock_eps, pred_shock_eps, 10);
    res_out << eps << " = " << dif << '\n';
    Igor::Info("{:.5f} => {:.12f} ({:.12f})", eps, dif, dif_no_pred);
  }
  std::cout << "------------------------------------------------------------\n";
  Igor::Info("Saved shock fronts in `{}`.", OUTPUT_DIR);
  Igor::Info("Saved results in `{}`.", results_filename);
}
