#include <filesystem>
#include <numbers>

#include <AD/ad.hpp>

// #define ZAP_2ND_ORDER_CORRECTION
// #define ZAP_TANGENTIAL_CORRECTION
// #define ZAP_STATIC_CUT

#include "CellBased/ReconstructShock.hpp"
#include "CellBased/Solver.hpp"
#include "IO/IncCellWriter.hpp"
#include "IO/IncMatrixWriter.hpp"
#include "IO/VTKWriter.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Macros.hpp"
#include "Igor/Timer.hpp"

#include "Common.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "cell_based/"

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
  size_t nx                          = 25;
  size_t ny                          = 25;
  PassiveFloat tend                  = 1.0;
  PassiveFloat eps                   = 0.0;
  PassiveFloat CFL_safety_factor     = 0.25;
  std::string_view print_shock_curve = ""sv;
  size_t ns                          = 20;
  bool print_cuts                    = false;
};

// - Print usage -----------------------------------------------------------------------------------
template <Igor::detail::Level level>
void usage(std::string_view prog, std::ostream& out) noexcept {
  Args args{};

  out << Igor::detail::level_repr(level) << "Usage: " << prog
      << " [--nx nx] [--ny ny] [--tend tend] [--eps eps] [--CFL CFL-safety-factor] [--print-shock "
         "method] "
         "[--ns ns] [--print-cuts]\n";
  out << "\t--nx              Number of cells in x-direction, default is " << args.nx << '\n';
  out << "\t--ny              Number of cells in y-direction, default is " << args.ny << '\n';
  out << "\t--tend            Final time for simulation, default is " << args.tend << '\n';
  out << "\t--eps             Perturbation of initial condition, default is " << args.eps << '\n';
  out << "\t--CFL             Safety factor for CFL condition, must be in (0, 1), default is "
      << args.CFL_safety_factor << '\n';
  out << "\t--print-shock     Print the final shock curve and its derivative. Decide which "
         "reconstruction is used. Available options are 'ALL', 'LINEAR', 'CUBIC', and 'SMOOTHSTEP' "
         "default is "
      << (args.print_shock_curve.empty() ? "not printing the shock" : args.print_shock_curve)
      << '\n';
  out << "\t--ns              Number of points at which the reconstruction of the shock is "
         "evaluated, default is "
      << args.ns << '\n';
  out << "\t--print-cuts      Print the final cut values and their derivative, default is "
      << std::boolalpha << args.print_cuts << '\n';
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
    } else if (arg == "--eps"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --eps.\n";
        return std::nullopt;
      }
      args.eps = parse_double(*arg);
    } else if (arg == "--CFL"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --CFL.\n";
        return std::nullopt;
      }
      args.CFL_safety_factor = parse_double(*arg);
    } else if (arg == "--print-shock"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide method for option --print-shock.\n";
        return std::nullopt;
      }
      args.print_shock_curve = *arg;
    } else if (arg == "--print-cuts"sv) {
      args.print_cuts = true;
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
auto main(int argc, char** argv) -> int {
  const std::string_view prog = *argv;
  const auto args             = parse_args(argc, argv);
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
  Igor::Info("eps  = {}", args->eps);
  Igor::Info("CFL  = {}", args->CFL_safety_factor);

  Zap::CellBased::UniformGrid<ActiveFloat, PassiveFloat> grid(
      X_MIN, X_MAX, args->nx, Y_MIN, Y_MAX, args->ny);
  // grid.same_value_boundary();
  grid.periodic_boundary();

  ActiveFloat eps     = args->eps;
  ad::derivative(eps) = 1;

  // #define RAMP_X
#define QUARTER_CIRCLE
#ifdef QUARTER_CIRCLE
  const PassiveFloat r = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;

  auto u0 = [=](ActiveFloat x, ActiveFloat y) -> ActiveFloat {
    return (1 + eps) * (std::pow(x - X_MIN, 2) + std::pow(y - Y_MIN, 2)) *
           static_cast<ActiveFloat>((std::pow(x - X_MIN, 2) + std::pow(y - Y_MIN, 2)) <=
                                    std::pow(r, 2));
  };

  [[maybe_unused]] auto init_shock = [=]<typename T>(T t) -> Zap::CellBased::SimCoord<T> {
    // assert(t >= 0 && t <= 1);
    return {
        r * std::cos(std::numbers::pi_v<PassiveFloat> / 2 * t),
        r * std::sin(std::numbers::pi_v<PassiveFloat> / 2 * t),
    };
  };
#elif defined(RAMP_X)
  auto u0 = [=](ActiveFloat x, ActiveFloat /*y*/) -> ActiveFloat {
    // return (1+eps) * (x - x_min) * static_cast<ActiveFloat>((x - X_MIN) < (X_MAX - X_MIN) / 2);

    // TODO: Inverse ramp is not solved correctly
    return (X_MAX - x) * (1.0 - static_cast<ActiveFloat>((x - X_MIN) < (X_MAX - X_MIN) / 2));
  };

  [[maybe_unused]] auto init_shock = [=](PassiveFloat t) -> Zap::CellBased::SimCoord<PassiveFloat> {
    assert(t >= 0 && t <= 1);
    return {
        (X_MAX - X_MIN) / 2,
        t * (Y_MAX - Y_MIN) + Y_MIN,
    };
  };
#elif defined(HAT_X)
  auto u0 = [=](Float x, Float /*y*/) -> Float {
    if ((x - X_MIN) < 0.25 * (X_MAX - X_MIN)) {
      return x;
    } else if ((x - X_MIN) < 0.5 * (X_MAX - X_MIN)) {
      return 0.25 * (X_MAX - X_MIN) - (x - 0.25 * (X_MAX - X_MIN));
    } else {
      return 0;
    }
  };

  [[maybe_unused]] auto init_shock = [=]<typename T>(T t) -> Zap::CellBased::SimCoord<T> {
    static_assert(false, "Not implemented yet.");
    return Zap::CellBased::SimCoord<T>::Zero();
  };
#else
  static_assert(false, "No initial condition defined.");
#endif

  IGOR_TIME_SCOPE("Cutting the grid") {
    if (!grid.cut_curve(init_shock)) { return 1; }
  }
  grid.fill_four_point(u0);

  constexpr auto u_file = OUTPUT_DIR "u_1d.grid";
  Zap::IO::IncCellWriter<ActiveFloat, PassiveFloat> grid_writer{u_file, grid};

  constexpr auto t_file = OUTPUT_DIR "t_1d.mat";
  Zap::IO::IncMatrixWriter<PassiveFloat, 1, 1, 0> t_writer(t_file, 1, 1, 0);

// #define SAVE_ONLY_INITIAL_STATE
#ifdef SAVE_ONLY_INITIAL_STATE
  if (!grid_writer.write_data(grid)) { return 1; }
  if (!t_writer.write_data(PassiveFloat{0})) { return 1; }
#else
  IGOR_TIME_SCOPE("Solver") {
    Zap::CellBased::Solver<Zap::CellBased::ExtendType::NEAREST> solver;
    const auto res = solver.solve(grid, args->tend, grid_writer, t_writer, args->CFL_safety_factor);
    if (!res.has_value()) {
      Igor::Warn("Solver failed.");
      return 1;
    }

    if (!args->print_shock_curve.empty()) {
      Igor::Info("ns = {}", args->ns);
      Igor::Info("Reconsturction method = {}", args->print_shock_curve);

      const auto final_shock = [&] {
        if (args->print_shock_curve == "ALL"sv) { return res->get_shock_curve(); }

        std::vector<PassiveFloat> ts(args->ns);
        std::generate(std::begin(ts), std::end(ts), [i = 0UZ, &args]() mutable {
          return static_cast<PassiveFloat>(i++) / static_cast<PassiveFloat>(args->ns - 1);
        });
        if (args->print_shock_curve == "LINEAR"sv) {
          return Zap::CellBased::piecewise_linear(ts, res->get_shock_curve());
        }
        if (args->print_shock_curve == "CUBIC"sv) {
          return Zap::CellBased::natural_cubic_spline(ts, res->get_shock_curve());
        }
        if (args->print_shock_curve == "SMOOTHSTEP"sv) {
          return Zap::CellBased::smoothstep(ts, res->get_shock_curve());
        }

        usage<Igor::detail::Level::PANIC>(prog, std::cerr);
        Igor::Panic("Invaild method `{}` to reconstruction the shock", args->print_shock_curve);
        std::unreachable();
      }();

      std::vector<Zap::CellBased::SimCoord<PassiveFloat>> x_shock(final_shock.size());
      std::transform(std::cbegin(final_shock),
                     std::cend(final_shock),
                     std::begin(x_shock),
                     [](const Zap::CellBased::SimCoord<ActiveFloat>& p) {
                       return Zap::CellBased::SimCoord<PassiveFloat>{ad::value(p.x),
                                                                     ad::value(p.y)};
                     });

      std::vector<Zap::CellBased::SimCoord<PassiveFloat>> xi_shock(final_shock.size());
      std::transform(std::cbegin(final_shock),
                     std::cend(final_shock),
                     std::begin(xi_shock),
                     [](const Zap::CellBased::SimCoord<ActiveFloat>& p) {
                       return Zap::CellBased::SimCoord<PassiveFloat>{ad::derivative(p.x),
                                                                     ad::derivative(p.y)};
                     });

      std::cout << "==============================================================================="
                   "==========\n";
      Igor::Info("x_shock = {}", x_shock);
      std::cout << "==============================================================================="
                   "==========\n";
      Igor::Info("xi_shock = {}", xi_shock);
      std::cout << "==============================================================================="
                   "==========\n";
    }

    if (args->print_cuts) {
      const auto final_cuts = res->get_cuts();

      std::vector<Zap::CellBased::SimCoord<PassiveFloat>> cuts(final_cuts.size());
      std::transform(std::cbegin(final_cuts),
                     std::cend(final_cuts),
                     std::begin(cuts),
                     [](const Zap::CellBased::SimCoord<ActiveFloat>& p) {
                       return Zap::CellBased::SimCoord<PassiveFloat>{ad::value(p.x),
                                                                     ad::value(p.y)};
                     });

      std::vector<Zap::CellBased::SimCoord<PassiveFloat>> dcuts(final_cuts.size());
      std::transform(std::cbegin(final_cuts),
                     std::cend(final_cuts),
                     std::begin(dcuts),
                     [](const Zap::CellBased::SimCoord<ActiveFloat>& p) {
                       return Zap::CellBased::SimCoord<PassiveFloat>{ad::derivative(p.x),
                                                                     ad::derivative(p.y)};
                     });

      std::cout << "==============================================================================="
                   "==========\n";
      Igor::Info("cuts = {}", cuts);
      std::cout << "==============================================================================="
                   "==========\n";
      Igor::Info("dcuts = {}", dcuts);
      std::cout << "==============================================================================="
                   "==========\n";
    }
  }
  Igor::Info("Solver finished successfully.");
#endif
  Igor::Info("Saved grid to {}.", u_file);
  Igor::Info("Saved time steps to {}.", t_file);
}
