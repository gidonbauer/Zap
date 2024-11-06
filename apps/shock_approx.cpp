#include <filesystem>
#include <numbers>

#include <AD/ad.hpp>

#define ZAP_TANGENTIAL_CORRECTION
// #define ZAP_STATIC_CUT

#include "CellBased/EigenDecomp.hpp"
#include "CellBased/Solver.hpp"

#include "IO/IncCellWriter.hpp"
#include "IO/IncMatrixWriter.hpp"
#include "IO/NoopWriter.hpp"
#include "IO/ToNpy.hpp"

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
  size_t nx                      = 25;
  size_t ny                      = 25;
  PassiveFloat tend              = 1.0;
  PassiveFloat CFL_safety_factor = 0.25;
  PassiveFloat eps               = 0.0;
};

// - Print usage -----------------------------------------------------------------------------------
template <Igor::detail::Level level>
void usage(std::string_view prog, std::ostream& out) noexcept {
  Args args{};

  out << Igor::detail::level_repr(level) << "Usage: " << prog
      << " [--nx nx] [--ny ny] [--tend tend] [--CFL CFL-safety-factor] [--eps eps]\n";
  out << "\t--nx              Number of cells in x-direction, default is " << args.nx << '\n';
  out << "\t--ny              Number of cells in y-direction, default is " << args.ny << '\n';
  out << "\t--tend            Final time for simulation, default is " << args.tend << '\n';
  out << "\t--CFL             Safety factor for CFL condition, must be in (0, 1), default is "
      << args.CFL_safety_factor << '\n';
  out << "\t--eps             Perturbation of initial condition, default is " << args.eps << '\n';
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
    } else if (arg == "--eps"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --eps.\n";
        return std::nullopt;
      }
      args.eps = parse_double(*arg);
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

  UniformGrid<ActiveFloat, PassiveFloat, DIM> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.periodic_boundary();

  if (!grid.cut_curve(init_shock)) { return std::nullopt; }
  grid.fill_four_point(u0);

  Zap::IO::NoopWriter grid_writer;
  Zap::IO::NoopWriter t_writer;
  auto solver = make_solver<ExtendType::NEAREST>(SingleEq::A{}, SingleEq::B{});
  return solver.solve(grid, tend, grid_writer, t_writer, CFL_safety_factor);
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

  Igor::Info("nx   = {}", args->nx);
  Igor::Info("ny   = {}", args->ny);
  Igor::Info("tend = {}", args->tend);
  Igor::Info("CFL  = {}", args->CFL_safety_factor);
  Igor::Info("eps  = {}", args->eps);

  // 0.0,
  // 1e-5,
  // 1e-4,
  // 1e-3,
  // 2e-3,
  // 5e-3,
  // 1e-2,
  // 2e-2,
  // 5e-2,
  // 1e-1,
  // 2e-1,
  // 5e-1,
  const auto res = run(args->nx, args->ny, args->tend, args->eps, args->CFL_safety_factor);
  if (!res.has_value()) {
    Igor::Warn("Solver for eps={} failed.", args->eps);
    return 1;
  } else {
    Igor::Info("Solver for eps={} finished successfully.", args->eps);
  }

  //  - Save final solution ------------------------------------------------------------------------
  const std::string u_filename = fmt::format("{}u_{:.6f}.grid", OUTPUT_DIR, args->eps);
  const std::string v_filename = fmt::format("{}v_{:.6f}.grid", OUTPUT_DIR, args->eps);
  Zap::IO::IncCellWriter<ActiveFloat, PassiveFloat, DIM> u_writer(u_filename, v_filename, *res);
  if (!u_writer.write_data(*res)) { Igor::Warn("Could not save grid data."); }

  const std::string t_filename = fmt::format("{}t_{:.6f}.mat", OUTPUT_DIR, args->eps);
  Zap::IO::IncMatrixWriter<PassiveFloat, 1, 1, 0> t_writer(t_filename, 1, 1, 0);
  if (!t_writer.write_data(args->tend)) { Igor::Warn("Could not save time data."); }

  // - Save cut points -----------------------------------------------------------------------------
  const auto cut_points = res->get_shock_curve();
  std::vector<SimCoord<PassiveFloat>> x_shock(cut_points.size());
  std::transform(std::cbegin(cut_points),
                 std::cend(cut_points),
                 std::begin(x_shock),
                 [](const SimCoord<ActiveFloat>& p) {
                   return SimCoord<PassiveFloat>{ad::value(p.x), ad::value(p.y)};
                 });
  const std::string x_filename = fmt::format("{}x_{:.6f}.npy", OUTPUT_DIR, args->eps);
  if (!Zap::IO::points_to_npy(x_filename, x_shock)) {
    return 1;
  } else {
    Igor::Info("Saved cut points to `{}`.", x_filename);
  }

  // - Save cut derivatives ------------------------------------------------------------------------
  std::vector<SimCoord<PassiveFloat>> xi_shock(cut_points.size());
  std::transform(std::cbegin(cut_points),
                 std::cend(cut_points),
                 std::begin(xi_shock),
                 [](const SimCoord<ActiveFloat>& p) {
                   return SimCoord<PassiveFloat>{ad::derivative(p.x), ad::derivative(p.y)};
                 });
  const std::string xi_filename = fmt::format("{}xi_{:.6f}.npy", OUTPUT_DIR, args->eps);
  if (!Zap::IO::points_to_npy(xi_filename, xi_shock)) {
    return 1;
  } else {
    Igor::Info("Saved derivative of cut points to `{}`.", xi_filename);
  }
}
