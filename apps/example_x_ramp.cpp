#include <filesystem>
#include <iomanip>
#include <optional>

// #define ZAP_SERIAL
#define ZAP_NO_TANGENTIAL_CORRECTION

// TODO: Why does this only work if t is active?
// #define ZAP_T_ACTIVE

#include <AD/ad.hpp>

#include "CellBased/Solver.hpp"
#include "IO/IncCellWriter.hpp"
#include "IO/IncMatrixWriter.hpp"
#include "IO/ToNpy.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Macros.hpp"
#include "Igor/Timer.hpp"

#include "Common.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "example_x_ramp/"

using namespace Zap::CellBased;

// - Setup -----------------------------------------------------------------------------------------
using PassiveFloat           = double;
using ActiveFloat            = ad::gt1s<PassiveFloat>::type;
constexpr PassiveFloat X_MIN = 0.0;
constexpr PassiveFloat X_MAX = 2.0;
constexpr PassiveFloat Y_MIN = 0.0;
constexpr PassiveFloat Y_MAX = 2.0;

// - Available command line options ----------------------------------------------------------------
struct Args {
  bool run_benchmark    = false;
  bool square_benchmark = false;
  size_t nx             = 25;
  size_t ny             = 25;
  PassiveFloat tend     = 1.0;
  bool dump_npy         = false;
};

// - Print usage -----------------------------------------------------------------------------------
template <Igor::detail::Level level>
void usage(std::string_view prog, std::ostream& out) noexcept {
  Args args{};

  out << Igor::detail::level_repr(level) << "Usage: " << prog
      << " [--run-benchmark] [--square-benchmark] [--nx nx] [--ny ny] [--tend tend] [--dump-npy]\n";
  out << "\t--run-benchmark      Run the solver for different grid sizes and compare result "
         "against analytical solution, default is "
      << std::boolalpha << args.run_benchmark << '\n';
  out << "\t--square-benchmark   Run the benchmark with nxn grid instead of nx3 grid, default is "
      << std::boolalpha << args.square_benchmark << '\n';
  out << "\t--nx                 Number of cells in x-direction, default is " << args.nx << '\n';
  out << "\t--ny                 Number of cells in y-direction, default is " << args.ny << '\n';
  out << "\t--tend               Final time for simulation, default is " << args.tend << '\n';
  out << "\t--dump-npy           Write the averaged one-dimensional solution along x-axis to a "
         ".npy file, default is "
      << std::boolalpha << args.dump_npy << '\n';
}

// - Parse command line arguments ------------------------------------------------------------------
[[nodiscard]] auto parse_args(int& argc, char**& argv) noexcept -> std::optional<Args> {
  Args args{};

  const auto prog = next_arg(argc, argv);
  IGOR_ASSERT(prog.has_value(), "Could not get program name from command line arguments.");

  for (auto arg = next_arg(argc, argv); arg.has_value(); arg = next_arg(argc, argv)) {
    using namespace std::string_view_literals;
    if (arg == "--run-benchmark"sv || arg == "--run-benchmarks"sv) {
      args.run_benchmark = true;
    } else if (arg == "--square-benchmark"sv || arg == "--square-benchmarks"sv) {
      args.square_benchmark = true;
    } else if (arg == "--nx"sv) {
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
    } else if (arg == "--dump-npy") {
      args.dump_npy = true;
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
[[nodiscard]] constexpr auto chi(PassiveFloat x, PassiveFloat x_min, PassiveFloat x_max) noexcept
    -> PassiveFloat {
  return static_cast<PassiveFloat>(x >= x_min && x <= x_max);
}

// -------------------------------------------------------------------------------------------------
template <typename AT>
[[nodiscard]] constexpr auto x_shock(PassiveFloat t, AT eps) noexcept -> AT {
  return std::sqrt(1 + (1 + eps) * t);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto xi_shock(PassiveFloat t) noexcept -> PassiveFloat {
  return t / (2 * std::sqrt(1 + t));
}

// -------------------------------------------------------------------------------------------------
template <typename AT>
[[nodiscard]] constexpr auto u_eps(AT x, PassiveFloat t, AT eps) noexcept -> AT {
  return (1 + eps) * x / (1 + (1 + eps) * t) *
         chi(ad::value(x), PassiveFloat{0}, ad::value(x_shock(t, eps)));
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto v(PassiveFloat x, PassiveFloat t) noexcept -> PassiveFloat {
  return x / ((1 + t) * (1 + t)) * chi(x, 0, x_shock(t, PassiveFloat{0}));
}

// -------------------------------------------------------------------------------------------------
void print_shock_error(const UniformGrid<ActiveFloat, PassiveFloat>& numerical_solution,
                       PassiveFloat tend,
                       std::ostream& out) noexcept {
  const auto final_shock = numerical_solution.get_shock_curve();

  const auto avg_x_shock = std::transform_reduce(std::cbegin(final_shock),
                                                 std::cend(final_shock),
                                                 ActiveFloat{0},
                                                 std::plus<>{},
                                                 [](const auto& p) { return p.x; }) /
                           static_cast<PassiveFloat>(final_shock.size());

  const auto std_dev_x_shock = std::sqrt(
      std::transform_reduce(std::cbegin(final_shock),
                            std::cend(final_shock),
                            ActiveFloat{0},
                            std::plus<>{},
                            [=](const auto& p) { return std::pow(p.x - avg_x_shock, 2); }) /
      static_cast<PassiveFloat>(final_shock.size()));

  const auto avg_xi_shock =
      std::transform_reduce(std::cbegin(final_shock),
                            std::cend(final_shock),
                            ActiveFloat{0},
                            std::plus<>{},
                            [](const auto& p) { return ad::derivative(p.x); }) /
      static_cast<PassiveFloat>(final_shock.size());

  const auto std_dev_xi_shock =
      std::sqrt(std::transform_reduce(std::cbegin(final_shock),
                                      std::cend(final_shock),
                                      ActiveFloat{0},
                                      std::plus<>{},
                                      [=](const auto& p) {
                                        return std::pow(ad::derivative(p.x) - avg_xi_shock, 2);
                                      }) /
                static_cast<PassiveFloat>(final_shock.size()));

  const auto expected_x_shock = x_shock(tend, PassiveFloat{0});
  const auto abs_err_x_shock  = std::abs(avg_x_shock - expected_x_shock);
  const auto rel_err_x_shock  = abs_err_x_shock / expected_x_shock;

  const auto expected_xi_shock = xi_shock(tend);
  const auto abs_err_xi_shock  = std::abs(avg_xi_shock - expected_xi_shock);
  const auto rel_err_xi_shock  = abs_err_xi_shock / expected_xi_shock;

  out << "avg_x_shock       = " << std::setprecision(8) << avg_x_shock << '\n';
  out << "std_dev_x_shock   = " << std::setprecision(8) << std_dev_x_shock << '\n';
  out << "expected_x_shock  = " << std::setprecision(8) << expected_x_shock << '\n';
  out << "abs_err_x_shock   = " << std::setprecision(8) << abs_err_x_shock << '\n';
  out << "rel_err_x_shock   = " << std::setprecision(8) << rel_err_x_shock << '\n';
  out << "avg_xi_shock      = " << std::setprecision(8) << avg_xi_shock << '\n';
  out << "std_dev_xi_shock  = " << std::setprecision(8) << std_dev_xi_shock << '\n';
  out << "expected_xi_shock = " << std::setprecision(8) << expected_xi_shock << '\n';
  out << "abs_err_xi_shock  = " << std::setprecision(8) << abs_err_xi_shock << '\n';
  out << "rel_err_xi_shock  = " << std::setprecision(8) << rel_err_xi_shock << '\n';
}

// -------------------------------------------------------------------------------------------------
void print_solution_error(const UniformGrid<ActiveFloat, PassiveFloat>& numerical_solution,
                          PassiveFloat tend,
                          size_t n,
                          std::ostream& out) noexcept {

  auto L1_func_u = [&](PassiveFloat x) -> PassiveFloat {
    return std::abs(ad::value(numerical_solution.eval(SimCoord<PassiveFloat>{x, 1.0})) -
                    u_eps(x, tend, PassiveFloat{0}));
  };
  const auto L1_err_u = simpsons_rule_1d(L1_func_u, X_MIN, X_MAX, n);

  auto L1_func_v = [&](PassiveFloat x) -> PassiveFloat {
    return std::abs(ad::derivative(numerical_solution.eval(SimCoord<PassiveFloat>{x, 1.0})) -
                    v(x, tend));
  };
  const auto L1_err_v = simpsons_rule_1d(L1_func_v, X_MIN, X_MAX, n);

  out << "L1_err_u          = " << std::setprecision(8) << L1_err_u << '\n';
  out << "L1_err_v          = " << std::setprecision(8) << L1_err_v << '\n';
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto dump_to_npy(const UniformGrid<ActiveFloat, PassiveFloat>& grid) noexcept
    -> bool {
  const auto final_shock = grid.get_shock_curve();
  const auto avg_x_shock = std::transform_reduce(std::cbegin(final_shock),
                                                 std::cend(final_shock),
                                                 PassiveFloat{0},
                                                 std::plus<>{},
                                                 [](const auto& p) { return ad::value(p.x); }) /
                           static_cast<PassiveFloat>(final_shock.size());
  const auto avg_xi_shock =
      std::transform_reduce(std::cbegin(final_shock),
                            std::cend(final_shock),
                            PassiveFloat{0},
                            std::plus<>{},
                            [](const auto& p) { return ad::derivative(p.x); }) /
      static_cast<PassiveFloat>(final_shock.size());

  std::vector<PassiveFloat> xs;
  std::vector<PassiveFloat> us;
  std::vector<PassiveFloat> vs;

  for (size_t xi = 0; xi < grid.nx(); ++xi) {
    const size_t yi = grid.ny() / 2UZ;

    const auto& cell = grid[grid.to_vec_idx(xi, yi)];
    if (cell.is_cartesian()) {
      xs.push_back(cell.x_min<SIM_C>() + cell.dx<SIM_C>() / 2);
      us.push_back(ad::value(cell.get_cartesian().value));
      vs.push_back(ad::derivative(cell.get_cartesian().value));
    } else {
      xs.push_back((cell.x_min<SIM_C>() + avg_x_shock) / 2);
      us.push_back(ad::value(cell.get_cut().left_value));
      vs.push_back(ad::derivative(cell.get_cut().left_value));

      xs.push_back((cell.x_min<SIM_C>() + cell.dx<SIM_C>() + avg_x_shock) / 2);
      us.push_back(ad::value(cell.get_cut().right_value));
      vs.push_back(ad::derivative(cell.get_cut().right_value));
    }
  }

  constexpr auto x_file   = OUTPUT_DIR "x.npy";
  constexpr auto x1_file  = OUTPUT_DIR "x1.npy";
  constexpr auto xi1_file = OUTPUT_DIR "xi1.npy";
  constexpr auto u_file   = OUTPUT_DIR "u.npy";
  constexpr auto v_file   = OUTPUT_DIR "v.npy";

  if (!Zap::IO::vector_to_npy(x_file, xs)) { return false; }
  Igor::Info("Saved one-dimensional grid to {}", x_file);
  if (!Zap::IO::vector_to_npy(x1_file, std::vector{avg_x_shock})) { return false; }
  Igor::Info("Saved shock location to {}", x1_file);
  if (!Zap::IO::vector_to_npy(xi1_file, std::vector{avg_xi_shock})) { return false; }
  Igor::Info("Saved change of shock location to {}", xi1_file);
  if (!Zap::IO::vector_to_npy(u_file, us)) { return false; }
  Igor::Info("Saved one-dimensional solution to {}", u_file);
  if (!Zap::IO::vector_to_npy(v_file, vs)) { return false; }
  Igor::Info("Saved one-dimensional derivative to {}", v_file);

  return true;
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto
run(size_t nx, size_t ny, PassiveFloat tend, std::ostream& out, bool dump_npy) noexcept -> bool {
  UniformGrid<ActiveFloat, PassiveFloat> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.periodic_boundary();

  // if (!grid.cut_curve(init_shock)) { return false; }
  {
    std::vector<SimCoord<PassiveFloat>> points{
        {(X_MAX - X_MIN) / 2, 0.0},
        {(X_MAX - X_MIN) / 2, (Y_MAX - Y_MIN) + Y_MIN},
    };
    if (!grid.cut_piecewise_linear<ExtendType::MAX>(points)) { return false; }
  }

  ActiveFloat eps     = 0.0;
  ad::derivative(eps) = 1.0;

  auto u0 = [&eps](ActiveFloat x, ActiveFloat /*y*/) { return u_eps(x, PassiveFloat{0}, eps); };
  grid.fill_four_point(u0);

  const auto u_file = OUTPUT_DIR "u_1d_" + std::to_string(nx) + "x" + std::to_string(ny) + ".grid";
  Zap::IO::IncCellWriter<ActiveFloat, PassiveFloat> grid_writer{u_file, grid};

  const auto t_file = OUTPUT_DIR "t_1d_" + std::to_string(nx) + "x" + std::to_string(ny) + ".mat";
  Zap::IO::IncMatrixWriter<PassiveFloat, 1, 1, 0> t_writer(t_file, 1, 1, 0);

  const auto res = solve_2d_burgers<ExtendType::MAX>(grid, tend, grid_writer, t_writer, 0.25);
  if (!res.has_value()) {
    Igor::Warn("Solver for {}x{}-grid failed.", nx, ny);
    return false;
  }

  out << nx << 'x' << ny << ":\n";
  print_shock_error(*res, tend, out);
  print_solution_error(*res, tend, 10'000, out);
  out << "--------------------------------------------------------------------------------\n";

  if (dump_npy) {
    if (!dump_to_npy(*res)) { return false; }
  }

  Igor::Info("Solver for {}x{}-grid finished successfully.", nx, ny);
  return true;
}

// -------------------------------------------------------------------------------------------------
auto main(int argc, char** argv) -> int {
  const auto args = parse_args(argc, argv);
  if (!args.has_value()) { return 1; }

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

  Igor::Info("Save results in `" OUTPUT_DIR "`.");
  Igor::Info("run-benchmark    = {}", args->run_benchmark);
  Igor::Info("square-benchmark = {}", args->square_benchmark);
  Igor::Info("tend             = {}", args->tend);

  if (args->run_benchmark) {
    constexpr auto output_file = "./x_ramp_results.txt";
    std::ofstream out(output_file);
    if (!out) {
      Igor::Warn("Could not open output file `{}`: {}", output_file, std::strerror(errno));
    }

    bool all_success = true;
    // constexpr std::array ns = {
    //     10UZ,  15UZ,  20UZ,  25UZ,  30UZ,  35UZ,  40UZ,  45UZ,  50UZ,  55UZ,  60UZ,  65UZ,
    //     70UZ,  75UZ,  80UZ,  85UZ,  90UZ,  95UZ,  100UZ, 105UZ, 110UZ, 115UZ, 120UZ, 125UZ,
    //     130UZ, 135UZ, 140UZ, 145UZ, 150UZ, 155UZ, 160UZ, 165UZ, 170UZ, 175UZ, 180UZ, 185UZ,
    //     190UZ, 195UZ, 200UZ, 205UZ, 210UZ, 215UZ, 220UZ, 225UZ, 230UZ, 235UZ, 240UZ, 245UZ,
    //     250UZ, 255UZ, 260UZ, 265UZ, 270UZ, 275UZ, 280UZ, 285UZ, 290UZ, 295UZ, 300UZ,
    // };
    // for (size_t n : ns) {
    for (size_t n = 10; n <= 400UZ; n += 10) {
      const auto success = run(n, args->square_benchmark ? n : 3, args->tend, out, false);
      all_success        = all_success && success;
    }

    Igor::Info("Saved errors to `{}`.", output_file);
    return all_success ? 0 : 1;
  } else {
    Igor::Info("nx               = {}", args->nx);
    Igor::Info("ny               = {}", args->ny);

    const auto success = run(args->nx, args->ny, args->tend, std::cout, args->dump_npy);
    return success ? 0 : 1;
  }
}
