#include <filesystem>
#include <iomanip>
#include <optional>

// #define USE_FLUX_FOR_CARTESIAN

#include "CellBased/EigenDecomp.hpp"
#include "CellBased/Solver.hpp"
#include "IO/IncCellWriter.hpp"
#include "IO/IncMatrixWriter.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Macros.hpp"
#include "Igor/Timer.hpp"

#include "Common.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "example_x_ramp/"

// - Setup -----------------------------------------------------------------------------------------
using ActiveFloat            = double;
using PassiveFloat           = double;
constexpr size_t DIM         = 1;
constexpr PassiveFloat X_MIN = 0.0;
constexpr PassiveFloat X_MAX = 2.0;
constexpr PassiveFloat Y_MIN = 0.0;
constexpr PassiveFloat Y_MAX = 2.0;

// - Available command line options ----------------------------------------------------------------
struct Args {
  bool run_benchmark = false;
  size_t nx          = 25;
  size_t ny          = 25;
  PassiveFloat tend  = 1.0;
};

// - Print usage -----------------------------------------------------------------------------------
template <Igor::detail::Level level>
void usage(std::string_view prog, std::ostream& out) noexcept {
  Args args{};

  out << Igor::detail::level_repr(level) << "Usage: " << prog
      << " [--run-benchmark] [--nx nx] [--ny ny] [--tend tend]\n";
  out << "\t--run-benchmark   Run the solver for different grid sizes and compare result against "
         "high-resolution Godunov solver, default is "
      << std::boolalpha << args.run_benchmark << '\n';
  out << "\t--nx              Number of cells in x-direction, default is " << args.nx << '\n';
  out << "\t--ny              Number of cells in y-direction, default is " << args.ny << '\n';
  out << "\t--tend            Final time for simulation, default is " << args.tend << '\n';
}

// - Parse command line arguments ------------------------------------------------------------------
[[nodiscard]] auto parse_args(int& argc, char**& argv) noexcept -> std::optional<Args> {
  Args args{};

  const auto prog = next_arg(argc, argv);
  IGOR_ASSERT(prog.has_value(), "Could not get program name from command line arguments.");

  for (auto arg = next_arg(argc, argv); arg.has_value(); arg = next_arg(argc, argv)) {
    using namespace std::string_view_literals;
    if (arg == "--run-benchmark"sv) {
      args.run_benchmark = true;
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
[[nodiscard]] auto init_shock(PassiveFloat t) -> Zap::CellBased::SimCoord<ActiveFloat> {
  assert(t >= 0 && t <= 1);
  return {
      .x = (X_MAX - X_MIN) / 2,
      .y = t * (Y_MAX - Y_MIN) + Y_MIN,
  };
};

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto chi(PassiveFloat x, PassiveFloat x_min, PassiveFloat x_max) noexcept
    -> PassiveFloat {
  return static_cast<PassiveFloat>(x >= x_min && x <= x_max);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto expected_shock_location(PassiveFloat t, ActiveFloat eps) noexcept
    -> ActiveFloat {
  return std::sqrt(1 + (1 + eps) * t);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto
analytical_quasi_1d(PassiveFloat x, PassiveFloat t, ActiveFloat eps) noexcept -> ActiveFloat {
  return (1 + eps) * x / (1 + (1 + eps) * t) *
         chi(x, PassiveFloat{0}, expected_shock_location(t, eps));
}

// -------------------------------------------------------------------------------------------------
void print_shock_error(
    const Zap::CellBased::UniformGrid<ActiveFloat, PassiveFloat, DIM>& numerical_solution,
    PassiveFloat tend,
    std::ostream& out) noexcept {
  const auto final_shock = numerical_solution.get_shock_curve();

  const auto mean_shock_x = std::transform_reduce(std::cbegin(final_shock),
                                                  std::cend(final_shock),
                                                  ActiveFloat{0},
                                                  std::plus<>{},
                                                  [](const auto& p) { return p.x; }) /
                            static_cast<PassiveFloat>(final_shock.size());

  const auto std_dev_shock_x = std::sqrt(
      std::transform_reduce(std::cbegin(final_shock),
                            std::cend(final_shock),
                            ActiveFloat{0},
                            std::plus<>{},
                            [=](const auto& p) { return std::pow(p.x - mean_shock_x, 2); }) /
      static_cast<PassiveFloat>(final_shock.size()));

  out << "mean_shock_x     = " << std::setprecision(8) << mean_shock_x << '\n';
  out << "std_dev_shock_x  = " << std::setprecision(8) << std_dev_shock_x << '\n';

  const auto expected_shock_x = expected_shock_location(tend, ActiveFloat{0});
  const auto abs_err          = std::abs(mean_shock_x - expected_shock_x);
  const auto rel_err          = abs_err / expected_shock_x;

  out << "expected_shock_x = " << std::setprecision(8) << expected_shock_x << '\n';
  out << "abs_error_pos    = " << std::setprecision(8) << abs_err << '\n';
  out << "rel_error_pos    = " << std::setprecision(8) << rel_err << '\n';
}

// -------------------------------------------------------------------------------------------------
void print_solution_error(
    const Zap::CellBased::UniformGrid<ActiveFloat, PassiveFloat, DIM>& numerical_solution,
    PassiveFloat tend,
    size_t n,
    std::ostream& out) noexcept {
  assert(n > 0);
  if (n % 2 == 1) { n += 1; }

  // Use simpsons rule
  ActiveFloat L1_error  = 0;
  const PassiveFloat dx = (X_MAX - X_MIN) / static_cast<PassiveFloat>(n);
  for (size_t i = 1; i <= n / 2; ++i) {
    {
      const PassiveFloat x = static_cast<PassiveFloat>(2 * i - 2) * dx + X_MIN;
      L1_error += std::abs(
          numerical_solution.eval(Zap::CellBased::SimCoord<PassiveFloat>{.x = x, .y = 1.0})(0) -
          analytical_quasi_1d(x, tend, ActiveFloat{0}));
    }
    {
      const PassiveFloat x = static_cast<PassiveFloat>(2 * i - 1) * dx + X_MIN;
      L1_error += 4 * std::abs(numerical_solution.eval(
                                   Zap::CellBased::SimCoord<PassiveFloat>{.x = x, .y = 1.0})(0) -
                               analytical_quasi_1d(x, tend, ActiveFloat{0}));
    }
    {
      const PassiveFloat x = static_cast<PassiveFloat>(2 * i) * dx + X_MIN;
      L1_error += std::abs(
          numerical_solution.eval(Zap::CellBased::SimCoord<PassiveFloat>{.x = x, .y = 1.0})(0) -
          analytical_quasi_1d(x, tend, ActiveFloat{0}));
    }
  }
  L1_error *= dx / 3;
  out << "L1_error         = " << std::setprecision(8) << L1_error << '\n';
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto run(size_t nx, size_t ny, PassiveFloat tend, std::ostream& out) noexcept
    -> bool {
  Zap::CellBased::UniformGrid<ActiveFloat, PassiveFloat, DIM> grid(
      X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.same_value_boundary();

  // if (!grid.cut_curve(init_shock)) { return false; }
  {
    std::vector<Zap::CellBased::SimCoord<PassiveFloat>> points{
        {.x = (X_MAX - X_MIN) / 2, .y = 0.0},
        {.x = (X_MAX - X_MIN) / 2, .y = (Y_MAX - Y_MIN) + Y_MIN},
    };
    if (!grid.cut_piecewise_linear(points)) { return false; }
  }

  constexpr auto u0 = [](PassiveFloat x, PassiveFloat /*y*/) {
    return analytical_quasi_1d(x, PassiveFloat{0}, ActiveFloat{0});
  };
  grid.fill_four_point(u0);

  const auto u_file = OUTPUT_DIR "u_1d_" + std::to_string(nx) + "x" + std::to_string(ny) + ".grid";
  Zap::IO::IncCellWriter<ActiveFloat, PassiveFloat, DIM> grid_writer{u_file, grid};

  const auto t_file = OUTPUT_DIR "t_1d_" + std::to_string(nx) + "x" + std::to_string(ny) + ".mat";
  Zap::IO::IncMatrixWriter<PassiveFloat, 1, 1, 0> t_writer(t_file, 1, 1, 0);

  auto solver = Zap::CellBased::make_solver<Zap::CellBased::ExtendType::MAX>(
      Zap::CellBased::SingleEq::A{}, Zap::CellBased::SingleEq::B{});
  const auto res = solver.solve(grid, tend, grid_writer, t_writer, 0.25);
  if (!res.has_value()) {
    Igor::Warn("Solver for {}x{}-grid failed.", nx, ny);
    return false;
  }

  out << nx << 'x' << ny << ":\n";
  print_shock_error(*res, tend, out);
  print_solution_error(*res, tend, 10'000, out);
  out << "--------------------------------------------------------------------------------\n";

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
  Igor::Info("run-benchmark = {}", args->run_benchmark);
  Igor::Info("tend          = {}", args->tend);

  if (args->run_benchmark) {
    constexpr auto output_file = "./x_ramp_results.txt";
    std::ofstream out(output_file);
    if (!out) {
      Igor::Warn("Could not open output file `{}`: {}", output_file, std::strerror(errno));
    }

    bool all_success        = true;
    constexpr std::array ns = {
        3UZ,   5UZ,   10UZ,  15UZ,  20UZ,  25UZ,  30UZ,  35UZ,  40UZ,  45UZ,
        50UZ,  55UZ,  60UZ,  65UZ,  70UZ,  75UZ,  80UZ,  85UZ,  90UZ,  95UZ,
        100UZ, 105UZ, 110UZ, 115UZ, 120UZ, 125UZ, 130UZ, 135UZ, 140UZ, 145UZ,
        150UZ, 155UZ, 160UZ, 165UZ, 170UZ, 175UZ, 180UZ, 185UZ, 190UZ, 195UZ,
        200UZ, 205UZ, 210UZ, 215UZ, 220UZ, 225UZ, 230UZ, 235UZ, 240UZ, /*245UZ,*/ 250UZ,
        255UZ, 260UZ, 265UZ, 270UZ, 275UZ, 280UZ, 285UZ, 290UZ, 295UZ, 300UZ,
    };
    for (size_t n : ns) {
      const auto success = run(n, 3, args->tend, out);
      all_success        = all_success && success;
    }

    Igor::Info("Saved errors to `{}`.", output_file);
    return all_success ? 0 : 1;
  } else {
    Igor::Info("nx            = {}", args->nx);
    Igor::Info("ny            = {}", args->ny);

    const auto success = run(args->nx, args->ny, args->tend, std::cout);
    return success ? 0 : 1;
  }
}
