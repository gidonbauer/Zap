#include <filesystem>
#include <iomanip>

// #define USE_FLUX_FOR_CARTESIAN

#include "CellBased/EigenDecomp.hpp"
#include "CellBased/Solver.hpp"
#include "IO/IncCellWriter.hpp"
#include "IO/IncMatrixWriter.hpp"

#include "Igor.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "example_x_ramp/"

// - Setup -----------------------------------------------------------------------------------------
using Float           = double;
constexpr size_t DIM  = 1;
constexpr Float X_MIN = 0.0;
constexpr Float X_MAX = 2.0;
constexpr Float Y_MIN = 0.0;
constexpr Float Y_MAX = 2.0;

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto parse_size_t(const char* cstr) -> size_t {
  char* end        = nullptr;
  const size_t val = std::strtoul(cstr, &end, 10);
  if (end != cstr + std::strlen(cstr)) {  // NOLINT
    Igor::Panic("String `{}` contains non-digits.", cstr);
  }
  if (val == 0UL) { Igor::Panic("Could not parse string `{}` to size_t.", cstr); }
  if (val == std::numeric_limits<unsigned long>::max()) {
    Igor::Panic("Could not parse string `{}` to size_t: {}", cstr, std::strerror(errno));
  }
  return val;
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto parse_double(const char* cstr) -> double {
  char* end        = nullptr;
  const double val = std::strtod(cstr, &end);
  if (val == HUGE_VAL) {
    Igor::Panic("Could not parse string `{}` to double: Out of range.", cstr);
  }
  if (cstr == end) { Igor::Panic("Could not parse string `{}` to double.", cstr); }
  return val;
}

[[nodiscard]] auto init_shock(Float t) -> Zap::CellBased::Point<Float> {
  assert(t >= 0 && t <= 1);
  return Zap::CellBased::Point<Float>{
      (X_MAX - X_MIN) / 2,
      t * (Y_MAX - Y_MIN) + Y_MIN,
  };
};

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto chi(Float x, Float x_min, Float x_max) noexcept -> Float {
  return static_cast<Float>(x >= x_min && x <= x_max);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto analytical_quasi_1d(Float x, Float t) noexcept -> Float {
  return x / (1 + t) * chi(x, Float{0}, std::sqrt(1 + t));
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto expected_shock_location(Float t) noexcept -> Float {
  return std::sqrt(1 + t);
}

// -------------------------------------------------------------------------------------------------
void print_shock_error(const Zap::CellBased::UniformGrid<Float, DIM>& numerical_solution,
                       Float tend,
                       std::ostream& out) noexcept {
  const auto final_shock = numerical_solution.get_shock_curve();

  const auto mean_shock_x = std::transform_reduce(std::cbegin(final_shock),
                                                  std::cend(final_shock),
                                                  Float{0},
                                                  std::plus<>{},
                                                  [](const Zap::CellBased::Point<Float>& p) {
                                                    return p(Zap::CellBased::X);
                                                  }) /
                            static_cast<Float>(final_shock.size());

  const auto std_dev_shock_x =
      std::sqrt(std::transform_reduce(std::cbegin(final_shock),
                                      std::cend(final_shock),
                                      Float{0},
                                      std::plus<>{},
                                      [=](const Zap::CellBased::Point<Float>& p) {
                                        return std::pow(p(Zap::CellBased::X) - mean_shock_x, 2);
                                      }) /
                static_cast<Float>(final_shock.size()));

  out << "mean_shock_x                                     = " << std::setprecision(8)
      << mean_shock_x << '\n';
  out << "std_dev_shock_x                                  = " << std::setprecision(8)
      << std_dev_shock_x << '\n';

  const auto expected_shock_x = expected_shock_location(tend);
  const auto abs_err          = std::abs(mean_shock_x - expected_shock_x);
  const auto rel_err          = abs_err / expected_shock_x;

  out << "expected_shock_x                                 = " << std::setprecision(8)
      << expected_shock_x << '\n';
  out << "|mean_shock_x - expected_shock_x|                = " << std::setprecision(8) << abs_err
      << '\n';
  out << "|mean_shock_x - expected_shock_x| / mean_shock_x = " << std::setprecision(8) << rel_err
      << " (≈" << std::setprecision(2) << rel_err * 100 << "%)\n";
}

// -------------------------------------------------------------------------------------------------
void print_solution_error(const Zap::CellBased::UniformGrid<Float, DIM>& numerical_solution,
                          Float tend,
                          size_t n,
                          std::ostream& out) noexcept {
  assert(n > 0);
  if (n % 2 == 1) { n += 1; }

  // Use simpsons rule
  Float L1_error = 0;
  const Float dx = (X_MAX - X_MIN) / static_cast<Float>(n);
  for (size_t i = 1; i <= n / 2; ++i) {
    {
      const Float x = static_cast<Float>(2 * i - 2) * dx + X_MIN;
      L1_error += std::abs(numerical_solution.eval(Zap::CellBased::Point<Float>{x, 1.0})(0) -
                           analytical_quasi_1d(x, tend));
    }
    {
      const Float x = static_cast<Float>(2 * i - 1) * dx + X_MIN;
      L1_error += 4 * std::abs(numerical_solution.eval(Zap::CellBased::Point<Float>{x, 1.0})(0) -
                               analytical_quasi_1d(x, tend));
    }
    {
      const Float x = static_cast<Float>(2 * i) * dx + X_MIN;
      L1_error += std::abs(numerical_solution.eval(Zap::CellBased::Point<Float>{x, 1.0})(0) -
                           analytical_quasi_1d(x, tend));
    }
  }
  L1_error *= dx / 3;
  out << "L1 error                                         = " << std::setprecision(8) << L1_error
      << '\n';
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto run(size_t nx, size_t ny, Float tend, std::ostream& out) noexcept -> bool {
  Zap::CellBased::UniformGrid<Float, DIM> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.same_value_boundary();

  if (!grid.cut_curve(init_shock)) { return false; }

  constexpr auto u0 = [](Float x, Float /*y*/) { return analytical_quasi_1d(x, Float{0}); };
  grid.fill_four_point(u0);

  const auto u_file = OUTPUT_DIR "u_1d_" + std::to_string(nx) + "_" + std::to_string(ny) + ".grid";
  Zap::IO::IncCellWriter<Float, DIM> grid_writer{u_file, grid};

  const auto t_file = OUTPUT_DIR "t_1d_" + std::to_string(nx) + "_" + std::to_string(ny) + ".mat";
  Zap::IO::IncMatrixWriter<Float, 1, 1, 0> t_writer(t_file, 1, 1, 0);

  Zap::CellBased::Solver solver(Zap::CellBased::SingleEq::A{}, Zap::CellBased::SingleEq::B{});
  const auto res = solver.solve(grid, static_cast<Float>(tend), grid_writer, t_writer, 0.5);
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
  if (argc < 2) {
    Igor::Warn("Usage: {} <tend>", *argv);
    return 1;
  }

  {
    std::error_code ec;
    std::filesystem::create_directories(OUTPUT_DIR, ec);
    if (ec) {
      Igor::Warn("Could not create directory `{}`: {}", OUTPUT_DIR, ec.message());
      return 1;
    }
  }

  const auto tend = static_cast<Float>(parse_double(argv[1]));  // NOLINT

  if (tend <= 0.0) {
    Igor::Warn("tend must be larger than 0, but is {}.", tend);
    return 1;
  }

  Igor::Info("tend = {}", tend);
  Igor::Info("Save results in `" OUTPUT_DIR "`.");

  constexpr auto output_file = "./x_ramp_results.txt";
  std::ofstream out(output_file);
  if (!out) {
    Igor::Warn("Could not open output file `{}`: {}", output_file, std::strerror(errno));
  }

  bool all_success        = true;
  constexpr std::array ns = {
      3UZ,  5UZ,  7UZ,   9UZ,   11UZ,  15UZ,  21UZ,  31UZ,  41UZ,  51UZ,  61UZ,  71UZ,
      81UZ, 91UZ, 101UZ, 111UZ, 121UZ, 131UZ, 141UZ, 151UZ, 161UZ, 171UZ, 181UZ, 191UZ,
  };
  for (size_t n : ns) {
    const auto success = run(n, n, tend, out);
    all_success        = all_success && success;
  }

  Igor::Info("Saved errors to `{}`.", output_file);
  return all_success ? 0 : 1;
}
