#include <filesystem>

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
                       Float tend) noexcept {
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

  Igor::Info("mean_shock_x                                     = {:.8f}", mean_shock_x);
  Igor::Info("std_dev_shock_x                                  = {:.8f}", std_dev_shock_x);

  const auto expected_shock_x = expected_shock_location(tend);
  const auto abs_err          = std::abs(mean_shock_x - expected_shock_x);
  const auto rel_err          = abs_err / expected_shock_x;

  Igor::Info("expected_shock_x                                 = {:.8f}", expected_shock_x);
  Igor::Info("|mean_shock_x - expected_shock_x|                = {:.8f}", abs_err);
  Igor::Info("|mean_shock_x - expected_shock_x| / mean_shock_x = {:.8f} (≈{:.2f}%)",
             rel_err,
             rel_err * 100);
}

// -------------------------------------------------------------------------------------------------
void print_solution_error(const Zap::CellBased::UniformGrid<Float, DIM>& numerical_solution,
                          Float tend,
                          size_t n) noexcept {
  assert(n > 0);

  Float L1_error = 0;
  const Float dx = (X_MAX - X_MIN) / static_cast<Float>(n - 1);
  for (size_t i = 0; i < n; ++i) {
    const Float x = static_cast<Float>(i) * dx + X_MIN;
    const Float f = std::abs(numerical_solution.eval(Zap::CellBased::Point<Float>{x, 1.0})(0) -
                             analytical_quasi_1d(x, tend));
    L1_error += f * dx;
  }
  Igor::Info("L1 error = {:.8f}", L1_error);
}

// -------------------------------------------------------------------------------------------------
auto main(int argc, char** argv) -> int {
  if (argc < 4) {
    Igor::Warn("Usage: {} <nx> <ny> <tend>", *argv);
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

  const auto nx   = parse_size_t(argv[1]);  // NOLINT
  const auto ny   = parse_size_t(argv[2]);  // NOLINT
  const auto tend = parse_double(argv[3]);  // NOLINT
  assert(nx < std::numeric_limits<int>::max());
  assert(ny < std::numeric_limits<int>::max());

  if (tend <= 0.0) {
    Igor::Warn("tend must be larger than 0, but is {}.", tend);
    return 1;
  }

  Igor::Info("nx = {}", nx);
  Igor::Info("ny = {}", ny);
  Igor::Info("tend = {}", tend);

  Zap::CellBased::UniformGrid<Float, DIM> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.same_value_boundary();

  auto u0 = [=](Float x, Float /*y*/) -> Float {
    static_assert(DIM == 1);
    return (x - X_MIN) * static_cast<Float>((x - X_MIN) < (X_MAX - X_MIN) / 2);
  };

  auto init_shock = [=]<typename T>(T t) -> Zap::CellBased::Point<T> {
    assert(t >= 0 && t <= 1);
    return Zap::CellBased::Point<T>{
        (X_MAX - X_MIN) / 2,
        t * (Y_MAX - Y_MIN) + Y_MIN,
    };
  };

  IGOR_TIME_SCOPE("Cutting the grid") {
    if (!grid.cut_curve(init_shock)) { return 1; }
  }

  grid.fill_four_point(u0);

  constexpr auto u_file = OUTPUT_DIR "u_1d.grid";
  Zap::IO::IncCellWriter<Float, DIM> grid_writer{u_file, grid};

  constexpr auto t_file = OUTPUT_DIR "t_1d.mat";
  Zap::IO::IncMatrixWriter<Float, 1, 1, 0> t_writer(t_file, 1, 1, 0);

  IGOR_TIME_SCOPE("Solver") {
    Zap::CellBased::Solver solver(Zap::CellBased::SingleEq::A{}, Zap::CellBased::SingleEq::B{});
    const auto res = solver.solve(grid, static_cast<Float>(tend), grid_writer, t_writer, 0.5);
    if (!res.has_value()) {
      Igor::Warn("Solver failed.");
      return 1;
    }

    std::cout
        << "--------------------------------------------------------------------------------\n";
    print_shock_error(*res, tend);
    std::cout
        << "--------------------------------------------------------------------------------\n";
    print_solution_error(*res, tend, 1000);
    std::cout
        << "--------------------------------------------------------------------------------\n";
  }
  Igor::Info("Solver finished successfully.");
  Igor::Info("Saved grid to {}.", u_file);
  Igor::Info("Saved time steps to {}.", t_file);
}
