#include <filesystem>

#include "CellBased/Cell.hpp"
#include "CellBased/Solver.hpp"
#include "IO/IncCellWriter.hpp"
#include "IO/IncMatrixWriter.hpp"
#include "IO/VTKWriter.hpp"

#include "Igor.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "cell_based/"

auto main(int argc, char** argv) -> int {
  using Float          = double;
  constexpr size_t DIM = 1;

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

  auto parse_size_t = [](const char* cstr) -> size_t {
    char* end        = nullptr;
    const size_t val = std::strtoul(cstr, &end, 10);
    if (end != cstr + std::strlen(cstr)) {  // NOLINT
      Igor::Panic("String `{}` contains non-digits.", cstr);
    }
    if (val == 0UL) {
      Igor::Panic("Could not parse string `{}` to size_t.", cstr);
    }
    if (val == std::numeric_limits<unsigned long>::max()) {
      Igor::Panic("Could not parse string `{}` to size_t: {}", cstr, std::strerror(errno));
    }
    return val;
  };
  auto parse_double = [](const char* cstr) -> double {
    char* end        = nullptr;
    const double val = std::strtod(cstr, &end);
    if (val == HUGE_VAL) {
      Igor::Panic("Could not parse string `{}` to double: Out of range.", cstr);
    }
    if (cstr == end) {
      Igor::Panic("Could not parse string `{}` to double.", cstr);
    }
    return val;
  };

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

  const Float x_min = 0.0;
  const Float x_max = 5.0;
  const Float y_min = 0.0;
  const Float y_max = 5.0;

  auto grid = Zap::CellBased::Grid<Float, DIM>::Uniform(x_min, x_max, nx, y_min, y_max, ny);
  grid.periodic_boundary(Zap::CellBased::LEFT | Zap::CellBased::RIGHT);
  grid.same_value_boundary(Zap::CellBased::BOTTOM | Zap::CellBased::TOP);

  // grid.dump_cells(std::cout);
  // return 0;

  auto u0 = [=](Float x, [[maybe_unused]] Float y) -> Float {
    // return (x + y) * static_cast<Float>((std::pow(x - x_min, 2) + std::pow(y - y_min, 2)) <=
    //                                     std::pow((x_min + x_max + y_min + y_max) / 4, 2));

    return (x - x_min) * static_cast<Float>((x - x_min) < 0.5 * (x_max - x_min));

    // if ((x - x_min) < 0.25 * (x_max - x_min)) {
    //   return x;
    // } else if ((x - x_min) < 0.5 * (x_max - x_min)) {
    //   return 0.25 * (x_max - x_min) - (x - 0.25 * (x_max - x_min));
    // } else {
    //   return 0;
    // }
  };
  // grid.fill_center(u0);
  grid.fill_four_point(u0);

  constexpr auto flux = [](const auto& u) constexpr noexcept {
    return static_cast<Float>(0.5) * u * u;
  };

  // From LeVeque: Numerical Methods for Conservation Laws 2nd edition (13.24)
  constexpr auto godunov_flux = [flux]<typename T>(const T& u_left,
                                                   const T& u_right) constexpr noexcept -> T {
    constexpr auto zero = static_cast<T>(0);
    if (u_left <= u_right) {
      // min u in [u_left, u_right] f(u) = 0.5 * u^2
      if (u_left <= zero && u_right >= zero) {
        return flux(zero);
      }
      return flux(std::min(std::abs(u_left), std::abs(u_right)));
    }
    // max u in [u_right, u_left] f(u) = 0.5 * u^2
    return flux(std::max(std::abs(u_left), std::abs(u_right)));
  };

  // Zap::IO::NoopWriter grid_writer{};
  // Zap::IO::VTKWriter<Zap::IO::VTKFormat::UNSTRUCTURED_GRID> grid_writer{OUTPUT_DIR "u_grid"};

  constexpr auto u_file = OUTPUT_DIR "u.grid";
  Zap::IO::IncCellWriter<Float, DIM> grid_writer{u_file, grid};

  // Zap::IO::NoopWriter t_writer{};

  constexpr auto t_file = OUTPUT_DIR "t.mat";
  Zap::IO::IncMatrixWriter<Float, 1, 1, 0> t_writer(t_file, 1, 1, 0);

  IGOR_TIME_SCOPE("Solver") {
    Zap::CellBased::Solver solver(godunov_flux, godunov_flux);
    if (!solver.solve(grid, static_cast<Float>(tend), grid_writer, t_writer).has_value()) {
      Igor::Warn("Solver failed.");
      return 1;
    }
  }
  Igor::Info("Solver finished successfully.");
  Igor::Info("Saved grid to {}.", u_file);
  Igor::Info("Saved time steps to {}.", t_file);
}
