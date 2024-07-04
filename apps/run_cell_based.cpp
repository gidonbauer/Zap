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

  if (argc < 3) {
    Igor::Warn("Usage: {} <nx> <ny>", *argv);
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

  const auto nx = parse_size_t(argv[1]);  // NOLINT
  const auto ny = parse_size_t(argv[2]);  // NOLINT
  assert(nx < std::numeric_limits<int>::max());
  assert(ny < std::numeric_limits<int>::max());

  Igor::Info("nx = {}", nx);
  Igor::Info("ny = {}", ny);

  const Float x_min = 0.0;
  const Float x_max = 5.0;
  const Float y_min = 0.0;
  const Float y_max = 5.0;

  auto grid = Zap::CellBased::Grid<Float, DIM>::Uniform(x_min, x_max, nx, y_min, y_max, ny);
  grid.make_periodic();

  auto u0 = [=](Float x, Float y) {
    return (x + y) * static_cast<Float>((std::pow(x - x_min, 2) + std::pow(y - y_min, 2)) <=
                                        std::pow((x_min + x_max + y_min + y_max) / 4, 2));
  };
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
  // Zap::IO::VTKWriter<Zap::IO::VTKFormat::UNSTRUCTURED_GRID> grid_writer{OUTPUT_DIR "grid"};
  Zap::IO::IncCellWriter<Float, DIM> grid_writer{OUTPUT_DIR "u.grid", grid};

  // Zap::IO::NoopWriter t_writer{};
  Zap::IO::IncMatrixWriter<Float, 1, 1, 0> t_writer(OUTPUT_DIR "t.mat", 1, 1, 0);

  IGOR_TIME_SCOPE("Solver") {
    Zap::CellBased::Solver solver(godunov_flux, godunov_flux);
    if (!solver.solve(grid, Float{1}, grid_writer, t_writer).has_value()) {
      Igor::Warn("Solver failed.");
      return 1;
    }
  }
  Igor::Info("Solver finished successfully.");
}
