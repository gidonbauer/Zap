#include <filesystem>
#include <optional>

#include "CellBased/Solver.hpp"
#include "IO/IncCellWriter.hpp"
#include "IO/IncMatrixWriter.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Macros.hpp"
#include "Igor/Timer.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "example_inverse_ramp/"

// - Setup -----------------------------------------------------------------------------------------
using Float           = double;
constexpr Float X_MIN = 0.0;
constexpr Float X_MAX = 1.0;
constexpr Float Y_MIN = 0.0;
constexpr Float Y_MAX = 1.0;

// - U0 --------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto u0(Float x, Float /*y*/) noexcept -> Float {
  return (X_MAX - x) * static_cast<Float>(x >= (X_MAX + X_MIN) / 2);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto run(size_t nx, size_t ny, Float tend, Float CFL_safety_factor) noexcept -> bool {
  Zap::CellBased::UniformGrid<Float, Float> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.same_value_boundary();

  {
    std::vector<Zap::CellBased::SimCoord<Float>> points{
        {(X_MAX - X_MIN) / 2, 0.0},
        {(X_MAX - X_MIN) / 2, Y_MAX},
    };
    if (!grid.cut_piecewise_linear<Zap::CellBased::ExtendType::MAX>(points)) {
      Igor::Warn("Could not cut the grid with points {}", points);
      return false;
    }
  }

  grid.fill_four_point(u0);

  const auto u_file = OUTPUT_DIR "u_" + std::to_string(nx) + "x" + std::to_string(ny) + ".grid";
  Zap::IO::IncCellWriter<Float, Float> grid_writer{u_file, grid};

  const auto t_file = OUTPUT_DIR "t_" + std::to_string(nx) + "x" + std::to_string(ny) + ".mat";
  Zap::IO::IncMatrixWriter<Float, 1, 1, 0> t_writer(t_file, 1, 1, 0);

  Zap::CellBased::Solver<Zap::CellBased::ExtendType::MAX> solver;
  const auto res = solver.solve(grid, tend, grid_writer, t_writer, CFL_safety_factor);
  if (!res.has_value()) {
    Igor::Warn("Solver for {}x{}-grid failed.", nx, ny);
    return false;
  }

  Igor::Info("Solver for {}x{}-grid finished successfully.", nx, ny);
  return true;
}

// -------------------------------------------------------------------------------------------------
auto main() -> int {
  {
    std::error_code ec;
    std::filesystem::create_directories(OUTPUT_DIR, ec);
    if (ec) {
      Igor::Warn("Could not create directory `{}`: {}", OUTPUT_DIR, ec.message());
      return 1;
    }
  }

  Igor::Info("Save results in `" OUTPUT_DIR "`.");

  return run(100, 100, 1.0, 0.1) ? 0 : 1;
}
