#if 1

#include "CellBased/Cell.hpp"
#include "IO/IncCellReader.hpp"
#include "IO/IncCellWriter.hpp"

auto main() -> int {
  using Float          = double;
  constexpr size_t DIM = 2;

  const Float x_min = -1.0;
  const Float x_max = 1.0;
  const size_t nx   = 10;
  const Float y_min = -2.5;
  const Float y_max = 2.5;
  const size_t ny   = 60;
  auto grid         = Zap::CellBased::Grid<Float, DIM>::Uniform(x_min, x_max, nx, y_min, y_max, ny);
  auto u0           = [=](Float x, Float y) {
    return Eigen::Vector<Float, 2>{
        (x + y) * static_cast<Float>((std::pow(x - x_min, 2) + std::pow(y - y_min, 2)) <=
                                     std::pow((x_min + x_max + y_min + y_max) / 4, 2)),
        (x - x_min) * static_cast<Float>(x <= (x_min + x_max) / 2),
    };
  };
  grid.fill_four_point(u0);
  for (size_t i = 0; i < 5; ++i) {
    const auto& cell = grid[i];
    Igor::Info("{} => {{ .x_min = {}, .dx = {}, .y_min = {}, .dy = {}, .value = {} }}",
               i,
               cell.x_min,
               cell.dx,
               cell.y_min,
               cell.dy,
               cell.value);
  }

  std::string filename = "test.grid";

  try {
    Zap::IO::IncCellWriter writer(filename, grid);

    if (!writer.write_data(grid)) {
      Igor::Warn("Could not write data to file.");
      return 1;
    }

    Igor::Info("Wrote grid to `{}`.", filename);
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }

  try {
    Zap::IO::IncCellReader<Float, DIM> reader(filename);

    if (!reader.read_next()) {
      Igor::Warn("Could not read next grid.");
      return 1;
    }

    for (size_t i = 0; i < 5; ++i) {
      const auto& cell = reader[i];
      Igor::Info("{} => {{ .x_min = {}, .dx = {}, .y_min = {}, .dy = {}, .value = {} }}",
                 i,
                 cell.x_min,
                 cell.dx,
                 cell.y_min,
                 cell.dy,
                 cell.value);
    }

    Igor::Info("Read grid from `{}`.", filename);
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}

#else

#include "Renderer/Canvas.hpp"

namespace Rd = Zap::Renderer;

auto main() -> int {
  const size_t width  = 800;
  const size_t height = 600;
  Rd::Canvas canvas{width, height};

  canvas.draw_rect(200, 100, 300, 400, Rd::RGB{.r = 0xFF, .g = 0x00, .b = 0x00});
  canvas.draw_rect(100, 300, 300, 200, Rd::RGB{.r = 0x00, .g = 0xFF, .b = 0x00});
  canvas.draw_rect(500, 0, 100, 100, Rd::RGB{.r = 0x00, .g = 0x00, .b = 0xFF});

  const std::string filename = "cell.ppm";
  if (!canvas.to_ppm(filename)) {
    Igor::Warn("Could not save canvas to file.");
    return 1;
  }
}
#endif
