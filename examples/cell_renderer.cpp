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
