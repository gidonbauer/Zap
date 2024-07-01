#include <algorithm>
#include <chrono>

#include "IO/IncMatrixReader.hpp"
#include "Igor.hpp"
#include "Renderer/Canvas.hpp"
#include "Renderer/FFmpeg.hpp"

namespace Rd = Zap::Renderer;

// -------------------------------------------------------------------------------------------------
auto main(int argc, char** argv) -> int {
  using Float = double;

  if (argc < 4) {
    Igor::Warn("Usage: {} <u input file> <t input file> <output_file>", *argv);
    return 1;
  }
  const auto u_input_file = argv[1];  // NOLINT
  const auto t_input_file = argv[2];  // NOLINT
  const auto output_file  = argv[3];  // NOLINT

  Zap::IO::IncMatrixReader<Float> u_reader{u_input_file};
  Zap::IO::IncMatrixReader<Float> t_reader{t_input_file};
  if (u_reader.num_elem() != t_reader.num_elem()) {
    Igor::Warn("Incompatible number of elements in u ({}) and t ({}).", u_input_file, t_input_file);
    return 1;
  }

  constexpr size_t TEXT_HEIGHT      = 100UZ;
  constexpr size_t MIN_CANVAS_WIDTH = 600UZ;
  Rd::Canvas canvas(std::max(MIN_CANVAS_WIDTH, static_cast<size_t>(u_reader.cols())),
                    static_cast<size_t>(u_reader.rows()) + TEXT_HEIGHT);
  Rd::Box text_box{
      .col    = 0,
      .row    = 0,
      .width  = canvas.width(),
      .height = TEXT_HEIGHT,
  };

  const size_t avail_padding = canvas.width() - static_cast<size_t>(u_reader.cols());
  Rd::Box graph_box{
      .col    = avail_padding / 2UZ,
      .row    = text_box.height,
      .width  = static_cast<size_t>(u_reader.cols()),
      .height = static_cast<size_t>(u_reader.rows()),
  };
  constexpr std::string_view font_file = "../assets/LinLibertine_R.ttf";
  if (!canvas.load_font(font_file)) {
    Igor::Warn("Could not load font from file `{}`.", font_file);
    return 1;
  }

  const Rd::FFmpeg ffmpeg(canvas.width(), canvas.height(), output_file);

  Igor::ScopeTimer timer{"Rendering"};
  while (u_reader.read_next<false>() && t_reader.read_next<false>()) {
    canvas.clear({});

    if (const std::string str = std::format("t = {:.6f}", t_reader(0, 0));
        !canvas.draw_text(str, text_box, true)) {
      Igor::Warn("Could not render string `{}` to canvas.", str);
      return 1;
    }

    const auto& data      = u_reader.data();
    const auto [min, max] = std::minmax_element(std::cbegin(data), std::cend(data));
    if (!canvas.draw_buffer(
            data,
            graph_box,
            [&](Float v) { return Rd::float_to_rgb<Rd::Float2RGB::COLORMAP>(v, *min, *max); },
            u_reader.is_row_major(),
            true)) {
      Igor::Warn("Could not draw graph to canvas.");
      return 1;
    }

    if (!canvas.to_raw_stream(ffmpeg.stream())) {
      Igor::Warn("Could not write canvas to FFmpeg.");
      return 1;
    }
  }
}
