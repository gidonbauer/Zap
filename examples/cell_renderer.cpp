#include "CellBased/Cell.hpp"
#include "IO/IncCellReader.hpp"
#include "IO/IncMatrixReader.hpp"
#include "Renderer/Canvas.hpp"
#include "Renderer/FFmpeg.hpp"

#define ZAP_PARALLEL_RENDERING

namespace Rd = Zap::Renderer;

template <Igor::detail::Level level>
constexpr void usage(std::string_view prog) {
  std::cerr << std::format("{}Usage: {} [--scale s] <u input file> <t input file> <output file>\n",
                           Igor::detail::level_repr(level),
                           prog);
};

auto main(int argc, char** argv) -> int {
  std::string u_input_file{};
  std::string t_input_file{};
  std::string output_file{};
  size_t scale = 1;

  argc -= 1;
  std::string_view prog = *argv;
  argv += 1;  // NOLINT
  auto parse_size_t = [](const char* cstr) -> size_t {
    char* end        = nullptr;
    const size_t val = std::strtoul(cstr, &end, 10);
    if (end != cstr + std::strlen(cstr)) {  // NOLINT
      Igor::Panic("String `{}` contains non-digits.", cstr);
    }
    if (val == 0UL) {
      Igor::Panic("Could not parse string `{}` to size_t or number is zero.", cstr);
    }
    if (val == std::numeric_limits<unsigned long>::max()) {
      Igor::Panic("Could not parse string `{}` to size_t: {}", cstr, std::strerror(errno));
    }
    return val;
  };

  while (argc > 0) {
    using namespace std::string_view_literals;
    if (*argv == "--help"sv || *argv == "-h"sv) {
      usage<Igor::detail::Level::INFO>(prog);
      return 0;
    } else if (*argv == "--scale"sv) {
      argc -= 1;
      argv += 1;  // NOLINT
      if (argc <= 0) {
        usage<Igor::detail::Level::WARN>(prog);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
                  << "Did not provide value for scale.\n";
        return 1;
      }
      scale = parse_size_t(*argv);
      if (scale == 0) {
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
                  << "Scale must be larger than zero but is " << scale << '\n';
      }
    } else if (u_input_file.empty()) {
      u_input_file = *argv;
    } else if (t_input_file.empty()) {
      t_input_file = *argv;
    } else if (output_file.empty()) {
      output_file = *argv;
    } else {
      usage<Igor::detail::Level::WARN>(prog);
      std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
                << "Provided too many arguments.\n";
      return 1;
    }
    argc -= 1;
    argv += 1;  // NOLINT
  }

  if (u_input_file.empty()) {
    usage<Igor::detail::Level::WARN>(prog);
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
              << "Did not provide u input file.\n";
    return 1;
  }
  if (t_input_file.empty()) {
    usage<Igor::detail::Level::WARN>(prog);
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
              << "Did not provide t input file.\n";
    return 1;
  }
  if (output_file.empty()) {
    usage<Igor::detail::Level::WARN>(prog);
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
              << "Did not provide output file.\n";
    return 1;
  }

  Igor::Info("scale        = {}", scale);
  Igor::Info("u input file = {}", u_input_file);
  Igor::Info("t input file = {}", t_input_file);
  Igor::Info("output file  = {}", output_file);

  // -----------------------------------------------------------------------------------------------

  using Float          = double;
  constexpr size_t DIM = 1;  // TODO: Implement renderer for 2D data

  try {
    Zap::IO::IncMatrixReader<Float> t_reader(t_input_file);
    Zap::IO::IncCellReader<Float, DIM> u_reader(u_input_file);

    constexpr size_t TEXT_HEIGHT      = 100UZ;
    constexpr size_t MIN_CANVAS_WIDTH = 600UZ;
    const size_t graph_width          = u_reader.nx() * scale;
    const size_t graph_height         = u_reader.ny() * scale;
    Rd::Canvas canvas(std::max(MIN_CANVAS_WIDTH, graph_width), graph_height + TEXT_HEIGHT);
    const Rd::Box text_box{
        .col    = 0,
        .row    = 0,
        .width  = canvas.width(),
        .height = TEXT_HEIGHT,
    };

    Igor::Info("Resolution: {}x{}", canvas.width(), canvas.height());

    const size_t avail_padding = canvas.width() - graph_width;
    const Rd::Box graph_box{
        .col    = avail_padding / 2UZ,
        .row    = text_box.height,
        .width  = graph_width,
        .height = graph_height,
    };
    constexpr std::string_view font_file = "../assets/LinLibertine_RI.ttf";
    if (!canvas.load_font(font_file)) {
      Igor::Warn("Could not load font from file `{}`.", font_file);
      return 1;
    }

    Igor::ScopeTimer timer{"Rendering"};

    const Rd::FFmpeg ffmpeg(canvas.width(), canvas.height(), output_file);
    while (true) {
      const auto got_next_t = t_reader.read_next<false>();
      const auto got_next_u = u_reader.read_next<false>();
      if (!got_next_u && !got_next_t) {
        break;
      } else if (!got_next_u) {
        Igor::Warn("Got more t-entries than u-entries.");
        return 1;
      } else if (!got_next_t) {
        Igor::Warn("Got more u-entries than t-entries.");
        return 1;
      }

      canvas.clear(Rd::RGB{.r = 0x00, .g = 0x00, .b = 0x00});

      if (const std::string str = std::format("t = {:.6f}", t_reader(0, 0));
          !canvas.draw_text(str, text_box, true)) {
        Igor::Warn("Could not render string `{}` to canvas.", str);
        return 1;
      }

      const auto& cells = u_reader.cells();
      assert(cells.size() > 0);
      auto min = cells[0].value;
      auto max = cells[0].value;
      for (const auto& cell : cells) {
        for (Eigen::Index i = 0; i < cell.value.rows(); ++i) {
          min(i) = std::min(min(i), cell.value(i));
          max(i) = std::max(max(i), cell.value(i));
        }
      }

      bool failed_drawing_rects = false;
#ifdef ZAP_PARALLEL_RENDERING
#pragma omp parallel for
#endif  // ZAP_PARALLEL_RENDERING
      for (const auto& cell : cells) {
        const Rd::Box rect{
            .col = static_cast<size_t>(
                std::round((cell.x_min - u_reader.x_min()) / (u_reader.x_max() - u_reader.x_min()) *
                           static_cast<Float>(u_reader.nx() * scale))),
            .row = static_cast<size_t>(
                std::round((cell.y_min - u_reader.y_min()) / (u_reader.y_max() - u_reader.y_min()) *
                           static_cast<Float>(u_reader.ny() * scale))),
            .width =
                static_cast<size_t>(std::round(cell.dx / (u_reader.x_max() - u_reader.x_min()) *
                                               static_cast<Float>(u_reader.nx() * scale))),
            .height =
                static_cast<size_t>(std::round(cell.dy / (u_reader.y_max() - u_reader.y_min()) *
                                               static_cast<Float>(u_reader.ny() * scale))),
        };
        const auto c = Rd::float_to_rgb(cell.value[0], min[0], max[0]);

        if (!canvas.draw_rect(rect, graph_box, c, true)) {
          Igor::Warn("Could not draw cell");
          failed_drawing_rects = true;
        }
      }
      if (failed_drawing_rects) {
        return 1;
      }

      if (!canvas.to_raw_stream(ffmpeg.stream())) {
        Igor::Warn("Could not write canvas to FFmpeg.");
        return 1;
      }
    }

  } catch (const std::exception& e) {
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC) << e.what() << '\n';
    return 1;
  }
}
