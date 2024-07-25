#include <cmath>
#include <optional>

#include "CellBased/Cell.hpp"
#include "IO/IncCellReader.hpp"
#include "IO/IncMatrixReader.hpp"
#include "Renderer/Canvas.hpp"
#include "Renderer/FFmpeg.hpp"

#define ZAP_PARALLEL_RENDERING

namespace Rd = Zap::Renderer;

[[maybe_unused]] constexpr Rd::RGB BLACK  = {.r = 0x00, .g = 0x00, .b = 0x00};
[[maybe_unused]] constexpr Rd::RGB RED    = {.r = 0xFF, .g = 0x00, .b = 0x00};
[[maybe_unused]] constexpr Rd::RGB PINK   = {.r = 0xFF, .g = 0x00, .b = 0xFF};
[[maybe_unused]] constexpr Rd::RGB ORANGE = {.r = 0xFF, .g = 0xA5, .b = 0x00};

// -------------------------------------------------------------------------------------------------
struct Args {
  std::string u_input_file;
  std::string t_input_file;
  std::string output_file;
  size_t scale           = 1;
  bool keep_range        = false;
  bool same_range        = false;
  bool save_frame_images = false;
  size_t fps             = 60;
  bool two_dim           = false;
};

// -------------------------------------------------------------------------------------------------
template <Igor::detail::Level level>
constexpr void usage(std::string_view prog) {
  Args args{};

  std::cerr << Igor::detail::level_repr(level) << "Usage: " << prog
            << " [--scale s] [--keep-range] [--same-range] [--save-frame-img] [--fps f] [--2d] "
               "<u input file> <t input file> <output file>\n";
  std::cerr << "\t--scale           Number of pixels for single cell, default is " << args.scale
            << '\n';
  std::cerr << "\t--keep-range      Keep min and max value of first frame for colormap, default is "
            << std::boolalpha << args.keep_range << '\n';
  std::cerr << "\t--same-range      Same min and max value for both plots in 2d, default is "
            << std::boolalpha << args.same_range << '\n';
  std::cerr << "\t--save-frame-img  Save each individual frame as .ppm image, default is "
            << std::boolalpha << args.save_frame_images << '\n';
  std::cerr << "\t--fps             Frames per second for rendered video, default is " << args.fps
            << '\n';
  std::cerr
      << "\t--2d              Render solution of 2-dimensional solution (u in R^2), default is "
      << std::boolalpha << args.two_dim << '\n';
  std::cerr << "\tu input file      Solution field for each iteration\n";
  std::cerr << "\tt input file      Time for each iteration\n";
  std::cerr << "\toutput file       File for the rendered video\n";
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto parse_args(int argc, char** argv) noexcept -> std::optional<Args> {
  Args args{};

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
      return std::nullopt;
    } else if (*argv == "--scale"sv) {
      argc -= 1;
      argv += 1;  // NOLINT
      if (argc <= 0) {
        usage<Igor::detail::Level::WARN>(prog);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
                  << "Did not provide value for scale.\n";
        return std::nullopt;
      }
      args.scale = parse_size_t(*argv);
      if (args.scale == 0) {
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
                  << "Scale must be larger than zero but is " << args.scale << '\n';
      }
    } else if (*argv == "--fps"sv) {
      argc -= 1;
      argv += 1;  // NOLINT
      if (argc <= 0) {
        usage<Igor::detail::Level::WARN>(prog);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
                  << "Did not provide value for fps.\n";
        return std::nullopt;
      }
      args.fps = parse_size_t(*argv);
      if (args.fps == 0) {
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
                  << "FPS must be larger than zero but is " << args.fps << '\n';
      }
    } else if (*argv == "--keep-range"sv) {
      args.keep_range = true;
    } else if (*argv == "--same-range"sv) {
      args.same_range = true;
    } else if (*argv == "--save-frame-img"sv) {
      args.save_frame_images = true;
    } else if (*argv == "--2d"sv) {
      args.two_dim = true;
    } else if (args.u_input_file.empty()) {
      args.u_input_file = *argv;
    } else if (args.t_input_file.empty()) {
      args.t_input_file = *argv;
    } else if (args.output_file.empty()) {
      args.output_file = *argv;
    } else {
      usage<Igor::detail::Level::WARN>(prog);
      std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
                << "Provided too many arguments.\n";
      return std::nullopt;
    }
    argc -= 1;
    argv += 1;  // NOLINT
  }

  if (args.u_input_file.empty()) {
    usage<Igor::detail::Level::WARN>(prog);
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
              << "Did not provide u input file.\n";
    return std::nullopt;
  }
  if (args.t_input_file.empty()) {
    usage<Igor::detail::Level::WARN>(prog);
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
              << "Did not provide t input file.\n";
    return std::nullopt;
  }
  if (args.output_file.empty()) {
    usage<Igor::detail::Level::WARN>(prog);
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
              << "Did not provide output file.\n";
    return std::nullopt;
  }

  return args;
}

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
[[nodiscard]] auto min_max_cell_value(
    const std::vector<typename Zap::IO::IncCellReader<Float, DIM>::ReducedCell>& cells) noexcept
    -> std::pair<Eigen::Vector<Float, DIM>, Eigen::Vector<Float, DIM>> {
  Eigen::Vector<Float, DIM> min = [cell = cells.front()] {
    assert(cell.is_cartesian() || cell.is_cut());
    if (cell.is_cartesian()) {
      return cell.get_cartesian().value;
    } else {
      Eigen::Vector<Float, DIM> res;
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        res(i) = std::min(cell.get_cut().left_value(i), cell.get_cut().right_value(i));
      }
      return res;
    }
  }();

  Eigen::Vector<Float, DIM> max = [cell = cells.front()] {
    assert(cell.is_cartesian() || cell.is_cut());
    if (cell.is_cartesian()) {
      return cell.get_cartesian().value;
    } else {
      Eigen::Vector<Float, DIM> res;
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        res(i) = std::min(cell.get_cut().left_value(i), cell.get_cut().right_value(i));
      }
      return res;
    }
  }();

  for (const auto& cell : cells) {
    if (cell.is_cartesian()) {
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        min(i) = std::min(min(i), cell.get_cartesian().value(i));
        max(i) = std::max(max(i), cell.get_cartesian().value(i));
      }
    } else {
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        min(i) =
            std::min(min(i), std::min(cell.get_cut().left_value(i), cell.get_cut().right_value(i)));
        max(i) =
            std::max(max(i), std::max(cell.get_cut().left_value(i), cell.get_cut().right_value(i)));
      }
    }
  }

  return std::make_pair(min, max);
}

// -------------------------------------------------------------------------------------------------
template <typename Float>
[[nodiscard]] constexpr auto norm_point(Float x, Float min, Float max) noexcept -> Float {
  return (x - min) / (max - min);
}

// -------------------------------------------------------------------------------------------------
template <typename Float>
[[nodiscard]] constexpr auto norm_length(Float dx, Float min, Float max) noexcept -> Float {
  return dx / (max - min);
}

// -------------------------------------------------------------------------------------------------
template <typename Float>
[[nodiscard]] constexpr auto
to_pixel_coord(Float norm_value, size_t num_cells, size_t scale) noexcept -> size_t {
  return static_cast<size_t>(std::round(std::clamp(norm_value, Float{0}, Float{1}) *
                                        static_cast<Float>(num_cells * scale)));
}

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
[[nodiscard]] auto render_graph(Rd::Canvas& canvas,
                                Eigen::Index solution_idx,
                                const Zap::IO::IncCellReader<Float, DIM>& u_reader,
                                const Eigen::Vector<Float, DIM>& min,
                                const Eigen::Vector<Float, DIM>& max,
                                Rd::Box graph_box,
                                Rd::Box colorbar_box,
                                const Args& args) noexcept -> bool {
  const auto& cells = u_reader.cells();

  bool failed_rendering = false;
#ifdef ZAP_PARALLEL_RENDERING
#pragma omp parallel for
#endif  // ZAP_PARALLEL_RENDERING
  for (const auto& cell : cells) {
    if (cell.is_cartesian()) {
      const Rd::Box rect{
          .col   = to_pixel_coord(norm_point(cell.x_min, u_reader.x_min(), u_reader.x_max()),
                                u_reader.nx(),
                                args.scale),
          .row   = to_pixel_coord(norm_point(cell.y_min, u_reader.y_min(), u_reader.y_max()),
                                u_reader.ny(),
                                args.scale),
          .width = to_pixel_coord(
              norm_length(cell.dx, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), args.scale),
          .height = to_pixel_coord(
              norm_length(cell.dy, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), args.scale),
      };
      Rd::RGB c{};
      if (args.same_range) {
        const auto total_min = std::min_element(std::cbegin(min), std::cend(min));
        const auto total_max = std::max_element(std::cbegin(max), std::cend(max));
        c = Rd::float_to_rgb(cell.get_cartesian().value(solution_idx), *total_min, *total_max);
      } else {
        c = Rd::float_to_rgb(
            cell.get_cartesian().value(solution_idx), min(solution_idx), max(solution_idx));
      }

      if (!canvas.draw_rect(rect, graph_box, c, true)) {
        Igor::Warn("Could not draw cell");
        failed_rendering = true;
      }
    } else if (cell.is_cut()) {
      const auto& cut_value = cell.get_cut();

      // Left sub-cell
      {
        const auto points = Zap::CellBased::get_left_points<decltype(cell), Float>(cell);
        std::vector<Rd::Point> canvas_points(points.size());
        for (size_t i = 0; i < points.size(); ++i) {
          canvas_points[i] = Rd::Point{
              .col = to_pixel_coord(norm_point(points[i](0), u_reader.x_min(), u_reader.x_max()),
                                    u_reader.nx(),
                                    args.scale),
              .row = to_pixel_coord(norm_point(points[i](1), u_reader.y_min(), u_reader.y_max()),
                                    u_reader.ny(),
                                    args.scale),
          };
        }

        Rd::RGB c{};
        if (args.same_range) {
          const auto total_min = std::min_element(std::cbegin(min), std::cend(min));
          const auto total_max = std::max_element(std::cbegin(max), std::cend(max));
          c = Rd::float_to_rgb(cell.get_cut().left_value(solution_idx), *total_min, *total_max);
        } else {
          c = Rd::float_to_rgb(
              cell.get_cut().left_value(solution_idx), min(solution_idx), max(solution_idx));
        }

        if (!canvas.draw_polygon(canvas_points, graph_box, c, true)) {
          Igor::Warn("Could not draw left polygon.");
          failed_rendering = true;
        }
      }

      // Right sub-cell
      {
        const auto points = Zap::CellBased::get_right_points<decltype(cell), Float>(cell);
        std::vector<Rd::Point> canvas_points(points.size());
        for (size_t i = 0; i < points.size(); ++i) {
          canvas_points[i] = Rd::Point{
              .col = to_pixel_coord(norm_point(points[i](0), u_reader.x_min(), u_reader.x_max()),
                                    u_reader.nx(),
                                    args.scale),
              .row = to_pixel_coord(norm_point(points[i](1), u_reader.y_min(), u_reader.y_max()),
                                    u_reader.ny(),
                                    args.scale),
          };
        }

        Rd::RGB c{};
        if (args.same_range) {
          const auto total_min = std::min_element(std::cbegin(min), std::cend(min));
          const auto total_max = std::max_element(std::cbegin(max), std::cend(max));
          c = Rd::float_to_rgb(cell.get_cut().right_value(solution_idx), *total_min, *total_max);
        } else {
          c = Rd::float_to_rgb(
              cell.get_cut().right_value(solution_idx), min(solution_idx), max(solution_idx));
        }

        if (!canvas.draw_polygon(canvas_points, graph_box, c, true)) {
          Igor::Warn("Could not draw right polygon.");
          failed_rendering = true;
        }
      }

      // Draw shock curve
      const Rd::Point p1{
          .col = to_pixel_coord(norm_point(cut_value.x1_cut, u_reader.x_min(), u_reader.x_max()),
                                u_reader.nx(),
                                args.scale),
          .row = to_pixel_coord(norm_point(cut_value.y1_cut, u_reader.y_min(), u_reader.y_max()),
                                u_reader.ny(),
                                args.scale),
      };
      const Rd::Point p2{
          .col = to_pixel_coord(norm_point(cut_value.x2_cut, u_reader.x_min(), u_reader.x_max()),
                                u_reader.nx(),
                                args.scale),
          .row = to_pixel_coord(norm_point(cut_value.y2_cut, u_reader.y_min(), u_reader.y_max()),
                                u_reader.ny(),
                                args.scale),
      };

      if (!canvas.draw_line(p1, p2, graph_box, RED, true)) {
        Igor::Warn("Could not draw line");
        failed_rendering = true;
      }
    } else {
      Igor::Warn("Unknown cell type with variant index {}", cell.value.index());
      failed_rendering = true;
    }
  }

  if (!canvas.set_font_size(6)) { failed_rendering = true; }
  constexpr size_t NUM_COLORBAR_BLOCKS = 10;
  for (size_t i = 0; i < NUM_COLORBAR_BLOCKS; ++i) {
    const Float frac = static_cast<Float>(i) / static_cast<Float>(NUM_COLORBAR_BLOCKS - 1);
    Rd::RGB c{};
    Float value{};
    if (args.same_range) {
      const auto total_min = std::min_element(std::cbegin(min), std::cend(min));
      const auto total_max = std::max_element(std::cbegin(max), std::cend(max));
      value                = std::lerp(*total_min, *total_max, frac);
      c                    = Rd::float_to_rgb(value, *total_min, *total_max);
    } else {
      value = std::lerp(min(solution_idx), max(solution_idx), frac);
      c     = Rd::float_to_rgb(value, min(solution_idx), max(solution_idx));
    }

    const Rd::Box block = {
        .col    = 0,
        .row    = i * colorbar_box.height / NUM_COLORBAR_BLOCKS,
        .width  = colorbar_box.width / 4,
        .height = colorbar_box.height / NUM_COLORBAR_BLOCKS,
    };
    if (!canvas.draw_rect(block, colorbar_box, c, true)) { failed_rendering = true; }

    const auto str_value     = std::format("{:.6f}", value);
    const Rd::Box text_block = {
        .col = colorbar_box.col + block.width,
        .row = colorbar_box.row +
               (NUM_COLORBAR_BLOCKS - i - 1) * colorbar_box.height / NUM_COLORBAR_BLOCKS,
        .width  = colorbar_box.width - block.width,
        .height = colorbar_box.height / NUM_COLORBAR_BLOCKS,
    };

    if (!canvas.draw_text(str_value, text_block, true)) { failed_rendering = true; }
  }

  return !failed_rendering;
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto render_1d(const Args& args) noexcept -> bool {
  using Float          = double;
  constexpr size_t DIM = 1;
  const auto scale     = args.scale;

  try {
    Zap::IO::IncMatrixReader<Float> t_reader(args.t_input_file);
    Zap::IO::IncCellReader<Float, DIM> u_reader(args.u_input_file);

    constexpr size_t TEXT_HEIGHT    = 100UZ;
    constexpr size_t GRAPH_PADDING  = 25UZ;
    constexpr size_t COLORBAR_WIDTH = 200UZ;
    const size_t graph_width        = u_reader.nx() * scale;
    const size_t graph_height       = u_reader.ny() * scale;
    Rd::Canvas canvas(graph_width + 4 * GRAPH_PADDING + COLORBAR_WIDTH, graph_height + TEXT_HEIGHT);
    const Rd::Box text_box{
        .col    = 0,
        .row    = 0,
        .width  = canvas.width(),
        .height = TEXT_HEIGHT,
    };

    if (canvas.width() < 600) {
      Igor::Warn("Canvas is very narrow, redering might fail. Consider increasing scale with the "
                 "command line switch `--scale`.");
    }
    Igor::Info("Resolution: {}x{}", canvas.width(), canvas.height());

    const Rd::Box graph_box{
        .col    = GRAPH_PADDING,
        .row    = text_box.height,
        .width  = graph_width,
        .height = graph_height,
    };

    const Rd::Box colorbar_box{
        .col    = graph_box.width + 3 * GRAPH_PADDING,
        .row    = text_box.height,
        .width  = COLORBAR_WIDTH,
        .height = graph_height,
    };

    constexpr std::string_view font_file = "../assets/LinLibertine_RI.ttf";
    if (!canvas.load_font(font_file)) {
      Igor::Warn("Could not load font from file `{}`.", font_file);
      return false;
    }

    Igor::ScopeTimer timer{"Rendering"};

    const Rd::FFmpeg ffmpeg(canvas.width(), canvas.height(), args.output_file, args.fps);
    Eigen::Vector<Float, DIM> min{};
    Eigen::Vector<Float, DIM> max{};
    for (size_t iter = 0; true; ++iter) {
      const auto got_next_t = t_reader.read_next<false>();
      const auto got_next_u = u_reader.read_next<false>();
      if (!got_next_u && !got_next_t) {
        break;
      } else if (!got_next_u) {
        Igor::Warn("Got more t-entries than u-entries.");
        return false;
      } else if (!got_next_t) {
        Igor::Warn("Got more u-entries than t-entries.");
        return false;
      }

      canvas.clear(BLACK);

      if (!canvas.set_font_size(16)) { return false; }
      if (const std::string str = std::format("t = {:.6f}", t_reader(0, 0));
          !canvas.draw_text(str, text_box, true)) {
        Igor::Warn("Could not render string `{}` to canvas.", str);
        return false;
      }

      const auto& cells = u_reader.cells();
      assert(cells.size() > 0);
      if (iter == 0 || !args.keep_range) {
        std::tie(min, max) = min_max_cell_value<Float, DIM>(cells);
      }

      if (!render_graph<Float, DIM>(canvas, 0, u_reader, min, max, graph_box, colorbar_box, args)) {
        return false;
      }

      if (args.save_frame_images) {
        if (std::string filename = std::format("tmp_{:0>4}.ppm", iter); !canvas.to_ppm(filename)) {
          Igor::Warn("Could not write canvas to file `{}`: {}.", filename, std::strerror(errno));
        }
      }

      if (!canvas.to_raw_stream(ffmpeg.stream())) {
        Igor::Warn("Could not write canvas to FFmpeg.");
        return false;
      }
    }
  } catch (const std::exception& e) {
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC) << e.what() << '\n';
    return false;
  }
  return true;
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto render_2d(const Args& args) noexcept -> bool {
  using Float          = double;
  constexpr size_t DIM = 2;
  const auto scale     = args.scale;

  try {
    Zap::IO::IncMatrixReader<Float> t_reader(args.t_input_file);
    Zap::IO::IncCellReader<Float, DIM> u_reader(args.u_input_file);

    constexpr size_t TEXT_HEIGHT    = 100UZ;
    constexpr size_t GRAPH_PADDING  = 25UZ;
    constexpr size_t COLORBAR_WIDTH = 200UZ;
    const size_t graph_width        = u_reader.nx() * scale;
    const size_t graph_height       = u_reader.ny() * scale;
    Rd::Canvas canvas(2 * graph_width + 2 * COLORBAR_WIDTH + 8 * GRAPH_PADDING,
                      graph_height + TEXT_HEIGHT);
    const Rd::Box text_box{
        .col    = 0,
        .row    = 0,
        .width  = canvas.width(),
        .height = TEXT_HEIGHT,
    };

    Igor::Info("Resolution: {}x{}", canvas.width(), canvas.height());

    const Rd::Box left_graph_box{
        .col    = GRAPH_PADDING,
        .row    = text_box.height,
        .width  = graph_width,
        .height = graph_height,
    };
    const Rd::Box left_colorbar_box{
        .col    = left_graph_box.col + left_graph_box.width + GRAPH_PADDING,
        .row    = text_box.height,
        .width  = COLORBAR_WIDTH,
        .height = graph_height,
    };

    const Rd::Box right_graph_box{
        .col    = left_colorbar_box.col + left_colorbar_box.width + GRAPH_PADDING,
        .row    = text_box.height,
        .width  = graph_width,
        .height = graph_height,
    };
    const Rd::Box right_colorbar_box{
        .col    = right_graph_box.col + right_graph_box.width + GRAPH_PADDING,
        .row    = text_box.height,
        .width  = COLORBAR_WIDTH,
        .height = graph_height,
    };

    constexpr std::string_view font_file = "../assets/LinLibertine_RI.ttf";
    if (!canvas.load_font(font_file)) {
      Igor::Warn("Could not load font from file `{}`.", font_file);
      return false;
    }

    Igor::ScopeTimer timer{"Rendering"};

    const Rd::FFmpeg ffmpeg(canvas.width(), canvas.height(), args.output_file, args.fps);
    Eigen::Vector<Float, DIM> min{};
    Eigen::Vector<Float, DIM> max{};
    for (size_t iter = 0; true; ++iter) {
      const auto got_next_t = t_reader.read_next<false>();
      const auto got_next_u = u_reader.read_next<false>();
      if (!got_next_u && !got_next_t) {
        break;
      } else if (!got_next_u) {
        Igor::Warn("Got more t-entries than u-entries.");
        return false;
      } else if (!got_next_t) {
        Igor::Warn("Got more u-entries than t-entries.");
        return false;
      }

      canvas.clear(BLACK);

      if (!canvas.set_font_size(16)) { return false; }
      if (const std::string str = std::format("t = {:.6f}", t_reader(0, 0));
          !canvas.draw_text(str, text_box, true)) {
        Igor::Warn("Could not render string `{}` to canvas.", str);
        return false;
      }

      const auto& cells = u_reader.cells();
      assert(cells.size() > 0);
      if (iter == 0 || !args.keep_range) {
        std::tie(min, max) = min_max_cell_value<Float, DIM>(cells);
      }

      if (!render_graph<Float, DIM>(
              canvas, 0, u_reader, min, max, left_graph_box, left_colorbar_box, args)) {
        return false;
      }
      if (!render_graph<Float, DIM>(
              canvas, 1, u_reader, min, max, right_graph_box, right_colorbar_box, args)) {
        return false;
      }

      //       bool failed_drawing_rects = false;
      // #ifdef ZAP_PARALLEL_RENDERING
      // #pragma omp parallel for
      // #endif  // ZAP_PARALLEL_RENDERING
      //       for (const auto& cell : cells) {
      //         if (cell.is_cartesian()) {
      //           const Rd::Box rect{
      //               .col = to_pixel_coord(
      //                   norm_point(cell.x_min, u_reader.x_min(), u_reader.x_max()),
      //                   u_reader.nx(), scale),
      //               .row = to_pixel_coord(
      //                   norm_point(cell.y_min, u_reader.y_min(), u_reader.y_max()),
      //                   u_reader.ny(), scale),
      //               .width = to_pixel_coord(
      //                   norm_length(cell.dx, u_reader.x_min(), u_reader.x_max()), u_reader.nx(),
      //                   scale),
      //               .height = to_pixel_coord(
      //                   norm_length(cell.dy, u_reader.y_min(), u_reader.y_max()), u_reader.ny(),
      //                   scale),
      //           };
      //           Rd::RGB left_color{};
      //           Rd::RGB right_color{};
      //           if (args.same_range) {
      //             left_color = Rd::float_to_rgb(
      //                 cell.get_cartesian().value(0), std::min(min(0), min(1)), std::max(max(0),
      //                 max(1)));
      //             right_color = Rd::float_to_rgb(
      //                 cell.get_cartesian().value(1), std::min(min(0), min(1)), std::max(max(0),
      //                 max(1)));
      //           } else {
      //             left_color  = Rd::float_to_rgb(cell.get_cartesian().value(0), min(0), max(0));
      //             right_color = Rd::float_to_rgb(cell.get_cartesian().value(1), min(1), max(1));
      //           }
      //
      //           if (!canvas.draw_rect(rect, left_graph_box, left_color, true)) {
      //             Igor::Warn("Could not draw cell");
      //             failed_drawing_rects = true;
      //           }
      //           if (!canvas.draw_rect(rect, right_graph_box, right_color, true)) {
      //             Igor::Warn("Could not draw cell");
      //             failed_drawing_rects = true;
      //           }
      //         } else if (cell.is_cut()) {
      //           const auto& cut_value = cell.get_cut();
      //
      //           {
      //             const auto points = Zap::CellBased::get_left_points<decltype(cell),
      //             Float>(cell); std::vector<Rd::Point> canvas_points(points.size()); for (size_t
      //             i = 0; i < points.size(); ++i) {
      //               canvas_points[i] = Rd::Point{
      //                   .col =
      //                       to_pixel_coord(norm_point(points[i](0), u_reader.x_min(),
      //                       u_reader.x_max()),
      //                                      u_reader.nx(),
      //                                      scale),
      //                   .row =
      //                       to_pixel_coord(norm_point(points[i](1), u_reader.y_min(),
      //                       u_reader.y_max()),
      //                                      u_reader.ny(),
      //                                      scale),
      //               };
      //             }
      //
      //             Rd::RGB left_color{};
      //             Rd::RGB right_color{};
      //             if (args.same_range) {
      //               left_color = Rd::float_to_rgb(
      //                   cell.get_cut().left_value(0), std::min(min(0), min(1)), std::max(max(0),
      //                   max(1)));
      //               right_color = Rd::float_to_rgb(
      //                   cell.get_cut().left_value(1), std::min(min(0), min(1)), std::max(max(0),
      //                   max(1)));
      //             } else {
      //               left_color  = Rd::float_to_rgb(cell.get_cut().left_value(0), min(0), max(0));
      //               right_color = Rd::float_to_rgb(cell.get_cut().left_value(1), min(1), max(1));
      //             }
      //             if (!canvas.draw_polygon(canvas_points, left_graph_box, left_color, true)) {
      //               Igor::Warn("Could not draw left polygon.");
      //               failed_drawing_rects = true;
      //             }
      //             if (!canvas.draw_polygon(canvas_points, right_graph_box, right_color, true)) {
      //               Igor::Warn("Could not draw left polygon.");
      //               failed_drawing_rects = true;
      //             }
      //           }
      //
      //           {
      //             const auto points = Zap::CellBased::get_right_points<decltype(cell),
      //             Float>(cell); std::vector<Rd::Point> canvas_points(points.size()); for (size_t
      //             i = 0; i < points.size(); ++i) {
      //               canvas_points[i] = Rd::Point{
      //                   .col =
      //                       to_pixel_coord(norm_point(points[i](0), u_reader.x_min(),
      //                       u_reader.x_max()),
      //                                      u_reader.nx(),
      //                                      scale),
      //                   .row =
      //                       to_pixel_coord(norm_point(points[i](1), u_reader.y_min(),
      //                       u_reader.y_max()),
      //                                      u_reader.ny(),
      //                                      scale),
      //               };
      //             }
      //
      //             Rd::RGB left_color{};
      //             Rd::RGB right_color{};
      //             if (args.same_range) {
      //               left_color  = Rd::float_to_rgb(cell.get_cut().right_value(0),
      //                                             std::min(min(0), min(1)),
      //                                             std::max(max(0), max(1)));
      //               right_color = Rd::float_to_rgb(cell.get_cut().right_value(1),
      //                                              std::min(min(0), min(1)),
      //                                              std::max(max(0), max(1)));
      //             } else {
      //               left_color  = Rd::float_to_rgb(cell.get_cut().right_value(0), min(0),
      //               max(0)); right_color = Rd::float_to_rgb(cell.get_cut().right_value(1),
      //               min(1), max(1));
      //             }
      //             if (!canvas.draw_polygon(canvas_points, left_graph_box, left_color, true)) {
      //               Igor::Warn("Could not draw left polygon.");
      //               failed_drawing_rects = true;
      //             }
      //             if (!canvas.draw_polygon(canvas_points, right_graph_box, right_color, true)) {
      //               Igor::Warn("Could not draw left polygon.");
      //               failed_drawing_rects = true;
      //             }
      //           }
      //
      //           // Draw shock curve
      //           const Rd::Point p1{
      //               .col =
      //                   to_pixel_coord(norm_point(cut_value.x1_cut, u_reader.x_min(),
      //                   u_reader.x_max()),
      //                                  u_reader.nx(),
      //                                  scale),
      //               .row =
      //                   to_pixel_coord(norm_point(cut_value.y1_cut, u_reader.y_min(),
      //                   u_reader.y_max()),
      //                                  u_reader.ny(),
      //                                  scale),
      //           };
      //           const Rd::Point p2{
      //               .col =
      //                   to_pixel_coord(norm_point(cut_value.x2_cut, u_reader.x_min(),
      //                   u_reader.x_max()),
      //                                  u_reader.nx(),
      //                                  scale),
      //               .row =
      //                   to_pixel_coord(norm_point(cut_value.y2_cut, u_reader.y_min(),
      //                   u_reader.y_max()),
      //                                  u_reader.ny(),
      //                                  scale),
      //           };
      //
      //           if (!canvas.draw_line(p1, p2, left_graph_box, RED, true)) {
      //             Igor::Warn("Could not draw line");
      //             failed_drawing_rects = true;
      //           }
      //           if (!canvas.draw_line(p1, p2, right_graph_box, RED, true)) {
      //             Igor::Warn("Could not draw line");
      //             failed_drawing_rects = true;
      //           }
      //         } else {
      // #pragma omp critical
      //           {
      //             ffmpeg.~FFmpeg();
      //             Igor::Panic("Unknown cell type with variant index {}", cell.value.index());
      //           }
      //         }
      //       }
      // if (failed_drawing_rects) { return false; }

      if (args.save_frame_images) {
        if (std::string filename = std::format("tmp_{:0>4}.ppm", iter); !canvas.to_ppm(filename)) {
          Igor::Warn("Could not write canvas to file `{}`: {}.", filename, std::strerror(errno));
        }
      }

      if (!canvas.to_raw_stream(ffmpeg.stream())) {
        Igor::Warn("Could not write canvas to FFmpeg.");
        return false;
      }
    }
  } catch (const std::exception& e) {
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC) << e.what() << '\n';
    return false;
  }
  return true;
}

// -------------------------------------------------------------------------------------------------
auto main(int argc, char** argv) -> int {
  const auto args = parse_args(argc, argv);
  if (!args.has_value()) { return 1; }

  Igor::Info("scale        = {}", args->scale);
  Igor::Info("fps          = {}", args->fps);
  Igor::Info("keep-range   = {}", args->keep_range);
  Igor::Info("same-range   = {}", args->same_range);
  Igor::Info("dimension    = {}", !args->two_dim ? "1d" : "2d");
  Igor::Info("u input file = {}", args->u_input_file);
  Igor::Info("t input file = {}", args->t_input_file);
  Igor::Info("output file  = {}", args->output_file);

  if (!args->two_dim) {
    return render_1d(*args) ? 0 : 1;
  } else {
    return render_2d(*args) ? 0 : 1;
  }
}
