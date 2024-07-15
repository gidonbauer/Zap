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
template <Igor::detail::Level level>
constexpr void usage(std::string_view prog) {
  std::cerr << std::format("{}Usage: {} [--scale s] [--keep-range] [--save-frame-img] "
                           "<u input file> <t input file> <output file>\n",
                           Igor::detail::level_repr(level),
                           prog);
};

// -------------------------------------------------------------------------------------------------
struct Args {
  std::string u_input_file;
  std::string t_input_file;
  std::string output_file;
  size_t scale           = 1;
  bool keep_range        = false;
  bool save_frame_images = false;
};

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
    } else if (*argv == "--keep-range"sv) {
      args.keep_range = true;
    } else if (*argv == "--save-frame-img"sv) {
      args.save_frame_images = true;
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
auto main(int argc, char** argv) -> int {
  const auto args = parse_args(argc, argv);
  if (!args.has_value()) {
    return 1;
  }

  Igor::Info("scale        = {}", args->scale);
  Igor::Info("keep-range   = {}", args->keep_range);
  Igor::Info("u input file = {}", args->u_input_file);
  Igor::Info("t input file = {}", args->t_input_file);
  Igor::Info("output file  = {}", args->output_file);

  // -----------------------------------------------------------------------------------------------

  using Float          = double;
  constexpr size_t DIM = 1;  // TODO: Implement renderer for 2D data
  const auto scale     = args->scale;

  try {
    Zap::IO::IncMatrixReader<Float> t_reader(args->t_input_file);
    Zap::IO::IncCellReader<Float, DIM> u_reader(args->u_input_file);

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

    const Rd::FFmpeg ffmpeg(canvas.width(), canvas.height(), args->output_file);
    Eigen::Vector<Float, DIM> min{};
    Eigen::Vector<Float, DIM> max{};
    for (size_t iter = 0; true; ++iter) {
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

      canvas.clear(BLACK);

      if (const std::string str = std::format("t = {:.6f}", t_reader(0, 0));
          !canvas.draw_text(str, text_box, true)) {
        Igor::Warn("Could not render string `{}` to canvas.", str);
        return 1;
      }

      const auto& cells = u_reader.cells();
      assert(cells.size() > 0);
      if (iter == 0 || !args->keep_range) {
        std::tie(min, max) = min_max_cell_value<Float, DIM>(cells);
      }

      bool failed_drawing_rects = false;
#ifdef ZAP_PARALLEL_RENDERING
#pragma omp parallel for
#endif  // ZAP_PARALLEL_RENDERING
      for (const auto& cell : cells) {
        if (cell.is_cartesian()) {
          const Rd::Box rect{
              .col = to_pixel_coord(
                  norm_point(cell.x_min, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
              .row = to_pixel_coord(
                  norm_point(cell.y_min, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
              .width = to_pixel_coord(
                  norm_length(cell.dx, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
              .height = to_pixel_coord(
                  norm_length(cell.dy, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
          };
          const auto c = Rd::float_to_rgb(cell.get_cartesian().value(0), min(0), max(0));

          if (!canvas.draw_rect(rect, graph_box, c, true)) {
            Igor::Warn("Could not draw cell");
            failed_drawing_rects = true;
          }
        } else if (cell.is_cut()) {
          const auto& cut_value = cell.get_cut();

          const auto left_points = Zap::CellBased::get_left_points<decltype(cell), Float>(cell);
          std::vector<Rd::Point> left_canvas_points(left_points.size());
          for (size_t i = 0; i < left_points.size(); ++i) {
            left_canvas_points[i] = Rd::Point{
                .col = to_pixel_coord(
                    norm_point(left_points[i](0), u_reader.x_min(), u_reader.x_max()),
                    u_reader.nx(),
                    scale),
                .row = to_pixel_coord(
                    norm_point(left_points[i](1), u_reader.y_min(), u_reader.y_max()),
                    u_reader.ny(),
                    scale),
            };
          }
          if (!canvas.draw_polygon(left_canvas_points,
                                   graph_box,
                                   Rd::float_to_rgb(cut_value.left_value(0), min(0), max(0)),
                                   true)) {
            Igor::Warn("Could not draw left polygon.");
            failed_drawing_rects = true;
          }

          const auto right_points = Zap::CellBased::get_right_points<decltype(cell), Float>(cell);
          std::vector<Rd::Point> right_canvas_points(right_points.size());
          for (size_t i = 0; i < right_points.size(); ++i) {
            right_canvas_points[i] = Rd::Point{
                .col = to_pixel_coord(
                    norm_point(right_points[i](0), u_reader.x_min(), u_reader.x_max()),
                    u_reader.nx(),
                    scale),
                .row = to_pixel_coord(
                    norm_point(right_points[i](1), u_reader.y_min(), u_reader.y_max()),
                    u_reader.ny(),
                    scale),
            };
          }
          if (!canvas.draw_polygon(right_canvas_points,
                                   graph_box,
                                   Rd::float_to_rgb(cut_value.right_value(0), min(0), max(0)),
                                   true)) {
            Igor::Warn("Could not draw right polygon.");
            failed_drawing_rects = true;
          }

          // Draw shock curve
          const Rd::Point p1{
              .col =
                  to_pixel_coord(norm_point(cut_value.x1_cut, u_reader.x_min(), u_reader.x_max()),
                                 u_reader.nx(),
                                 scale),
              .row =
                  to_pixel_coord(norm_point(cut_value.y1_cut, u_reader.y_min(), u_reader.y_max()),
                                 u_reader.ny(),
                                 scale),
          };
          const Rd::Point p2{
              .col =
                  to_pixel_coord(norm_point(cut_value.x2_cut, u_reader.x_min(), u_reader.x_max()),
                                 u_reader.nx(),
                                 scale),
              .row =
                  to_pixel_coord(norm_point(cut_value.y2_cut, u_reader.y_min(), u_reader.y_max()),
                                 u_reader.ny(),
                                 scale),
          };

          if (!canvas.draw_line(p1, p2, graph_box, RED, true)) {
            Igor::Warn("Could not draw line");
            failed_drawing_rects = true;
          }
        } else {
#pragma omp critical
          {
            ffmpeg.~FFmpeg();
            Igor::Panic("Unknown cell type with variant index {}", cell.value.index());
          }
        }
      }
      if (failed_drawing_rects) {
        return 1;
      }

      if (args->save_frame_images) {
        if (std::string filename = std::format("tmp_{}.ppm", iter); !canvas.to_ppm(filename)) {
          Igor::Warn("Could not write canvas to file `{}`: {}.", filename, std::strerror(errno));
        }
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
