#include <cmath>
#include <optional>

#include "CellBased/Cell.hpp"
#include "IO/IncCellReader.hpp"
#include "IO/IncMatrixReader.hpp"
#include "Renderer/Canvas.hpp"
#include "Renderer/FFmpeg.hpp"
#include "Renderer/Modules.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Timer.hpp"

#define ZAP_PARALLEL_RENDERING

// TODO: Fix rendering of cut cells, sometimes cell goes out of bounds
// TODO: Display total mass -> can see if it is conservative

namespace Rd = Zap::Renderer;

[[maybe_unused]] constexpr Rd::RGB BLACK  = {.r = 0x00, .g = 0x00, .b = 0x00};
[[maybe_unused]] constexpr Rd::RGB RED    = {.r = 0xFF, .g = 0x00, .b = 0x00};
[[maybe_unused]] constexpr Rd::RGB PINK   = {.r = 0xFF, .g = 0x00, .b = 0xFF};
[[maybe_unused]] constexpr Rd::RGB ORANGE = {.r = 0xFF, .g = 0xA5, .b = 0x00};

// -------------------------------------------------------------------------------------------------
struct Args {
  enum ImageFormat { NO_SAVE, PPM, JPEG };

  std::string u_input_file;
  std::string t_input_file;
  std::string output_file;
  size_t scale                  = 1;
  bool keep_range               = false;
  bool same_range               = false;
  ImageFormat save_frame_images = NO_SAVE;
  size_t fps                    = 60;
  bool two_dim                  = false;
  bool no_mass                  = false;
};

// -------------------------------------------------------------------------------------------------
template <Igor::detail::Level level>
constexpr void usage(std::string_view prog) {
  Args args{};

  std::cerr << Igor::detail::level_repr(level) << "Usage: " << prog
            << " [--scale s] [--keep-range] [--same-range] [--save-frame-img <format>] [--fps f] "
               "[--2d] [--no-mass] <u input file> <t input file> <output file>\n";
  std::cerr << "\t--scale           Number of pixels for single cell, default is " << args.scale
            << '\n';
  std::cerr << "\t--keep-range      Keep min and max value of first frame for colormap, default is "
            << std::boolalpha << args.keep_range << '\n';
  std::cerr << "\t--same-range      Same min and max value for both plots in 2d, default is "
            << std::boolalpha << args.same_range << '\n';
  std::cerr << "\t--save-frame-img  Save each individual frame as .ppm or .jpeg image, default is "
               "to not save the frames\n";
  std::cerr << "\t--fps             Frames per second for rendered video, default is " << args.fps
            << '\n';
  std::cerr
      << "\t--2d              Render solution of 2-dimensional solution (u in R^2), default is "
      << std::boolalpha << args.two_dim << '\n';
  std::cerr << "\t--no-mass         Do not display the mass ||u||, this flag saves a lot of "
               "computation time, default is "
            << args.no_mass << '\n';
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
        return std::nullopt;
      }
    } else if (*argv == "--keep-range"sv) {
      args.keep_range = true;
    } else if (*argv == "--same-range"sv) {
      args.same_range = true;
    } else if (*argv == "--save-frame-img"sv) {
      argc -= 1;
      argv += 1;  // NOLINT
      if (argc <= 0) {
        usage<Igor::detail::Level::WARN>(prog);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN)
                  << "Did not provide format to save the frame images.\n";
        return std::nullopt;
      }
      if (*argv == "ppm"sv) {
        args.save_frame_images = Args::PPM;
      } else if (*argv == "jpeg"sv || *argv == "jpg"sv) {
        args.save_frame_images = Args::JPEG;
      } else {
        usage<Igor::detail::Level::WARN>(prog);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::WARN) << "Invalid format `"
                  << *argv
                  << "` to save the frames as image. Available formats are `ppm` and `jpeg`\n";
        return std::nullopt;
      }
    } else if (*argv == "--2d"sv) {
      args.two_dim = true;
    } else if (*argv == "--no-mass"sv) {
      args.no_mass = true;
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
[[nodiscard]] auto
min_max_cell_value(const std::vector<Zap::IO::ReducedCell<Float, DIM>>& cells) noexcept
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
template <typename ResultType, typename FUNC, typename Float>
void simpsons_rule_2d(FUNC f,
                      ResultType& res,
                      Float x_min,
                      Float x_max,
                      Float y_min,
                      Float y_max,
                      size_t nx,
                      size_t ny) noexcept {
  auto get_w = [](size_t idx, size_t end_idx) -> Float {
    if (idx == 0 || idx == end_idx) {
      return 1;
    } else {
      if (idx % 2 == 1) {
        return 4;
      } else {
        return 2;
      }
    }
  };

  if (nx % 2 == 1) { nx += 1; }
  if (ny % 2 == 1) { ny += 1; }

  const Float dx = (x_max - x_min) / static_cast<Float>(nx);
  const Float dy = (y_max - y_min) / static_cast<Float>(ny);

#pragma omp parallel for reduction(+ : res)
  for (size_t i = 0; i < (nx + 1) * (ny + 1); ++i) {
    const size_t ix = i / (nx + 1);
    const size_t iy = i % (nx + 1);
    const Float x   = x_min + static_cast<Float>(ix) * dx;
    const Float y   = y_min + static_cast<Float>(iy) * dy;
    const Float wx  = get_w(ix, nx);
    const Float wy  = get_w(iy, ny);
    res += wx * wy * f(x, y);
  }
  res *= dx * dy / 9;
}

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto
calculate_mass(const std::variant<Zap::IO::IncCellReader<Float, DIM>,
                                  Zap::IO::IncMatrixReader<Float>>& u_reader,
               size_t nx = 1000,
               size_t ny = 1000) noexcept -> Eigen::Vector<Float, DIM> {
  const bool is_cell_reader = std::holds_alternative<Zap::IO::IncCellReader<Float, DIM>>(u_reader);

  if (is_cell_reader) {
    const auto& cell_reader = std::get<Zap::IO::IncCellReader<Float, DIM>>(u_reader);
    const Float x_min       = cell_reader.x_min();
    const Float x_max       = cell_reader.x_max();
    const Float y_min       = cell_reader.y_min();
    const Float y_max       = cell_reader.y_max();

    auto eval_at = [&](Float x, Float y) -> Eigen::Vector<Float, DIM> {
      const auto& cells = cell_reader.cells();

      auto cell_contains_xy = [x, y](const auto& cell) {
        constexpr Float eps = std::is_same_v<std::remove_cvref_t<Float>, float> ? 1e-6f : 1e-8;
        return x - cell.x_min() >= -eps && x - (cell.x_min() + cell.dx()) <= eps &&
               y - cell.y_min() >= -eps && y - (cell.y_min() + cell.dy()) <= eps;
      };

      const auto cell_it = std::find_if(std::cbegin(cells), std::cend(cells), cell_contains_xy);
      IGOR_ASSERT(cell_it != std::cend(cells),
                  "Point {{ .x={}, .y={} }} is not in grid [{}, {}]x[{}, {}]",
                  x,
                  y,
                  x_min,
                  x_max,
                  y_min,
                  y_max);

      if (cell_it->is_cartesian()) {
        return cell_it->get_cartesian().value;
      } else if (cell_it->is_cut()) {
        const auto point_in_left = Zap::CellBased::Geometry::Polygon(cell_it->get_left_points())
                                       .point_in_polygon({.x = x, .y = y});
        const auto point_in_right = Zap::CellBased::Geometry::Polygon(cell_it->get_right_points())
                                        .point_in_polygon({.x = x, .y = y});

        if (point_in_left) {
          return cell_it->get_cut().left_value;
        } else if (point_in_right) {
          return cell_it->get_cut().right_value;
        } else {
          return (cell_it->get_cut().left_value + cell_it->get_cut().right_value) / 2;
        }
      } else {
        Igor::Panic("Unknown cell type with variant index {}", cell_it->value.index());
        std::unreachable();
      }
    };

    Eigen::Vector<Float, DIM> res = Eigen::Vector<Float, DIM>::Zero();
    simpsons_rule_2d(eval_at, res, x_min, x_max, y_min, y_max, nx, ny);
    return res;
  } else {
    IGOR_ASSERT(DIM == 1,
                "The matrix reader only supports one-dimensional output and not {}-dimensional, "
                "i.e. u : R^2 -> R.",
                DIM);

    const auto& mat_reader = std::get<Zap::IO::IncMatrixReader<Float>>(u_reader);
    const auto x_min       = static_cast<Float>(0);
    const auto x_max       = static_cast<Float>(mat_reader.cols());
    const auto y_min       = static_cast<Float>(0);
    const auto y_max       = static_cast<Float>(mat_reader.rows());

    auto eval_at = [&](Float x, Float y) -> Eigen::Vector<Float, DIM> {
      constexpr Float eps = std::is_same_v<std::remove_cvref_t<Float>, float> ? 1e-6f : 1e-8;
      const int64_t xi    = [&] {
        if (std::abs(x - static_cast<Float>(mat_reader.cols())) <= eps) {
          return mat_reader.cols() - 1;
        } else {
          return static_cast<int64_t>(std::floor(x));
        }
      }();
      const int64_t yi = [&] {
        if (std::abs(y - static_cast<Float>(mat_reader.rows())) <= eps) {
          return mat_reader.rows() - 1;
        } else {
          return static_cast<int64_t>(std::floor(y));
        }
      }();

      IGOR_ASSERT(xi >= 0 && xi < mat_reader.cols(),
                  "xi={} is out of bounds must be in [0, {})",
                  xi,
                  mat_reader.cols());
      IGOR_ASSERT(yi >= 0 && yi < mat_reader.rows(),
                  "yi={} is out of bounds must be in [0, {})",
                  yi,
                  mat_reader.rows());
      return Eigen::Vector<Float, DIM>{mat_reader(yi, xi)};
    };

    Eigen::Vector<Float, DIM> res = Eigen::Vector<Float, DIM>::Zero();
    simpsons_rule_2d(eval_at, res, x_min, x_max, y_min, y_max, nx, ny);
    return res / (mat_reader.rows() * mat_reader.cols());
  }
  Igor::Panic("Unreachable");
  std::unreachable();
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto render_1d(const Args& args) noexcept -> bool {
  using Float          = double;
  constexpr size_t DIM = 1;
  const auto scale     = args.scale;

  try {
    using CellReader   = Zap::IO::IncCellReader<Float, DIM>;
    using MatrixReader = Zap::IO::IncMatrixReader<Float>;
    using Readers      = std::variant<CellReader, MatrixReader>;

    MatrixReader t_reader(args.t_input_file);
    auto u_reader = [&args]() -> Readers {
      if (args.u_input_file.ends_with(".grid")) {
        return CellReader(args.u_input_file);
      } else if (args.u_input_file.ends_with(".mat")) {
        return MatrixReader(args.u_input_file);
      } else {
        Igor::Panic("Unknown file format, use file extension `.grid` for cell-based solution or "
                    "`.mat` for matrix-based solution. Filename is {}",
                    args.u_input_file);
        std::unreachable();
      }
    }();

    constexpr size_t TEXT_HEIGHT    = 100UZ;
    constexpr size_t GRAPH_PADDING  = 25UZ;
    constexpr size_t COLORBAR_WIDTH = 200UZ;
    size_t graph_width              = 0;
    size_t graph_height             = 0;
    if (std::holds_alternative<CellReader>(u_reader)) {
      graph_width  = std::get<CellReader>(u_reader).nx() * scale;
      graph_height = std::get<CellReader>(u_reader).ny() * scale;
    } else {
      const auto n_rows = std::get<MatrixReader>(u_reader).rows();
      const auto n_cols = std::get<MatrixReader>(u_reader).cols();
      assert(n_rows > 0);
      assert(n_cols > 0);

      graph_width  = static_cast<size_t>(n_rows) * scale;
      graph_height = static_cast<size_t>(n_cols) * scale;
    }

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

    const Rd::FFmpeg ffmpeg(canvas.width(),
                            canvas.height(),
                            args.output_file,
                            args.fps,
                            "100000k",
                            args.save_frame_images != Args::NO_SAVE);
    Eigen::Vector<Float, DIM> min{};
    Eigen::Vector<Float, DIM> max{};
    for (size_t iter = 0; true; ++iter) {
      const auto got_next_t = t_reader.read_next<false>();
      const auto got_next_u = [&u_reader] {
        if (std::holds_alternative<CellReader>(u_reader)) {
          return std::get<CellReader>(u_reader).read_next<false>();
        } else {
          return std::get<MatrixReader>(u_reader).read_next<false>();
        }
      }();

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
      if (args.no_mass) {
        if (const std::string str = std::format("t = {:.6f}", t_reader(0, 0));
            !canvas.draw_text(str, text_box, true)) {
          Igor::Warn("Could not render string `{}` to canvas.", str);
          return false;
        }
      } else {
        if (const std::string str = std::format("t = {:.6f}, ||u|| = {:.6f}",
                                                t_reader(0, 0),
                                                calculate_mass<Float, DIM>(u_reader)(0));
            !canvas.draw_text(str, text_box, true)) {
          Igor::Warn("Could not render string `{}` to canvas.", str);
          return false;
        }
      }

      if (std::holds_alternative<CellReader>(u_reader)) {
        const auto& cells = std::get<CellReader>(u_reader).cells();
        assert(cells.size() > 0);
        if (iter == 0 || !args.keep_range) {
          std::tie(min, max) = min_max_cell_value<Float, DIM>(cells);
        }
      } else {
        const auto& data = std::get<MatrixReader>(u_reader).data();
        if (iter == 0 || !args.keep_range) {
          const auto min_max = std::minmax_element(std::cbegin(data), std::cend(data));
          min(0)             = *min_max.first;
          max(0)             = *min_max.second;
        }
      }

      if (std::holds_alternative<CellReader>(u_reader)) {
        if (!Rd::render_graph<Float, DIM>(canvas,
                                          std::get<CellReader>(u_reader),
                                          0,
                                          min,
                                          max,
                                          graph_box,
                                          args.scale,
                                          args.same_range)) {
          // return false;
        }
      } else {
        if (!Rd::render_graph<Float, DIM>(canvas,
                                          std::get<MatrixReader>(u_reader),
                                          0,
                                          min,
                                          max,
                                          graph_box,
                                          args.scale,
                                          args.same_range)) {
          // return false;
        }
      }

      constexpr size_t NUM_COLORBAR_BLOCKS = 10;
      if (!Rd::render_color_bar<Float, DIM>(
              canvas, colorbar_box, 0, min, max, args.same_range, NUM_COLORBAR_BLOCKS)) {
        // return false;
      }

      switch (args.save_frame_images) {
        case Args::PPM:
          if (std::string filename = std::format("tmp_{:0>4}.ppm", iter);
              !canvas.to_ppm(filename)) {
            Igor::Warn("Could not write canvas to file `{}`: {}.", filename, std::strerror(errno));
          }
          break;
        case Args::JPEG:
          if (std::string filename = std::format("tmp_{:0>4}.jpeg", iter);
              !canvas.to_jpeg(filename)) {
            Igor::Warn("Could not write canvas to file `{}`: {}.", filename, std::strerror(errno));
          }
          break;
        case Args::NO_SAVE:
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

      if (!Rd::render_graph<Float, DIM>(
              canvas, u_reader, 0, min, max, left_graph_box, args.scale, args.same_range)) {
        // return false;
      }
      if (!Rd::render_graph<Float, DIM>(
              canvas, u_reader, 1, min, max, right_graph_box, args.scale, args.same_range)) {
        // return false;
      }

      constexpr size_t NUM_COLORBAR_BLOCKS = 10;
      if (!Rd::render_color_bar<Float, DIM>(
              canvas, left_colorbar_box, 0, min, max, args.same_range, NUM_COLORBAR_BLOCKS)) {
        // return false;
      }
      if (!Rd::render_color_bar<Float, DIM>(
              canvas, right_colorbar_box, 1, min, max, args.same_range, NUM_COLORBAR_BLOCKS)) {
        // return false;
      }

      switch (args.save_frame_images) {
        case Args::PPM:
          if (std::string filename = std::format("tmp_{:0>4}.ppm", iter);
              !canvas.to_ppm(filename)) {
            Igor::Warn("Could not write canvas to file `{}`: {}.", filename, std::strerror(errno));
          }
          break;
        case Args::JPEG:
          if (std::string filename = std::format("tmp_{:0>4}.jpeg", iter);
              !canvas.to_jpeg(filename)) {
            Igor::Warn("Could not write canvas to file `{}`: {}.", filename, std::strerror(errno));
          }
          break;
        case Args::NO_SAVE:
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
