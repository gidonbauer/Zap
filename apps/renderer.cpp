#include <optional>

#include "IO/IncCellReader.hpp"
#include "IO/IncMatrixReader.hpp"

#include "Renderer/Canvas.hpp"
#include "Renderer/FFmpeg.hpp"
#include "Renderer/Modules.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Timer.hpp"

#include "Common.hpp"

#define ZAP_PARALLEL_RENDERING

using namespace std::string_view_literals;
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

  size_t scale = 1;
  size_t fps   = 60;
  bool no_mass = false;

  bool keep_range = false;

  std::string v_input_file;

  bool symmetry_error = false;

  ImageFormat save_frame_images = NO_SAVE;
};

// -------------------------------------------------------------------------------------------------
template <Igor::detail::Level level>
constexpr void usage(std::string_view prog) {
  Args args{};

  std::cerr << Igor::detail::level_repr(level) << "Usage: " << prog
            << " [--scale s] [--fps f] [--no-mass] [--keep-range] [--der <v input file>] "
               "[--symmetry-error] [--save-frame-img <format>] "
               "<u input file> <t input file> <output file>\n";
  std::cerr << "\t--scale           Number of pixels for single cell, default is " << args.scale
            << '\n';
  std::cerr << "\t--fps             Frames per second for rendered video, default is " << args.fps
            << '\n';
  std::cerr << "\t--no-mass         Do not display the mass ||u||, this flag saves some "
               "computation time, default is "
            << std::boolalpha << args.no_mass << '\n';
  std::cerr << "\t--keep-range      Keep min and max value of first frame for colormap, default is "
            << std::boolalpha << args.keep_range << '\n';
  std::cerr
      << "\t--der             Input file for the derivative of u, default is no derivative.\n";
  std::cerr << "\t--symmetry-error  Show the symmetry error with symmetry axis [1, 1], default is "
            << std::boolalpha << args.symmetry_error << '\n';
  std::cerr << "\t--save-frame-img  Save each individual frame as .ppm or .jpeg image, default is "
               "to not save the frames\n";
  std::cerr << "\tu input file      Solution field for each iteration\n";
  std::cerr << "\tt input file      Time for each iteration\n";
  std::cerr << "\toutput file       File for the rendered video\n";
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto parse_args(int argc, char** argv) noexcept -> std::optional<Args> {
  Args args{};

  const auto prog = next_arg(argc, argv);
  IGOR_ASSERT(prog.has_value(), "Missing program name, this should never happen.");

  for (auto arg = next_arg(argc, argv); arg.has_value(); arg = next_arg(argc, argv)) {
    // ImageFormat save_frame_images = NO_SAVE;
    if (*arg == "--scale"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide value for scale.\n";
        return std::nullopt;
      }
      args.scale = parse_size_t(*arg);
      if (args.scale == 0) {
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Scale must be larger than zero but is " << args.scale << '\n';
        return std::nullopt;
      }
    } else if (*arg == "--fps"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide value for fps.\n";
        return std::nullopt;
      }
      args.fps = parse_size_t(*arg);
      if (args.fps == 0) {
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "FPS must be larger than zero but is " << args.fps << '\n';
        return std::nullopt;
      }
    } else if (*arg == "--no-mass"sv) {
      args.no_mass = true;
    } else if (*arg == "--keep-range"sv) {
      args.keep_range = true;
    } else if (*arg == "--der"sv || *arg == "--derivative"sv || *arg == "-v"sv || *arg == "--v"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide input file for derivative data.\n";
        return std::nullopt;
      }
      args.v_input_file = *arg;
    } else if (*arg == "--symmetry-error"sv) {
      args.symmetry_error = true;
    } else if (*arg == "--save-frame-img"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide format to save the frame images.\n";
        return std::nullopt;
      }
      if (*arg == "ppm"sv) {
        args.save_frame_images = Args::PPM;
      } else if (*arg == "jpeg"sv || *arg == "jpg"sv) {
        args.save_frame_images = Args::JPEG;
      } else {
        usage<Igor::detail::Level::PANIC>(*prog);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC) << "Invalid format `"
                  << *arg
                  << "` to save the frames as image. Available formats are `ppm` and `jpeg`\n";
        return std::nullopt;
      }
    } else if (*arg == "--help"sv || *arg == "-h"sv) {
      usage<Igor::detail::Level::INFO>(*prog);
      return std::nullopt;
    } else if (args.u_input_file.empty()) {
      args.u_input_file = *arg;
    } else if (args.t_input_file.empty()) {
      args.t_input_file = *arg;
    } else if (args.output_file.empty()) {
      args.output_file = *arg;
    } else {
      usage<Igor::detail::Level::PANIC>(*prog);
      std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                << "Provided too many arguments or invalid switch.\n";
      return std::nullopt;
    }
  }

  if (args.u_input_file.empty()) {
    usage<Igor::detail::Level::PANIC>(*prog);
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
              << "Did not provide u input file.\n";
    return std::nullopt;
  }
  if (args.t_input_file.empty()) {
    usage<Igor::detail::Level::PANIC>(*prog);
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
              << "Did not provide t input file.\n";
    return std::nullopt;
  }
  if (args.output_file.empty()) {
    usage<Igor::detail::Level::PANIC>(*prog);
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
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
template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto calculate_mass(
    const std::variant<Zap::IO::IncCellReader<Float, DIM>, Zap::IO::IncMatrixReader<Float>>&
        u_reader) noexcept -> Eigen::Vector<Float, DIM> {
  // Note: We are doing a piecewise constant reconstruction of u, therefore the rectangle rule is
  //       appropriate
  const bool is_cell_reader = std::holds_alternative<Zap::IO::IncCellReader<Float, DIM>>(u_reader);

  if (is_cell_reader) {
    const auto& cell_reader = std::get<Zap::IO::IncCellReader<Float, DIM>>(u_reader);
    const Float x_min       = cell_reader.x_min();
    const Float x_max       = cell_reader.x_max();
    const Float y_min       = cell_reader.y_min();
    const Float y_max       = cell_reader.y_max();

    Eigen::Vector<Float, DIM> mass = Eigen::Vector<Float, DIM>::Zero();
    for (const auto& cell : cell_reader.cells()) {
      assert(cell.is_cartesian() || cell.is_cut());
      if (cell.is_cartesian()) {
        mass += cell.get_cartesian().value * cell.dx() * cell.dy();
      } else {
        const auto left_area  = Zap::CellBased::Geometry::Polygon(cell.get_left_points()).area();
        const auto right_area = Zap::CellBased::Geometry::Polygon(cell.get_right_points()).area();
        mass += cell.get_cut().left_value * left_area + cell.get_cut().right_value * right_area;
      }
    }

    // Scale down to 1x1 domain to compare against the matrix-based solution
    return mass / ((x_max - x_min) * (y_max - y_min));
  } else {
    IGOR_ASSERT(DIM == 1,
                "The matrix reader only supports one-dimensional output and not {}-dimensional, "
                "i.e. u : R^2 -> R.",
                DIM);

    const auto& mat_reader = std::get<Zap::IO::IncMatrixReader<Float>>(u_reader);

    Eigen::Vector<Float, DIM> mass = Eigen::Vector<Float, DIM>::Zero();
    for (const auto& value : mat_reader.data()) {
      mass(0) += value;
    }

    // Scale down to 1x1 domain
    return mass / (mat_reader.rows() * mat_reader.cols());
  }
  Igor::Panic("Unreachable");
  std::unreachable();
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto render_u(const Args& args) noexcept -> bool {
  using Float          = double;
  constexpr size_t DIM = 1;
  const auto scale     = args.scale;

  try {
    using CellReader   = Zap::IO::IncCellReader<Float, DIM>;
    using MatrixReader = Zap::IO::IncMatrixReader<Float>;
    using Readers      = std::variant<CellReader, MatrixReader>;

    MatrixReader t_reader(args.t_input_file);
    Igor::Info("#frames      = {}", t_reader.num_elem());

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

    const size_t canvas_width  = graph_width + 4 * GRAPH_PADDING + COLORBAR_WIDTH;
    const size_t canvas_height = graph_height + TEXT_HEIGHT;
    Rd::Canvas canvas(canvas_width + canvas_width % 2, canvas_height + canvas_height % 2);
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
    Igor::on_death.emplace_back([child = ffmpeg.child()]() {
      if (kill(child, SIGKILL) == -1) {
        Igor::Warn("Could not kill FFmpeg: {}", std::strerror(errno));
      }
    });

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

      if (args.symmetry_error) {
        std::visit([](auto& reader) { reader.calc_symmetry_error(); }, u_reader);
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

      std::visit(
          [&](auto& reader) {
            bool _ =
                Rd::render_graph<Float, DIM>(canvas, reader, 0, min, max, graph_box, args.scale);
          },
          u_reader);

      constexpr size_t NUM_COLORBAR_BLOCKS = 10;
      if (!Rd::render_color_bar<Float, DIM>(
              canvas, colorbar_box, 0, min, max, NUM_COLORBAR_BLOCKS)) {
        // return false;
      }

      switch (args.save_frame_images) {
        case Args::PPM:
          if (std::string filename = std::format("frame_{:0>4}.ppm", iter);
              !canvas.to_ppm(filename)) {
            Igor::Warn("Could not write canvas to file `{}`: {}.", filename, std::strerror(errno));
          }
          break;
        case Args::JPEG:
          if (std::string filename = std::format("frame_{:0>4}.jpeg", iter);
              !canvas.to_jpeg(filename)) {
            Igor::Warn("Could not write canvas to file `{}`: {}.", filename, std::strerror(errno));
          }
          break;
        case Args::NO_SAVE: break;
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
[[nodiscard]] auto render_uv(const Args& args) noexcept -> bool {
  using Float          = double;
  constexpr size_t DIM = 1;
  const auto scale     = args.scale;

  try {
    using CellReader   = Zap::IO::IncCellReader<Float, DIM>;
    using MatrixReader = Zap::IO::IncMatrixReader<Float>;
    using Readers      = std::variant<CellReader, MatrixReader>;

    MatrixReader t_reader(args.t_input_file);
    Igor::Info("#frames      = {}", t_reader.num_elem());

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

    auto v_reader = [&args]() -> Readers {
      if (args.v_input_file.ends_with(".grid")) {
        return CellReader(args.v_input_file);
      } else if (args.v_input_file.ends_with(".mat")) {
        return MatrixReader(args.v_input_file);
      } else {
        Igor::Panic("Unknown file format, use file extension `.grid` for cell-based solution or "
                    "`.mat` for matrix-based solution. Filename is {}",
                    args.v_input_file);
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

    const size_t canvas_width  = 2 * graph_width + 2 * COLORBAR_WIDTH + 8 * GRAPH_PADDING;
    const size_t canvas_height = graph_height + TEXT_HEIGHT;
    Rd::Canvas canvas(canvas_width + canvas_width % 2, canvas_height + canvas_height % 2);
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
    Igor::on_death.emplace_back([child = ffmpeg.child()]() {
      if (kill(child, SIGKILL) == -1) {
        Igor::Warn("Could not kill FFmpeg: {}", std::strerror(errno));
      }
    });

    Eigen::Vector<Float, DIM> u_min{};
    Eigen::Vector<Float, DIM> u_max{};
    Eigen::Vector<Float, DIM> v_min{};
    Eigen::Vector<Float, DIM> v_max{};
    for (size_t iter = 0; true; ++iter) {
      const auto got_next_t = t_reader.read_next<false>();
      const auto got_next_u =
          std::visit([](auto& r) { return r.template read_next<false>(); }, u_reader);
      const auto got_next_v =
          std::visit([](auto& r) { return r.template read_next<false>(); }, v_reader);
      if (!got_next_u && !got_next_v && !got_next_t) {
        break;
      } else if (!got_next_u) {
        Igor::Warn("No more u-entries.");
        return false;
      } else if (!got_next_v) {
        Igor::Warn("No more v-entries.");
        return false;
      } else if (!got_next_t) {
        Igor::Warn("No more t-entries.");
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
        if (const std::string str = std::format("t = {:.6f}, ||u|| = {:.6f}, ||v|| = {:.6f}",
                                                t_reader(0, 0),
                                                calculate_mass<Float, DIM>(u_reader)(0),
                                                calculate_mass<Float, DIM>(v_reader)(0));
            !canvas.draw_text(str, text_box, true)) {
          Igor::Warn("Could not render string `{}` to canvas.", str);
          return false;
        }
      }

      if (std::holds_alternative<CellReader>(u_reader)) {
        const auto& cells = std::get<CellReader>(u_reader).cells();
        assert(cells.size() > 0);
        if (iter == 0 || !args.keep_range) {
          std::tie(u_min, u_max) = min_max_cell_value<Float, DIM>(cells);
        }
      } else {
        const auto& data = std::get<MatrixReader>(u_reader).data();
        if (iter == 0 || !args.keep_range) {
          const auto min_max = std::minmax_element(std::cbegin(data), std::cend(data));
          u_min(0)           = *min_max.first;
          u_max(0)           = *min_max.second;
        }
      }
      if (std::holds_alternative<CellReader>(v_reader)) {
        const auto& cells = std::get<CellReader>(v_reader).cells();
        assert(cells.size() > 0);
        if (iter == 0 || !args.keep_range) {
          std::tie(v_min, v_max) = min_max_cell_value<Float, DIM>(cells);
        }
      } else {
        const auto& data = std::get<MatrixReader>(v_reader).data();
        if (iter == 0 || !args.keep_range) {
          const auto min_max = std::minmax_element(std::cbegin(data), std::cend(data));
          v_min(0)           = *min_max.first;
          v_max(0)           = *min_max.second;
        }
      }

      std::visit(
          [&](auto& reader) {
            bool _ = Rd::render_graph<Float, DIM>(
                canvas, reader, 0, u_min, u_max, left_graph_box, args.scale);
          },
          u_reader);

      std::visit(
          [&](auto& reader) {
            bool _ = Rd::render_graph<Float, DIM>(
                canvas, reader, 0, v_min, v_max, right_graph_box, args.scale);
          },
          v_reader);

      constexpr size_t NUM_COLORBAR_BLOCKS = 10;
      if (!Rd::render_color_bar<Float, DIM>(
              canvas, left_colorbar_box, 0, u_min, u_max, NUM_COLORBAR_BLOCKS)) {
        // return false;
      }
      if (!Rd::render_color_bar<Float, DIM>(
              canvas, right_colorbar_box, 0, v_min, v_max, NUM_COLORBAR_BLOCKS)) {
        // return false;
      }

      switch (args.save_frame_images) {
        case Args::PPM:
          if (std::string filename = std::format("frame_{:0>4}.ppm", iter);
              !canvas.to_ppm(filename)) {
            Igor::Warn("Could not write canvas to file `{}`: {}.", filename, std::strerror(errno));
          }
          break;
        case Args::JPEG:
          if (std::string filename = std::format("frame_{:0>4}.jpeg", iter);
              !canvas.to_jpeg(filename)) {
            Igor::Warn("Could not write canvas to file `{}`: {}.", filename, std::strerror(errno));
          }
          break;
        case Args::NO_SAVE: break;
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

  // std::string v_input_file;
  // bool symmetry_error = false;
  // ImageFormat save_frame_images = NO_SAVE;

  Igor::Info("scale        = {}", args->scale);
  Igor::Info("fps          = {}", args->fps);
  Igor::Info("no-mass      = {}", args->no_mass);
  Igor::Info("keep-range   = {}", args->keep_range);
  Igor::Info("v input file = {}", args->v_input_file.empty() ? "None" : args->v_input_file);
  Igor::Info("u input file = {}", args->u_input_file);
  Igor::Info("t input file = {}", args->t_input_file);
  Igor::Info("output file  = {}", args->output_file);

  if (args->v_input_file.empty()) {
    return render_u(*args) ? 0 : 1;
  } else {
    return render_uv(*args) ? 0 : 1;
  }
}
