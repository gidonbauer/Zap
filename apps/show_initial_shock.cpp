#include "IO/IncCellReader.hpp"
#include "IO/IncMatrixReader.hpp"
#include "Igor.hpp"
#include "Renderer/Canvas.hpp"

namespace Rd = Zap::Renderer;

[[maybe_unused]] constexpr Rd::RGB BLACK  = {.r = 0x00, .g = 0x00, .b = 0x00};
[[maybe_unused]] constexpr Rd::RGB RED    = {.r = 0xFF, .g = 0x00, .b = 0x00};
[[maybe_unused]] constexpr Rd::RGB WHITE  = {.r = 0xFF, .g = 0xFF, .b = 0xFF};
[[maybe_unused]] constexpr Rd::RGB PINK   = {.r = 0xFF, .g = 0x00, .b = 0xFF};
[[maybe_unused]] constexpr Rd::RGB ORANGE = {.r = 0xFF, .g = 0xA5, .b = 0x00};

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
auto main() -> int {
  using Float            = double;
  constexpr size_t DIM   = 1;
  constexpr size_t scale = 100;

  try {
    // - Setup canvas ------------------------------------------------------------------------------
    constexpr auto t_input_file = "../output/cell_based/t.mat";
    constexpr auto u_input_file = "../output/cell_based/u.grid";
    constexpr auto output_file  = "test.ppm";

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
    IGOR_DEBUG_PRINT(graph_box);
    constexpr std::string_view font_file = "../assets/LinLibertine_RI.ttf";
    if (!canvas.load_font(font_file)) {
      Igor::Warn("Could not load font from file `{}`.", font_file);
      return 1;
    }
    // - Setup canvas ------------------------------------------------------------------------------

    if (!t_reader.read_next()) {
      return 1;
    }
    if (!u_reader.read_next()) {
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
      assert(cell.is_cartesian() || cell.is_cut());
      if (cell.is_cartesian()) {
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
          min(i) = std::min(min(i), cell.get_cartesian().value(i));
          max(i) = std::max(max(i), cell.get_cartesian().value(i));
        }
      } else {
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
          min(i) = std::min(min(i),
                            std::min(cell.get_cut().left_value(i), cell.get_cut().right_value(i)));
          max(i) = std::max(max(i),
                            std::max(cell.get_cut().left_value(i), cell.get_cut().right_value(i)));
        }
      }
    }

    bool failed_drawing_rects = false;
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
        const auto c = Rd::float_to_rgb(cell.get_cartesian().value[0], min[0], max[0]);

        if (!canvas.draw_rect(rect, graph_box, c, true)) {
          Igor::Warn("Could not draw cell");
          failed_drawing_rects = true;
        }
      } else {
        const auto& cut_value = cell.get_cut();

        std::vector<Eigen::Vector<Float, 2>> left_points =
            Zap::CellBased::get_left_points<decltype(cell), Float>(cell);
        std::vector<Rd::Point> left_canvas_points(left_points.size());
        for (size_t i = 0; i < left_points.size(); ++i) {
          left_canvas_points[i] = Rd::Point{
              .col =
                  to_pixel_coord(norm_point(left_points[i](0), u_reader.x_min(), u_reader.x_max()),
                                 u_reader.nx(),
                                 scale),
              .row =
                  to_pixel_coord(norm_point(left_points[i](1), u_reader.y_min(), u_reader.y_max()),
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

        std::vector<Eigen::Vector<Float, 2>> right_points =
            Zap::CellBased::get_right_points<decltype(cell), Float>(cell);
        std::vector<Rd::Point> right_canvas_points(right_points.size());
        for (size_t i = 0; i < right_points.size(); ++i) {
          right_canvas_points[i] = Rd::Point{
              .col =
                  to_pixel_coord(norm_point(right_points[i](0), u_reader.x_min(), u_reader.x_max()),
                                 u_reader.nx(),
                                 scale),
              .row =
                  to_pixel_coord(norm_point(right_points[i](1), u_reader.y_min(), u_reader.y_max()),
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
            .col = to_pixel_coord(
                norm_point(cell.get_cut().x1_cut, u_reader.x_min(), u_reader.x_max()),
                u_reader.nx(),
                scale),
            .row = to_pixel_coord(
                norm_point(cell.get_cut().y1_cut, u_reader.y_min(), u_reader.y_max()),
                u_reader.ny(),
                scale),
        };
        const Rd::Point p2{
            .col = to_pixel_coord(
                norm_point(cell.get_cut().x2_cut, u_reader.x_min(), u_reader.x_max()),
                u_reader.nx(),
                scale),
            .row = to_pixel_coord(
                norm_point(cell.get_cut().y2_cut, u_reader.y_min(), u_reader.y_max()),
                u_reader.ny(),
                scale),
        };

        if (!canvas.draw_line(p1, p2, graph_box, RED, true)) {
          Igor::Warn("Could not draw line");
          failed_drawing_rects = true;
        }
      }
    }
    if (failed_drawing_rects) {
      return 1;
    }

    if (!canvas.to_ppm(output_file)) {
      Igor::Warn("Could not write canvas to file.");
      return 1;
    }
  } catch (const std::exception& e) {
    std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC) << e.what() << '\n';
    return 1;
  }
}
