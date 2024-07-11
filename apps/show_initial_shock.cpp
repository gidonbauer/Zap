#include "IO/IncCellReader.hpp"
#include "IO/IncMatrixReader.hpp"
#include "Igor.hpp"
#include "Renderer/Canvas.hpp"

namespace Rd = Zap::Renderer;

constexpr Rd::RGB BLACK  = {.r = 0x00, .g = 0x00, .b = 0x00};
constexpr Rd::RGB RED    = {.r = 0xFF, .g = 0x00, .b = 0x00};
constexpr Rd::RGB WHITE  = {.r = 0xFF, .g = 0xFF, .b = 0xFF};
constexpr Rd::RGB PINK   = {.r = 0xFF, .g = 0x00, .b = 0xFF};
constexpr Rd::RGB ORANGE = {.r = 0xFF, .g = 0xA5, .b = 0x00};

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
[[nodiscard]] constexpr auto
draw_cut_top_right(Rd::Canvas& canvas,
                   const Zap::IO::IncCellReader<Float, DIM>& u_reader,
                   const typename Zap::IO::IncCellReader<Float, DIM>::ReducedCell& cell,
                   const Rd::Box& graph_box,
                   Eigen::Vector<Float, DIM> min,
                   Eigen::Vector<Float, DIM> max,
                   size_t scale) noexcept -> bool {
  const auto& cut_value = cell.get_cut();
  // Left side
  {
    // const auto c = Rd::float_to_rgb(cut_value.right_value(0), min(0), max(0));
    const auto c = ORANGE;

    const Rd::Box rect1{
        .col = to_pixel_coord(
            norm_point(cell.x_min, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cell.y_min, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
        .width = to_pixel_coord(
            norm_length(cell.dx, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .height = to_pixel_coord(norm_length(cell.dy - (cut_value.y1_cut - cell.y_min),
                                             u_reader.y_min(),
                                             u_reader.y_max()),
                                 u_reader.ny(),
                                 scale),
    };

    if (!canvas.draw_rect(rect1, graph_box, c, true)) {
      return false;
    }

    const Rd::Box rect2{
        .col = to_pixel_coord(
            norm_point(cell.x_min, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cell.y_min, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
        .width  = to_pixel_coord(norm_length(cell.dx - (cut_value.x2_cut - cell.x_min),
                                            u_reader.x_min(),
                                            u_reader.x_max()),
                                u_reader.nx(),
                                scale),
        .height = to_pixel_coord(
            norm_length(cell.dy, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };

    if (!canvas.draw_rect(rect2, graph_box, c, true)) {
      return false;
    }

    const Rd::Point p1{
        .col = to_pixel_coord(
            norm_point(cut_value.x1_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y1_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p2{
        .col = to_pixel_coord(
            norm_point(cut_value.x2_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y2_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p3{
        .col = to_pixel_coord(norm_point(std::min(cut_value.x1_cut, cut_value.x2_cut),
                                         u_reader.x_min(),
                                         u_reader.x_max()),
                              u_reader.nx(),
                              scale),
        .row = to_pixel_coord(norm_point(std::min(cut_value.y1_cut, cut_value.y2_cut),
                                         u_reader.y_min(),
                                         u_reader.y_max()),
                              u_reader.ny(),
                              scale),
    };

    if (!canvas.draw_triangle(p1, p2, p3, graph_box, c, true)) {
      return false;
    }
  }

  // Right side
  {
    // const auto c = Rd::float_to_rgb(cut_value.left_value(0), min(0), max(0));
    const auto c = PINK;

    const Rd::Point p1{
        .col = to_pixel_coord(
            norm_point(cut_value.x1_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y1_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p2{
        .col = to_pixel_coord(
            norm_point(cut_value.x2_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y2_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p3{
        .col = to_pixel_coord(norm_point(cell.x_min + cell.dx, u_reader.x_min(), u_reader.x_max()),
                              u_reader.nx(),
                              scale),
        .row = to_pixel_coord(norm_point(cell.y_min + cell.dy, u_reader.y_min(), u_reader.y_max()),
                              u_reader.ny(),
                              scale),
    };

    if (!canvas.draw_triangle(p1, p2, p3, graph_box, c, true)) {
      return false;
    }
  }

  return true;
}

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto
draw_cut_bottom_left(Rd::Canvas& canvas,
                     const Zap::IO::IncCellReader<Float, DIM>& u_reader,
                     const typename Zap::IO::IncCellReader<Float, DIM>::ReducedCell& cell,
                     const Rd::Box& graph_box,
                     Eigen::Vector<Float, DIM> min,
                     Eigen::Vector<Float, DIM> max,
                     size_t scale) noexcept -> bool {
  const auto& cut_value = cell.get_cut();
  // Left side
  {
    // const auto c = Rd::float_to_rgb(cut_value.left_value(0), min(0), max(0));
    const auto c = ORANGE;

    const Rd::Point p1{
        .col = to_pixel_coord(
            norm_point(cut_value.x1_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y1_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p2{
        .col = to_pixel_coord(
            norm_point(cut_value.x2_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y2_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p3{
        .col = to_pixel_coord(
            norm_point(cell.x_min, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cell.y_min, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };

    if (!canvas.draw_triangle(p1, p2, p3, graph_box, c, true)) {
      return false;
    }
  }

  // Right side
  {
    // const auto c = Rd::float_to_rgb(cut_value.right_value(0), min(0), max(0));
    const auto c = PINK;

    const Rd::Box rect1{
        .col = to_pixel_coord(
            norm_point(cut_value.x1_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y1_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
        .width  = to_pixel_coord(norm_length(cell.dx - (cut_value.x1_cut - cell.x_min),
                                            u_reader.x_min(),
                                            u_reader.x_max()),
                                u_reader.nx(),
                                scale),
        .height = to_pixel_coord(
            norm_length(cell.dy, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };

    if (!canvas.draw_rect(rect1, graph_box, c, true)) {
      return false;
    }

    const Rd::Box rect2{
        .col = to_pixel_coord(
            norm_point(cut_value.x2_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y2_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
        .width = to_pixel_coord(
            norm_length(cell.dx, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .height = to_pixel_coord(norm_length(cell.dy - (cut_value.y1_cut - cell.y_min),
                                             u_reader.y_min(),
                                             u_reader.y_max()),
                                 u_reader.ny(),
                                 scale),
    };

    if (!canvas.draw_rect(rect2, graph_box, c, true)) {
      return false;
    }

    const Rd::Point p1{
        .col = to_pixel_coord(
            norm_point(cut_value.x1_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y1_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p2{
        .col = to_pixel_coord(
            norm_point(cut_value.x2_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y2_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p3{
        .col = to_pixel_coord(norm_point(std::max(cut_value.x1_cut, cut_value.x2_cut),
                                         u_reader.x_min(),
                                         u_reader.x_max()),
                              u_reader.nx(),
                              scale),
        .row = to_pixel_coord(norm_point(std::max(cut_value.y1_cut, cut_value.y2_cut),
                                         u_reader.y_min(),
                                         u_reader.y_max()),
                              u_reader.ny(),
                              scale),
    };

    if (!canvas.draw_triangle(p1, p2, p3, graph_box, c, true)) {
      return false;
    }
  }

  return true;
}

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto
draw_cut_middle_hori(Rd::Canvas& canvas,
                     const Zap::IO::IncCellReader<Float, DIM>& u_reader,
                     const typename Zap::IO::IncCellReader<Float, DIM>::ReducedCell& cell,
                     const Rd::Box& graph_box,
                     Eigen::Vector<Float, DIM> min,
                     Eigen::Vector<Float, DIM> max,
                     size_t scale) noexcept -> bool {
  const auto& cut_value = cell.get_cut();
  // Left side
  {
    // const auto c = Rd::float_to_rgb(cut_value.left_value(0), min(0), max(0));
    const auto c = ORANGE;

    const Rd::Box rect{
        .col = to_pixel_coord(
            norm_point(cell.x_min, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cell.y_min, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
        .width = to_pixel_coord(
            norm_length(cell.dx, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .height =
            to_pixel_coord(norm_length(std::min(cut_value.y1_cut, cut_value.y2_cut) - cell.y_min,
                                       u_reader.y_min(),
                                       u_reader.y_max()),
                           u_reader.ny(),
                           scale),
    };

    if (!canvas.draw_rect(rect, graph_box, c, true)) {
      return false;
    }

    const Rd::Point p1{
        .col = to_pixel_coord(
            norm_point(cut_value.x1_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y1_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p2{
        .col = to_pixel_coord(
            norm_point(cut_value.x2_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y2_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p3{
        .col = to_pixel_coord(
            norm_point(cut_value.y1_cut < cut_value.y2_cut ? cut_value.x2_cut : cut_value.x1_cut,
                       u_reader.x_min(),
                       u_reader.x_max()),
            u_reader.nx(),
            scale),
        .row = to_pixel_coord(norm_point(std::min(cut_value.y1_cut, cut_value.y2_cut),
                                         u_reader.y_min(),
                                         u_reader.y_max()),
                              u_reader.ny(),
                              scale),
    };

    if (!canvas.draw_triangle(p1, p2, p3, graph_box, c, true)) {
      return false;
    }
  }

  // Right side
  {
    // const auto c = Rd::float_to_rgb(cut_value.right_value(0), min(0), max(0));
    const auto c = PINK;

    const Rd::Box rect{
        .col = to_pixel_coord(
            norm_point(cell.x_min, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row   = to_pixel_coord(norm_point(std::max(cut_value.y1_cut, cut_value.y2_cut),
                                         u_reader.y_min(),
                                         u_reader.y_max()),
                              u_reader.ny(),
                              scale),
        .width = to_pixel_coord(
            norm_length(cell.dx, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .height = to_pixel_coord(
            norm_length(cell.dy - (std::max(cut_value.y1_cut, cut_value.y2_cut) - cell.y_min),
                        u_reader.y_min(),
                        u_reader.y_max()),
            u_reader.ny(),
            scale),
    };

    if (!canvas.draw_rect(rect, graph_box, c, true)) {
      return false;
    }

    const Rd::Point p1{
        .col = to_pixel_coord(
            norm_point(cut_value.x1_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y1_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p2{
        .col = to_pixel_coord(
            norm_point(cut_value.x2_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y2_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p3{
        .col = to_pixel_coord(
            norm_point(cut_value.y1_cut > cut_value.y2_cut ? cut_value.x2_cut : cut_value.x1_cut,
                       u_reader.x_min(),
                       u_reader.x_max()),
            u_reader.nx(),
            scale),
        .row = to_pixel_coord(norm_point(std::max(cut_value.y1_cut, cut_value.y2_cut),
                                         u_reader.y_min(),
                                         u_reader.y_max()),
                              u_reader.ny(),
                              scale),
    };

    if (!canvas.draw_triangle(p1, p2, p3, graph_box, c, true)) {
      return false;
    }
  }

  return true;
}

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
[[nodiscard]] constexpr auto
draw_cut_middle_vert(Rd::Canvas& canvas,
                     const Zap::IO::IncCellReader<Float, DIM>& u_reader,
                     const typename Zap::IO::IncCellReader<Float, DIM>::ReducedCell& cell,
                     const Rd::Box& graph_box,
                     Eigen::Vector<Float, DIM> min,
                     Eigen::Vector<Float, DIM> max,
                     size_t scale) noexcept -> bool {
  const auto& cut_value = cell.get_cut();
  // Left side
  {
    // const auto c = Rd::float_to_rgb(cut_value.left_value(0), min(0), max(0));
    const auto c = ORANGE;

    const Rd::Box rect{
        .col = to_pixel_coord(
            norm_point(cell.x_min, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cell.y_min, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
        .width =
            to_pixel_coord(norm_length(std::min(cut_value.x1_cut, cut_value.x2_cut) - cell.x_min,
                                       u_reader.x_min(),
                                       u_reader.x_max()),
                           u_reader.nx(),
                           scale),
        .height = to_pixel_coord(
            norm_length(cell.dy, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };

    if (!canvas.draw_rect(rect, graph_box, c, true)) {
      return false;
    }

    const Rd::Point p1{
        .col = to_pixel_coord(
            norm_point(cut_value.x1_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y1_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p2{
        .col = to_pixel_coord(
            norm_point(cut_value.x2_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y2_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p3{
        .col = to_pixel_coord(norm_point(std::min(cut_value.x1_cut, cut_value.x2_cut),
                                         u_reader.x_min(),
                                         u_reader.x_max()),
                              u_reader.nx(),
                              scale),
        .row = to_pixel_coord(
            norm_point(cut_value.x1_cut < cut_value.x2_cut ? cut_value.y2_cut : cut_value.y1_cut,
                       u_reader.y_min(),
                       u_reader.y_max()),
            u_reader.ny(),
            scale),
    };

    if (!canvas.draw_triangle(p1, p2, p3, graph_box, c, true)) {
      return false;
    }
  }

  // Right side
  {
    // const auto c = Rd::float_to_rgb(cut_value.right_value(0), min(0), max(0));
    const auto c = PINK;

    const Rd::Box rect{
        .col = to_pixel_coord(norm_point(std::max(cut_value.x1_cut, cut_value.x2_cut),
                                         u_reader.x_min(),
                                         u_reader.x_max()),
                              u_reader.nx(),
                              scale),
        .row = to_pixel_coord(
            norm_point(cell.y_min, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
        .width = to_pixel_coord(
            norm_length(cell.dx - (std::max(cut_value.x1_cut, cut_value.x2_cut) - cell.x_min),
                        u_reader.x_min(),
                        u_reader.x_max()),
            u_reader.nx(),
            scale),
        .height = to_pixel_coord(
            norm_length(cell.dy, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };

    if (!canvas.draw_rect(rect, graph_box, c, true)) {
      return false;
    }

    const Rd::Point p1{
        .col = to_pixel_coord(
            norm_point(cut_value.x1_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y1_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p2{
        .col = to_pixel_coord(
            norm_point(cut_value.x2_cut, u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
        .row = to_pixel_coord(
            norm_point(cut_value.y2_cut, u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
    };
    const Rd::Point p3{
        .col = to_pixel_coord(norm_point(std::max(cut_value.x1_cut, cut_value.x2_cut),
                                         u_reader.x_min(),
                                         u_reader.x_max()),
                              u_reader.nx(),
                              scale),
        .row = to_pixel_coord(
            norm_point(cut_value.x1_cut > cut_value.x2_cut ? cut_value.y2_cut : cut_value.y1_cut,
                       u_reader.y_min(),
                       u_reader.y_max()),
            u_reader.ny(),
            scale),
    };

    if (!canvas.draw_triangle(p1, p2, p3, graph_box, c, true)) {
      return false;
    }
  }

  return true;
}

// -------------------------------------------------------------------------------------------------
auto main() -> int {
  using Float            = double;
  constexpr size_t DIM   = 1;
  constexpr size_t scale = 100;

  try {
    // - Setup canvas ------------------------------------------------------------------------------
    constexpr auto t_input_file = "../output/cell_based/t_21x21.mat";
    constexpr auto u_input_file = "../output/cell_based/u_21x21.grid";
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
        switch (cut_value.type) {
          case Zap::CellBased::CutType::BOTTOM_LEFT:
            {
              if (!draw_cut_bottom_left<Float, DIM>(
                      canvas, u_reader, cell, graph_box, min, max, scale)) {
                failed_drawing_rects = true;
              }
            }
            break;
          case Zap::CellBased::CutType::BOTTOM_RIGHT:
            {
              const Rd::Box rect{
                  .col = to_pixel_coord(norm_point(cell.x_min, u_reader.x_min(), u_reader.x_max()),
                                        u_reader.nx(),
                                        scale),
                  .row = to_pixel_coord(norm_point(cell.y_min, u_reader.y_min(), u_reader.y_max()),
                                        u_reader.ny(),
                                        scale),
                  .width  = to_pixel_coord(norm_length(cell.dx, u_reader.x_min(), u_reader.x_max()),
                                          u_reader.nx(),
                                          scale),
                  .height = to_pixel_coord(norm_length(cell.dy, u_reader.y_min(), u_reader.y_max()),
                                           u_reader.ny(),
                                           scale),
              };

              if (!canvas.draw_rect(rect, graph_box, WHITE, true)) {
                Igor::Warn("Could not draw cell");
                failed_drawing_rects = true;
              }
              Igor::Warn("Cell type BOTTOM_RIGHT is not implemented yet.");
            }
            break;
          case Zap::CellBased::CutType::TOP_RIGHT:
            {
              if (!draw_cut_top_right<Float, DIM>(
                      canvas, u_reader, cell, graph_box, min, max, scale)) {
                failed_drawing_rects = true;
              }
            }
            break;
          case Zap::CellBased::CutType::TOP_LEFT:
            {
              const Rd::Box rect{
                  .col = to_pixel_coord(norm_point(cell.x_min, u_reader.x_min(), u_reader.x_max()),
                                        u_reader.nx(),
                                        scale),
                  .row = to_pixel_coord(norm_point(cell.y_min, u_reader.y_min(), u_reader.y_max()),
                                        u_reader.ny(),
                                        scale),
                  .width  = to_pixel_coord(norm_length(cell.dx, u_reader.x_min(), u_reader.x_max()),
                                          u_reader.nx(),
                                          scale),
                  .height = to_pixel_coord(norm_length(cell.dy, u_reader.y_min(), u_reader.y_max()),
                                           u_reader.ny(),
                                           scale),
              };

              if (!canvas.draw_rect(rect, graph_box, WHITE, true)) {
                Igor::Warn("Could not draw cell");
                failed_drawing_rects = true;
              }
              Igor::Warn("Cell type TOP_LEFT is not implemented yet.");
            }
            break;
          case Zap::CellBased::CutType::MIDDLE_HORI:
            {
              if (!draw_cut_middle_hori<Float, DIM>(
                      canvas, u_reader, cell, graph_box, min, max, scale)) {
                failed_drawing_rects = true;
              }
            }
            break;
          case Zap::CellBased::CutType::MIDDLE_VERT:
            {
              if (!draw_cut_middle_vert<Float, DIM>(
                      canvas, u_reader, cell, graph_box, min, max, scale)) {
                failed_drawing_rects = true;
              }
            }
            break;
          default:
            {
              Igor::Warn("Unknown cut type with value `{}`", static_cast<int>(cell.get_cut().type));
              failed_drawing_rects = true;
              continue;
            }
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
        // Igor::Todo("Rendering cut cells is not implemented yet.");
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
