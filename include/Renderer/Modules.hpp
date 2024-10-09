#ifndef ZAP_RENDERER_MODULES_HPP_
#define ZAP_RENDERER_MODULES_HPP_

#include <cstddef>

#include "IO/IncCellReader.hpp"
#include "IO/IncMatrixReader.hpp"
#include "Renderer/Canvas.hpp"

namespace Zap::Renderer {

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
template <typename Float, int DIM>
[[nodiscard]] constexpr auto render_color_bar(Canvas& canvas,
                                              const Box& colorbar_box,
                                              Eigen::Index solution_idx,
                                              const Eigen::Vector<Float, DIM>& min,
                                              const Eigen::Vector<Float, DIM>& max,
                                              bool same_range,
                                              size_t num_colorbar_blocks = 10) noexcept -> bool {
  bool success = true;

  if (!canvas.set_font_size(6)) { success = false; }
  for (size_t i = 0; i < num_colorbar_blocks; ++i) {
    const Float frac = static_cast<Float>(i) / static_cast<Float>(num_colorbar_blocks - 1);
    RGB c{};
    Float value{};
    if (same_range) {
      const auto total_min = std::min_element(std::cbegin(min), std::cend(min));
      const auto total_max = std::max_element(std::cbegin(max), std::cend(max));
      value                = std::lerp(*total_min, *total_max, frac);
      c                    = float_to_rgb(value, *total_min, *total_max);
    } else {
      value = std::lerp(min(solution_idx), max(solution_idx), frac);
      c     = float_to_rgb(value, min(solution_idx), max(solution_idx));
    }

    const Box block = {
        .col    = 0,
        .row    = i * colorbar_box.height / num_colorbar_blocks,
        .width  = colorbar_box.width / 4,
        .height = colorbar_box.height / num_colorbar_blocks,
    };
    if (!canvas.draw_rect(block, colorbar_box, c, true)) { success = false; }

    const auto str_value = std::format("{:.6f}", value);
    const Box text_block = {
        .col = colorbar_box.col + block.width,
        .row = colorbar_box.row +
               (num_colorbar_blocks - i - 1) * colorbar_box.height / num_colorbar_blocks,
        .width  = colorbar_box.width - block.width,
        .height = colorbar_box.height / num_colorbar_blocks,
    };

    if (!canvas.draw_text(str_value, text_block, true)) { success = false; }
  }

  return success;
}

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
[[nodiscard]] auto render_graph(Canvas& canvas,
                                const Zap::IO::IncCellReader<Float, DIM>& u_reader,
                                Eigen::Index solution_idx,
                                const Eigen::Vector<Float, DIM>& min,
                                const Eigen::Vector<Float, DIM>& max,
                                const Box& graph_box,
                                size_t scale,
                                bool same_range) noexcept -> bool {
  assert(scale > 0);
  const auto& cells = u_reader.cells();

  bool success = true;
#ifdef ZAP_PARALLEL_RENDERING
#pragma omp parallel for
#endif  // ZAP_PARALLEL_RENDERING
  for (const auto& cell : cells) {
    if (cell.is_cartesian()) {
      const Box rect{
          .col = to_pixel_coord(
              norm_point(cell.x_min(), u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
          .row = to_pixel_coord(
              norm_point(cell.y_min(), u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
          .width = to_pixel_coord(
              norm_length(cell.dx(), u_reader.x_min(), u_reader.x_max()), u_reader.nx(), scale),
          .height = to_pixel_coord(
              norm_length(cell.dy(), u_reader.y_min(), u_reader.y_max()), u_reader.ny(), scale),
      };
      RGB c{};
      if (same_range) {
        const auto total_min = std::min_element(std::cbegin(min), std::cend(min));
        const auto total_max = std::max_element(std::cbegin(max), std::cend(max));
        c = float_to_rgb(cell.get_cartesian().value(solution_idx), *total_min, *total_max);
      } else {
        c = float_to_rgb(
            cell.get_cartesian().value(solution_idx), min(solution_idx), max(solution_idx));
      }

      if (!canvas.draw_rect(rect, graph_box, c, true)) {
        Igor::Warn("Could not draw cell");
        success = false;
      }
    } else if (cell.is_cut()) {
      const auto& cut_value = cell.get_cut();

      // Left sub-cell
      {
        const auto points = cell.get_left_points();
        std::vector<Point> canvas_points(points.size());
        for (size_t i = 0; i < points.size(); ++i) {
          canvas_points[i] = Point{
              .col = to_pixel_coord(norm_point(points[i].x, u_reader.x_min(), u_reader.x_max()),
                                    u_reader.nx(),
                                    scale),
              .row = to_pixel_coord(norm_point(points[i].y, u_reader.y_min(), u_reader.y_max()),
                                    u_reader.ny(),
                                    scale),
          };
        }

        RGB c{};
        if (same_range) {
          const auto total_min = std::min_element(std::cbegin(min), std::cend(min));
          const auto total_max = std::max_element(std::cbegin(max), std::cend(max));
          c = float_to_rgb(cell.get_cut().left_value(solution_idx), *total_min, *total_max);
        } else {
          c = float_to_rgb(
              cell.get_cut().left_value(solution_idx), min(solution_idx), max(solution_idx));
        }

        if (!canvas.draw_polygon(canvas_points, graph_box, c, true)) {
          Igor::Warn("Could not draw left polygon.");
          success = false;
        }
      }

      // Right sub-cell
      {
        const auto points = cell.get_right_points();
        std::vector<Point> canvas_points(points.size());
        for (size_t i = 0; i < points.size(); ++i) {
          canvas_points[i] = Point{
              .col = to_pixel_coord(norm_point(points[i].x, u_reader.x_min(), u_reader.x_max()),
                                    u_reader.nx(),
                                    scale),
              .row = to_pixel_coord(norm_point(points[i].y, u_reader.y_min(), u_reader.y_max()),
                                    u_reader.ny(),
                                    scale),
          };
        }

        RGB c{};
        if (same_range) {
          const auto total_min = std::min_element(std::cbegin(min), std::cend(min));
          const auto total_max = std::max_element(std::cbegin(max), std::cend(max));
          c = float_to_rgb(cell.get_cut().right_value(solution_idx), *total_min, *total_max);
        } else {
          c = float_to_rgb(
              cell.get_cut().right_value(solution_idx), min(solution_idx), max(solution_idx));
        }

        if (!canvas.draw_polygon(canvas_points, graph_box, c, true)) {
          Igor::Warn("Could not draw right polygon.");
          success = false;
        }
      }

      // Draw shock curve
      const Point p1{
          .col = to_pixel_coord(norm_point(cut_value.cut1.x, u_reader.x_min(), u_reader.x_max()),
                                u_reader.nx(),
                                scale),
          .row = to_pixel_coord(norm_point(cut_value.cut1.y, u_reader.y_min(), u_reader.y_max()),
                                u_reader.ny(),
                                scale),
      };
      const Point p2{
          .col = to_pixel_coord(norm_point(cut_value.cut2.x, u_reader.x_min(), u_reader.x_max()),
                                u_reader.nx(),
                                scale),
          .row = to_pixel_coord(norm_point(cut_value.cut2.y, u_reader.y_min(), u_reader.y_max()),
                                u_reader.ny(),
                                scale),
      };

      if (!canvas.draw_line(p1, p2, graph_box, RGB{.r = 0xFF}, true)) {
        Igor::Warn("Could not draw line");
        success = false;
      }
    } else {
      Igor::Warn("Unknown cell type with variant index {}", cell.value.index());
      success = false;
    }
  }

  return success;
}

// -------------------------------------------------------------------------------------------------
template <typename Float, size_t DIM>
[[nodiscard]] auto render_graph(Canvas& canvas,
                                const Zap::IO::IncMatrixReader<Float>& u_reader,
                                Eigen::Index solution_idx,
                                const Eigen::Vector<Float, DIM>& min,
                                const Eigen::Vector<Float, DIM>& max,
                                const Box& graph_box,
                                size_t scale,
                                bool same_range) noexcept -> bool {
  static_assert(DIM == 1, "2d matrix based solver is not implemented.");
  assert(scale > 0);
  assert(u_reader.rows() > 0);
  assert(u_reader.cols() > 0);
  assert(graph_box.width % static_cast<size_t>(u_reader.rows()) == 0);
  assert(graph_box.height % static_cast<size_t>(u_reader.cols()) == 0);
  assert(solution_idx == 0 && "2d matrix based solver is not implemented.");

  bool success = true;
  for (int64_t row = 0; row < u_reader.rows(); ++row) {
    for (int64_t col = 0; col < u_reader.cols(); ++col) {
      const Box rect{
          .col    = static_cast<size_t>(col) * scale,
          .row    = static_cast<size_t>(row) * scale,
          .width  = scale,
          .height = scale,
      };

      (void)same_range;
      RGB c = float_to_rgb(u_reader(row, col), min(0), max(0));

      if (!canvas.draw_rect(rect, graph_box, c, true)) { success = false; }
    }
  }

  return success;
}

}  // namespace Zap::Renderer

#endif  // ZAP_RENDERER_MODULES_HPP_
