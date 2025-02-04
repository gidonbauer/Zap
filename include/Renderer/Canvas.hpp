#ifndef ZAP_RENDERER_CANVAS_HPP_
#define ZAP_RENDERER_CANVAS_HPP_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>

#include <unistd.h>

#include <ft2build.h>
#include FT_FREETYPE_H
#include <freetype/ftglyph.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <STB/stb_image_write.h>

#include "Color.hpp"
#include "Igor/Logging.hpp"

namespace Zap::Renderer {

struct Point {
  size_t col{};
  size_t row{};
};

struct Box {
  size_t col{};
  size_t row{};
  size_t width{};
  size_t height{};

  [[nodiscard]] constexpr auto contains(const Point& p) const noexcept -> bool {
    return (p.col >= col) && (p.col < col + width) && (p.row >= row) && (p.row < row + height);
  }
};

}  // namespace Zap::Renderer

template <>
struct fmt::formatter<Zap::Renderer::Point, char> {
  template <typename ParseContext>
  static constexpr auto parse(ParseContext& ctx) noexcept {
    return ctx.begin();
  }
  template <typename FormatContext>
  static constexpr auto format(const Zap::Renderer::Point& p, FormatContext& ctx) noexcept {
    return fmt::format_to(ctx.out(), "{{.col = {}, .row = {}}}", p.col, p.row);
  }
};

template <>
struct fmt::formatter<Zap::Renderer::Box, char> {
  template <typename ParseContext>
  static constexpr auto parse(ParseContext& ctx) noexcept {
    return ctx.begin();
  }
  template <typename FormatContext>
  static constexpr auto format(const Zap::Renderer::Box& box, FormatContext& ctx) noexcept {
    return fmt::format_to(ctx.out(),
                          "{{.col = {}, .row = {}, .width = {}, .height = {}}}",
                          box.col,
                          box.row,
                          box.width,
                          box.height);
  }
};

namespace Zap::Renderer {

// -------------------------------------------------------------------------------------------------
class Canvas {
 public:
  using PixelType = RGB;

 private:
  std::vector<PixelType> m_data;
  size_t m_width;
  size_t m_height;

  FT_Library m_ft_library = nullptr;
  FT_Face m_ft_face       = nullptr;

 public:
  constexpr Canvas(size_t width, size_t height) noexcept
      : m_data(width * height),
        m_width(width),
        m_height(height) {}

  constexpr Canvas(const Canvas& other) noexcept                    = delete;
  constexpr Canvas(Canvas&& other) noexcept                         = delete;
  constexpr auto operator=(const Canvas& other) noexcept -> Canvas& = delete;
  constexpr auto operator=(Canvas&& other) noexcept -> Canvas&      = delete;

  constexpr ~Canvas() noexcept {
    if (m_ft_face != nullptr) {
      if (const auto err = FT_Done_Face(m_ft_face); err != FT_Err_Ok) {
        Igor::Warn("Could not clean up face: Exited with error code {}", err);
      }
    }

    if (m_ft_library != nullptr) {
      if (const auto err = FT_Done_FreeType(m_ft_library); err != FT_Err_Ok) {
        Igor::Warn("Could not clean up freetype library: Exited with error code {}", err);
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto width() const noexcept -> size_t { return m_width; }
  [[nodiscard]] constexpr auto height() const noexcept -> size_t { return m_height; }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto load_font(std::string_view font_file,
                                         FT_F26Dot6 font_size = 16) noexcept -> bool {
    if (const auto err = FT_Init_FreeType(&m_ft_library); err != FT_Err_Ok) {
      Igor::Warn("Could not initialize freetype library: Error code {}.", err);
      return false;
    }

    if (const auto err = FT_New_Face(m_ft_library, font_file.data(), 0, &m_ft_face)) {
      Igor::Warn("Could not load font from file `{}`: Error code {}.", font_file, err);
      return false;
    }

    return set_font_size(font_size);
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto set_font_size(FT_F26Dot6 font_size) noexcept -> bool {
    if (const auto err = FT_Set_Char_Size(m_ft_face, 0, font_size * 64L, 400, 400);
        err != FT_Err_Ok) {
      Igor::Warn("Could not set char size: Error code {}.", err);
      return false;
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto box_in_canvas(Box box) const noexcept -> bool {
    return box.col < m_width && box.row < m_height && (box.col + box.width) <= m_width &&
           (box.row + box.height) <= m_height;
  }

 private:
  // -------------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto load_letter(char letter) noexcept -> bool {
    if (m_ft_face == nullptr) {
      Igor::Warn("No font loaded, call `load_font` before drawing text.");
      return false;
    }

    if (const auto err = FT_Load_Char(m_ft_face, static_cast<FT_ULong>(letter), FT_LOAD_DEFAULT);
        err != FT_Err_Ok) {
      Igor::Warn("Could not load glyph for letter '{}'", letter);
      return false;
    }

    if (const auto err = FT_Render_Glyph(m_ft_face->glyph, FT_RENDER_MODE_NORMAL);
        err != FT_Err_Ok) {
      Igor::Warn("Could not render glyph for letter '{}'.", letter);
      return false;
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] static constexpr auto
  get_gylph_box(Box text_box, FT_GlyphSlot glyph, size_t& hori_pos, Box& glyph_box) noexcept
      -> bool {
    const auto above_line_height = 3UZ * text_box.height / 4UZ;

    assert(glyph->metrics.width >= 0);
    const auto g_width = static_cast<size_t>(glyph->metrics.width >> 6);
    assert(glyph->metrics.height >= 0);
    const auto g_height = static_cast<size_t>(glyph->metrics.height >> 6);

    const auto& metrics = glyph->metrics;
    assert(metrics.horiBearingX >= 0);
    const auto g_bear_x = static_cast<size_t>(metrics.horiBearingX >> 6);
    assert(metrics.horiBearingY >= 0);
    const auto g_bear_y = static_cast<size_t>(metrics.horiBearingY >> 6);

    if (g_height > text_box.height) {
      Igor::Warn(
          "Glyph with height {} is too large for line with height {}", g_height, text_box.height);
      return false;
    }
    if (g_bear_y > above_line_height) {
      Igor::Warn("Glyph with y-bearing {} is too large for line with above line height {}",
                 g_bear_y,
                 above_line_height);
      return false;
    }
    // if ((g_height - g_bear_y) > (line_height - above_line_height)) {
    //   Igor::Warn(
    //       "Glyph with height - y-bearing = {} is too large for line with below line height {}",
    //       g_height - g_bear_y,
    //       line_height - above_line_height);
    //   return false;
    // }

    assert(metrics.horiAdvance > 0);
    const auto g_advance = static_cast<size_t>(metrics.horiAdvance >> 6);
    // if (g_advance < g_width) {
    //   IGOR_DEBUG_PRINT(g_advance);
    //   IGOR_DEBUG_PRINT(g_width);
    // }
    // assert(g_advance >= g_width);

    // if (g_advance > m_width - hori_pos) {
    //   Igor::Warn("Glyph is too wide ({}) for canvas with width {} and horizontal position {}",
    //              g_advance,
    //              m_width,
    //              hori_pos);
    //   return false;
    // }

    glyph_box = {
        .col    = text_box.col + hori_pos + g_bear_x,
        .row    = text_box.row + above_line_height - g_bear_y,
        .width  = g_width,
        .height = g_height,
    };
    hori_pos += g_advance;
    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto get_idx(size_t col, size_t row) const noexcept -> size_t {
    assert(col < m_width);
    assert(row < m_height);
    return col + row * m_width;
  }

 public:
  // -----------------------------------------------------------------------------------------------
  constexpr void clear(RGB color) noexcept {
    std::fill(std::begin(m_data), std::end(m_data), color);
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto draw_text(std::string_view str, Box box, bool centered) noexcept
      -> bool {
    if (!box_in_canvas(box)) {
      Igor::Warn("Box {} is not within the canvas with dimension {}x{}", box, m_width, m_height);
      return false;
    }

    struct GlyphWrapper {
      FT_Glyph glyph = nullptr;

      constexpr GlyphWrapper() noexcept                                             = default;
      constexpr GlyphWrapper(const GlyphWrapper& other) noexcept                    = delete;
      constexpr GlyphWrapper(GlyphWrapper&& other) noexcept                         = delete;
      constexpr auto operator=(const GlyphWrapper& other) noexcept -> GlyphWrapper& = delete;
      constexpr auto operator=(GlyphWrapper&& other) noexcept -> GlyphWrapper&      = delete;
      constexpr ~GlyphWrapper() noexcept { FT_Done_Glyph(glyph); }
    };

    size_t length = 0UZ;
    std::vector<Box> glyph_boxes(str.size());
    std::vector<GlyphWrapper> glyphs(str.size());
    for (size_t idx = 0; idx < str.size(); ++idx) {
      if (!load_letter(str[idx])) { return false; }
      const auto& glyph = m_ft_face->glyph;

      if (!get_gylph_box(box, glyph, length, glyph_boxes[idx])) {
        Igor::Warn("Could not get box for glyph.");
        return false;
      }

      if (const auto err = FT_Get_Glyph(glyph, &(glyphs[idx].glyph)); err != FT_Err_Ok) {
        Igor::Warn("Could not copy glyph to local buffer: Error code {}", err);
        return false;
      }
    }

    if (length > box.width) {
      Igor::Warn("Text `{}` is too long ({}) for box with width {}", str, length, box.width);
      return false;
    }

    size_t offset = static_cast<size_t>(centered) * (box.width - length) / 2UZ;
    for (Box& b : glyph_boxes) {
      b.col += offset;
    }

    for (size_t idx = 0; idx < str.size(); ++idx) {
      const auto& glyph     = glyphs[idx].glyph;
      const auto& glyph_box = glyph_boxes[idx];

      if (glyph->format != ft_glyph_format_bitmap) {
        Igor::Warn("Invalid glyph format, needs to be bitmap format.");
        return false;
      }

      FT_BitmapGlyph bm_glyph = reinterpret_cast<FT_BitmapGlyph>(glyph);  // NOLINT
      if (!draw_buffer(bm_glyph->bitmap.buffer, glyph_box, greyscale_to_RGB)) {
        Igor::Warn("Could not draw glyph to canvas.");
        return false;
      }
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  template <typename Buffer, typename PixelConv>
  [[nodiscard]] constexpr auto draw_buffer(const Buffer& buffer,
                                           Box box,
                                           PixelConv to_pixel,
                                           bool is_row_major = true,
                                           bool is_y_upwards = false) noexcept -> bool {
    if (!box_in_canvas(box)) {
      Igor::Warn("Box {} is not within canvas with dimension {}x{}.", box, m_width, m_height);
      return false;
    }

    const auto buffer_idx = [&](size_t row, size_t col) constexpr noexcept -> size_t {
      if (is_y_upwards) { row = box.height - 1UZ - row; }

      if (is_row_major) {
        return col + row * box.width;
      } else {
        return row + col * box.height;
      }
    };

    for (size_t row = 0; row < box.height; ++row) {
      for (size_t col = 0; col < box.width; ++col) {
        const size_t c_col            = box.col + col;
        const size_t c_row            = box.row + row;
        m_data[get_idx(c_col, c_row)] = to_pixel(buffer[buffer_idx(row, col)]);  // NOLINT
      }
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto
  draw_rect(Box rect, Box bounding_box, PixelType color, bool is_y_upwards = false) noexcept
      -> bool {
    if (!box_in_canvas(bounding_box)) {
      Igor::Warn("Bounding box {} is not within canvas with dimension {}x{}.",
                 bounding_box,
                 m_width,
                 m_height);
      return false;
    }

    if ((rect.row + rect.height) > bounding_box.height) {
      Igor::Warn("Rect {} is not in bounding box {}.", rect, bounding_box);
      return false;
    }

    if ((rect.col + rect.width) > bounding_box.width) {
      Igor::Warn("Rect {} is not in bounding box {}.", rect, bounding_box);
      return false;
    }

    for (size_t row = 0; row < rect.height; ++row) {
      for (size_t col = 0; col < rect.width; ++col) {
        const size_t c_col = rect.col + bounding_box.col + col;
        assert(c_col < m_width);
        const size_t c_row = is_y_upwards
                                 ? bounding_box.row + (bounding_box.height - rect.row - row - 1UZ)
                                 : rect.row + bounding_box.row + row;
        assert(c_row < m_height);
        assert(bounding_box.contains({.col = c_col, .row = c_row}));
        m_data[get_idx(c_col, c_row)] = color;  // NOLINT
      }
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto draw_triangle(const Point& p1,
                                             const Point& p2,
                                             const Point& p3,
                                             const Box& bounding_box,
                                             PixelType color,
                                             bool is_y_upwards = false) noexcept -> bool {
    if (!box_in_canvas(bounding_box)) {
      Igor::Warn("Bounding box {} is not in canvas with width {} and height {}.",
                 bounding_box,
                 m_width,
                 m_height);
      return false;
    }
    if (bounding_box.width <= p1.col || bounding_box.height <= p1.row) {
      Igor::Warn("Bounding box {} does not contain point p1 {}.", bounding_box, p1);
      return false;
    }
    if (bounding_box.width <= p2.col || bounding_box.height <= p2.row) {
      Igor::Warn("Bounding box {} does not contain point p2 {}.", bounding_box, p2);
      return false;
    }
    if (bounding_box.width <= p3.col || bounding_box.height <= p3.row) {
      Igor::Warn("Bounding box {} does not contain point p2 {}.", bounding_box, p3);
      return false;
    }

    // From:
    // https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
    const auto normals_pointing_outwards =
        (static_cast<int64_t>(p2.col) - static_cast<int64_t>(p1.col)) *
                (static_cast<int64_t>(p2.row) + static_cast<int64_t>(p1.row)) +
            (static_cast<int64_t>(p3.col) - static_cast<int64_t>(p2.col)) *
                (static_cast<int64_t>(p3.row) + static_cast<int64_t>(p2.row)) +
            (static_cast<int64_t>(p1.col) - static_cast<int64_t>(p3.col)) *
                (static_cast<int64_t>(p1.row) + static_cast<int64_t>(p3.row)) >=
        0;

    const auto min_col = std::min(p1.col, std::min(p2.col, p3.col));
    const auto min_row = std::min(p1.row, std::min(p2.row, p3.row));
    const auto max_col = std::max(p1.col, std::max(p2.col, p3.col));
    const auto max_row = std::max(p1.row, std::max(p2.row, p3.row));

    struct Vec2 {
      int64_t col;
      int64_t row;
    };

    Vec2 n1 = {
        .col = static_cast<int64_t>(p1.row) - static_cast<int64_t>(p2.row),
        .row = -(static_cast<int64_t>(p1.col) - static_cast<int64_t>(p2.col)),
    };
    Vec2 n2 = {
        .col = static_cast<int64_t>(p2.row) - static_cast<int64_t>(p3.row),
        .row = -(static_cast<int64_t>(p2.col) - static_cast<int64_t>(p3.col)),
    };
    Vec2 n3 = {
        .col = static_cast<int64_t>(p3.row) - static_cast<int64_t>(p1.row),
        .row = -(static_cast<int64_t>(p3.col) - static_cast<int64_t>(p1.col)),
    };

    const auto is_in_triangle = [&](size_t row, size_t col) {
      // Line 1: p1 -> p2
      const auto v1 = n1.col * (static_cast<int64_t>(col) - static_cast<int64_t>(p1.col)) +
                      n1.row * (static_cast<int64_t>(row) - static_cast<int64_t>(p1.row));
      // Line 1: p2 -> p3
      const auto v2 = n2.col * (static_cast<int64_t>(col) - static_cast<int64_t>(p2.col)) +
                      n2.row * (static_cast<int64_t>(row) - static_cast<int64_t>(p2.row));
      // Line 1: p3 -> p1
      const auto v3 = n3.col * (static_cast<int64_t>(col) - static_cast<int64_t>(p3.col)) +
                      n3.row * (static_cast<int64_t>(row) - static_cast<int64_t>(p3.row));

      if (normals_pointing_outwards) {
        return v1 <= 0 && v2 <= 0 && v3 <= 0;
      } else {
        return v1 >= 0 && v2 >= 0 && v3 >= 0;
      }
    };

    for (auto row = min_row; row <= max_row; ++row) {
      for (auto col = min_col; col <= max_col; ++col) {
        if (is_in_triangle(row, col)) {
          const size_t c_col = col + bounding_box.col;
          assert(c_col < m_width);
          const size_t c_row = is_y_upwards ? bounding_box.row + (bounding_box.height - row - 1UZ)
                                            : row + bounding_box.row;
          assert(c_row < m_height);
          assert(bounding_box.contains({.col = c_col, .row = c_row}));
          m_data[get_idx(c_col, c_row)] = color;  // NOLINT
        }
      }
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto draw_polygon(std::vector<Point> points,
                                            const Box& bounding_box,
                                            PixelType color,
                                            bool is_y_upwards = false) noexcept -> bool {
    if (points.size() < 3) {
      Igor::Warn("Need at least tree points to draw a polygon, but got only {}", points.size());
      return false;
    }
    if (!box_in_canvas(bounding_box)) {
      Igor::Warn("Bounding box {} is not in canvas with width {} and height {}.",
                 bounding_box,
                 m_width,
                 m_height);
      return false;
    }
    for (size_t i = 0; i < points.size(); ++i) {
      const auto& p = points[i];
      if (bounding_box.width <= p.col || bounding_box.height <= p.row) {
        Igor::Warn("Bounding box {} does not contain point p[{}]={}.", bounding_box, i, p);
        return false;
      }
    }

    // Find centroid and bounding box of polygon
    Point center{0, 0};
    size_t min_col = points.front().col;
    size_t max_col = points.front().col;
    size_t min_row = points.front().row;
    size_t max_row = points.front().row;
    for (const auto& [col, row] : points) {
      center.col += col;
      center.row += row;
      min_col = std::min(min_col, col);
      max_col = std::max(max_col, col);
      min_row = std::min(min_row, row);
      max_row = std::max(max_row, row);
    }
    center.col /= points.size();
    center.row /= points.size();

    // Sort points in clockwise order
    std::sort(std::begin(points), std::end(points), [&center](const Point& p1, const Point& p2) {
      const auto a1 = std::atan2(static_cast<int64_t>(p1.col) - static_cast<int64_t>(center.col),
                                 static_cast<int64_t>(p1.row) - static_cast<int64_t>(center.row));
      const auto a2 = std::atan2(static_cast<int64_t>(p2.col) - static_cast<int64_t>(center.col),
                                 static_cast<int64_t>(p2.row) - static_cast<int64_t>(center.row));
      return a1 < a2;
    });

    // Find outwards pointing normals
    struct Vec2 {
      int64_t col;
      int64_t row;
    };

    std::vector<Vec2> normals(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
      const auto& p1 = points[i];
      const auto& p2 = points[(i + 1) % points.size()];
      auto& n        = normals[i];
      n              = Vec2{
                       .col = static_cast<int64_t>(p1.row) - static_cast<int64_t>(p2.row),
                       .row = -(static_cast<int64_t>(p1.col) - static_cast<int64_t>(p2.col)),
      };
    }

    const auto is_inside = [&](size_t row, size_t col) {
      for (size_t i = 0; i < points.size(); ++i) {
        const auto& n = normals[i];
        const auto& p = points[i];
        const auto v  = n.col * (static_cast<int64_t>(col) - static_cast<int64_t>(p.col)) +
                       n.row * (static_cast<int64_t>(row) - static_cast<int64_t>(p.row));
        if (v > 0) { return false; }
      }
      return true;
    };

    // Draw points inside of polygon
    for (auto row = min_row; row <= max_row; ++row) {
      for (auto col = min_col; col <= max_col; ++col) {
        if (is_inside(row, col)) {
          const size_t c_col = col + bounding_box.col;
          assert(c_col < m_width);
          const size_t c_row = is_y_upwards ? bounding_box.row + (bounding_box.height - row - 1UZ)
                                            : row + bounding_box.row;
          assert(c_row < m_height);
          assert(bounding_box.contains({.col = c_col, .row = c_row}));
          m_data[get_idx(c_col, c_row)] = color;  // NOLINT
        }
      }
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto draw_line(const Point& p1,
                                         const Point& p2,
                                         const Box& bounding_box,
                                         PixelType color,
                                         bool is_y_upwards = false) noexcept -> bool {
    if (!box_in_canvas(bounding_box)) {
      Igor::Warn("Bounding box {} is not in canvas with width {} and height {}.",
                 bounding_box,
                 m_width,
                 m_height);
      return false;
    }
    if (bounding_box.width <= p1.col || bounding_box.height <= p1.row) {
      Igor::Warn("Bounding box {} does not contain point p1 {}.", bounding_box, p1);
      return false;
    }
    if (bounding_box.width <= p2.col || bounding_box.height <= p2.row) {
      Igor::Warn("Bounding box {} does not contain point p2 {}.", bounding_box, p2);
      return false;
    }

    int dcol = static_cast<int>(p1.col > p2.col ? p1.col - p2.col : p2.col - p1.col);
    int drow = -static_cast<int>(p1.row > p2.row ? p1.row - p2.row : p2.row - p1.row);
    int scol = p1.col < p2.col ? 1 : -1;
    int srow = p1.row < p2.row ? 1 : -1;
    int err  = dcol + drow;

    assert(p1.col <= std::numeric_limits<int>::max());
    assert(p1.row <= std::numeric_limits<int>::max());
    int col = static_cast<int>(p1.col);
    int row = static_cast<int>(p1.row);
    while (true) {
      assert(col >= 0);
      assert(row >= 0);
      const size_t c_col = static_cast<size_t>(col) + bounding_box.col;
      const size_t c_row =
          is_y_upwards ? bounding_box.row + (bounding_box.height - static_cast<size_t>(row) - 1UZ)
                       : bounding_box.row + static_cast<size_t>(row);
      m_data[get_idx(c_col, c_row)] = color;
      if (static_cast<size_t>(col) == p2.col && static_cast<size_t>(row) == p2.row) { break; }

      if (static_cast<size_t>(col) == p2.col && static_cast<size_t>(row) == p2.row) { break; }
      int err2 = 2 * err;
      if (err2 > drow) {
        err += drow;
        col += scol;
      }
      if (err2 < dcol) {
        err += dcol;
        row += srow;
      }
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto to_raw_stream(pid_t stream) const noexcept -> bool {
    const auto data_size_bytes = m_data.size() * sizeof(PixelType);
    const auto bytes_writen =
        write(stream, reinterpret_cast<const char*>(m_data.data()), data_size_bytes);  // NOLINT

    if (bytes_writen == -1) {
      Igor::Warn("Writing data to pipe failed: {}", std::strerror(errno));
      return false;
    }
    if (static_cast<size_t>(bytes_writen) < data_size_bytes) {
      Igor::Warn(
          "Wrote only {} bytes to pipe, expected to write {}", bytes_writen, data_size_bytes);
      return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto to_ppm(std::string filename) const noexcept -> bool {
    if (!filename.ends_with(".ppm")) { filename += ".ppm"; }

    std::ofstream out(filename);
    if (!out) {
      Igor::Warn("Could not open file `{}` for writing: {}", filename, std::strerror(errno));
      return false;
    }

    out << "P6" << '\n' << m_width << ' ' << m_height << '\n' << 255 << '\n';
    if (!out.write(reinterpret_cast<const char*>(m_data.data()),  // NOLINT
                   static_cast<std::streamsize>(m_width * m_height * sizeof(PixelType)))) {
      Igor::Warn("Could not write data to file `{}`: {}", filename, std::strerror(errno));
      return false;
    }

    Igor::Info("Saved canvas to `{}`", filename);
    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto to_jpeg(std::string filename, int quality = 100) const noexcept -> bool {
    if (!filename.ends_with(".jpeg")) { filename += ".jpeg"; }

    assert(m_width <= std::numeric_limits<int>::max());
    assert(m_height <= std::numeric_limits<int>::max());
    static_assert(sizeof(PixelType) == 3, "Expect Pixel type to have 3 8-bit components.");
    const auto res = stbi_write_jpg(filename.c_str(),
                                    static_cast<int>(m_width),
                                    static_cast<int>(m_height),
                                    3,
                                    m_data.data(),
                                    quality);

    if (res == 0) {
      Igor::Warn("Could not write data to file `{}`", filename);
      return false;
    }

    Igor::Info("Saved canvas to `{}`", filename);
    return true;
  }
};

}  // namespace Zap::Renderer

#endif  // ZAP_RENDERER_CANVAS_HPP_
