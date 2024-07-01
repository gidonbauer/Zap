#ifndef ZAP_RENDERER_CANVAS_HPP_
#define ZAP_RENDERER_CANVAS_HPP_

#include <cassert>
#include <fstream>
#include <vector>

#include <unistd.h>

#include <ft2build.h>
#include FT_FREETYPE_H
#include <freetype/ftglyph.h>

#include "Color.hpp"
#include "Igor.hpp"

namespace Zap::Renderer {

struct Box {
  size_t col{};
  size_t row{};
  size_t width{};
  size_t height{};

  [[nodiscard]] constexpr auto contains(size_t x, size_t y) const noexcept -> bool {
    return (x >= col) && (x < col + width) && (y >= row) && (y < row + height);
  }
};

}  // namespace Zap::Renderer

template <>
struct std::formatter<Zap::Renderer::Box, char> {
  template <typename ParseContext>
  static constexpr auto parse(ParseContext& ctx) noexcept {
    return ctx.begin();
  }
  template <typename FormatContext>
  static constexpr auto format(const Zap::Renderer::Box& box, FormatContext& ctx) noexcept {
    return std::format_to(ctx.out(),
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
  [[nodiscard]] static constexpr auto get_gylph_box(Box text_box,
                                                    FT_GlyphSlot glyph,
                                                    size_t& hori_pos,
                                                    Box& glyph_box) noexcept -> bool {
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
        .col    = hori_pos + g_bear_x,
        .row    = text_box.row + above_line_height - g_bear_y,
        .width  = g_width,
        .height = g_height,
    };
    hori_pos += g_advance;
    return true;
  }

 public:
  // -----------------------------------------------------------------------------------------------
  constexpr void clear(RGB color) noexcept {
    std::fill(std::begin(m_data), std::end(m_data), color);
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto
  draw_text(std::string_view str, Box box, bool centered) noexcept -> bool {
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
      if (!load_letter(str[idx])) {
        return false;
      }
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
      Igor::Warn("Text is too long ({}) for box with width {}", length, box.width);
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
      if (is_y_upwards) {
        row = box.height - 1UZ - row;
      }

      if (is_row_major) {
        return col + row * box.width;
      } else {
        return row + col * box.height;
      }
    };

    for (size_t row = 0; row < box.height; ++row) {
      for (size_t col = 0; col < box.width; ++col) {
        const size_t c_col = box.col + col;
        assert(c_col < m_width);
        const size_t c_row = box.row + row;
        assert(c_row < m_height);

        m_data[c_col + c_row * m_width] = to_pixel(buffer[buffer_idx(row, col)]);  // NOLINT
      }
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  constexpr auto draw_rect(Box rect,
                           Box bounding_box,
                           PixelType color,
                           bool is_y_upwards = false) noexcept -> bool {
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
        assert(bounding_box.contains(c_col, c_row));
        m_data[c_col + c_row * m_width] = color;  // NOLINT
      }
    }

    return true;
  }

  // constexpr auto draw_rect(size_t x,
  //                          size_t y,
  //                          size_t w,
  //                          size_t h,
  //                          PixelType color,
  //                          bool is_y_upwards = false) noexcept -> bool {
  //   if (x >= m_width || x + w > m_width) {
  //     Igor::Warn("Rect outside of canvas: canvas width = {}, x = {}, x+w = {}", m_width, x, x +
  //     w); return false;
  //   }
  //   if (y >= m_height || y + h > m_height) {
  //     Igor::Warn(
  //         "Rect outside of canvas: canvas height = {}, y = {}, y+h = {}", m_height, y, y + h);
  //     return false;
  //   }

  //   for (size_t row = y; row < y + h; ++row) {
  //     for (size_t col = x; col < x + w; ++col) {
  //       if (is_y_upwards) {
  //         m_data[col + (m_height - 1UZ - row) * m_width] = color;
  //       } else {
  //         m_data[col + row * m_width] = color;
  //       }
  //     }
  //   }

  //   return true;
  // }

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
    if (!filename.ends_with(".ppm")) {
      filename += ".ppm";
    }

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
};

}  // namespace Zap::Renderer

#endif  // ZAP_RENDERER_CANVAS_HPP_
