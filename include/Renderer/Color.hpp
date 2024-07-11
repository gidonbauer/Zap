#ifndef ZAP_RENDERER_COLOR_HPP_
#define ZAP_RENDERER_COLOR_HPP_

#include <array>
#include <cassert>
#include <cstdint>

#include "Igor.hpp"

namespace Zap::Renderer {

using Greyscale = std::uint8_t;

struct RGB {
  std::uint8_t r = 0x00;
  std::uint8_t g = 0x00;
  std::uint8_t b = 0x00;
};
static_assert(sizeof(RGB) == 24UZ / 8UZ, "Expect struct RGB to not have any padding.");

constexpr auto greyscale_to_RGB(Greyscale c) noexcept -> RGB { return {.r = c, .g = c, .b = c}; }

enum class Float2RGB {
  GREYSCALE,
  COARSE_COLORMAP,
  COLORMAP,
};

template <Float2RGB conv = Float2RGB::COLORMAP, typename Float, bool WARN_ON_NAN = false>
[[nodiscard]] constexpr auto float_to_rgb(Float value, Float min, Float max) noexcept -> RGB {
  auto is_nan_or_inf = [](Float v) { return std::isnan(v) || std::isinf(v); };
  if (is_nan_or_inf(value) || is_nan_or_inf(min) || is_nan_or_inf(max)) {
    if constexpr (WARN_ON_NAN) {
      Igor::Warn("Got NaN or Inf values in provided data: value = {}, min = {}, max = {}",
                 value,
                 min,
                 max);
    }
    return RGB{.r = 0xFF, .g = 0x00, .b = 0x00};
  }

  if constexpr (conv == Float2RGB::GREYSCALE) {
    const auto c = static_cast<std::uint8_t>(
        std::clamp((value - min) / (max - min), Float{0}, Float{1}) * static_cast<Float>(255));
    return {.r = c, .g = c, .b = c};
  }

  const auto get_pixel = [=]([[maybe_unused]] const auto& xs, Float dx, const auto& cols) {
    const auto norm_value = std::clamp((value - min) / (max - min), Float{0}, Float{1});

    const auto i = static_cast<size_t>(std::round(norm_value / dx));
    assert(i <= cols.size());
    return cols[i];
  };

  if constexpr (conv == Float2RGB::COARSE_COLORMAP) {
    constexpr std::array xs = {
        static_cast<Float>(0.00), static_cast<Float>(0.05), static_cast<Float>(0.10),
        static_cast<Float>(0.15), static_cast<Float>(0.20), static_cast<Float>(0.25),
        static_cast<Float>(0.30), static_cast<Float>(0.35), static_cast<Float>(0.40),
        static_cast<Float>(0.45), static_cast<Float>(0.50), static_cast<Float>(0.55),
        static_cast<Float>(0.60), static_cast<Float>(0.65), static_cast<Float>(0.70),
        static_cast<Float>(0.75), static_cast<Float>(0.80), static_cast<Float>(0.85),
        static_cast<Float>(0.90), static_cast<Float>(0.95), static_cast<Float>(1.00),
    };
    constexpr auto dx         = xs[1] - xs[0];  // Assumes equal spacing
    constexpr std::array cols = {
        RGB{.r = 0x44, .g = 0x01, .b = 0x54}, RGB{.r = 0x47, .g = 0x13, .b = 0x65},
        RGB{.r = 0x48, .g = 0x24, .b = 0x75}, RGB{.r = 0x46, .g = 0x34, .b = 0x80},
        RGB{.r = 0x41, .g = 0x44, .b = 0x87}, RGB{.r = 0x3B, .g = 0x52, .b = 0x8B},
        RGB{.r = 0x35, .g = 0x5F, .b = 0x8D}, RGB{.r = 0x2F, .g = 0x6C, .b = 0x8E},
        RGB{.r = 0x2A, .g = 0x78, .b = 0x8E}, RGB{.r = 0x25, .g = 0x84, .b = 0x8E},
        RGB{.r = 0x21, .g = 0x91, .b = 0x8C}, RGB{.r = 0x1E, .g = 0x9C, .b = 0x89},
        RGB{.r = 0x22, .g = 0xA8, .b = 0x84}, RGB{.r = 0x2F, .g = 0xB4, .b = 0x7C},
        RGB{.r = 0x44, .g = 0xBF, .b = 0x70}, RGB{.r = 0x5E, .g = 0xC9, .b = 0x62},
        RGB{.r = 0x7A, .g = 0xD1, .b = 0x51}, RGB{.r = 0x9B, .g = 0xD9, .b = 0x3C},
        RGB{.r = 0xBD, .g = 0xDF, .b = 0x26}, RGB{.r = 0xDF, .g = 0xE3, .b = 0x18},
        RGB{.r = 0xFD, .g = 0xE7, .b = 0x25},
    };
    static_assert(xs.size() == cols.size());
    return get_pixel(xs, dx, cols);
  }

  if constexpr (conv == Float2RGB::COLORMAP) {
    constexpr std::array xs = {
        static_cast<Float>(0.00), static_cast<Float>(0.01), static_cast<Float>(0.02),
        static_cast<Float>(0.03), static_cast<Float>(0.04), static_cast<Float>(0.05),
        static_cast<Float>(0.06), static_cast<Float>(0.07), static_cast<Float>(0.08),
        static_cast<Float>(0.09), static_cast<Float>(0.10), static_cast<Float>(0.11),
        static_cast<Float>(0.12), static_cast<Float>(0.13), static_cast<Float>(0.14),
        static_cast<Float>(0.15), static_cast<Float>(0.16), static_cast<Float>(0.17),
        static_cast<Float>(0.18), static_cast<Float>(0.19), static_cast<Float>(0.20),
        static_cast<Float>(0.21), static_cast<Float>(0.22), static_cast<Float>(0.23),
        static_cast<Float>(0.24), static_cast<Float>(0.25), static_cast<Float>(0.26),
        static_cast<Float>(0.27), static_cast<Float>(0.28), static_cast<Float>(0.29),
        static_cast<Float>(0.30), static_cast<Float>(0.31), static_cast<Float>(0.32),
        static_cast<Float>(0.33), static_cast<Float>(0.34), static_cast<Float>(0.35),
        static_cast<Float>(0.36), static_cast<Float>(0.37), static_cast<Float>(0.38),
        static_cast<Float>(0.39), static_cast<Float>(0.40), static_cast<Float>(0.41),
        static_cast<Float>(0.42), static_cast<Float>(0.43), static_cast<Float>(0.44),
        static_cast<Float>(0.45), static_cast<Float>(0.46), static_cast<Float>(0.47),
        static_cast<Float>(0.48), static_cast<Float>(0.49), static_cast<Float>(0.50),
        static_cast<Float>(0.51), static_cast<Float>(0.52), static_cast<Float>(0.53),
        static_cast<Float>(0.54), static_cast<Float>(0.55), static_cast<Float>(0.56),
        static_cast<Float>(0.57), static_cast<Float>(0.58), static_cast<Float>(0.59),
        static_cast<Float>(0.60), static_cast<Float>(0.61), static_cast<Float>(0.62),
        static_cast<Float>(0.63), static_cast<Float>(0.64), static_cast<Float>(0.65),
        static_cast<Float>(0.66), static_cast<Float>(0.67), static_cast<Float>(0.68),
        static_cast<Float>(0.69), static_cast<Float>(0.70), static_cast<Float>(0.71),
        static_cast<Float>(0.72), static_cast<Float>(0.73), static_cast<Float>(0.74),
        static_cast<Float>(0.75), static_cast<Float>(0.76), static_cast<Float>(0.77),
        static_cast<Float>(0.78), static_cast<Float>(0.79), static_cast<Float>(0.80),
        static_cast<Float>(0.81), static_cast<Float>(0.82), static_cast<Float>(0.83),
        static_cast<Float>(0.84), static_cast<Float>(0.85), static_cast<Float>(0.86),
        static_cast<Float>(0.87), static_cast<Float>(0.88), static_cast<Float>(0.89),
        static_cast<Float>(0.90), static_cast<Float>(0.91), static_cast<Float>(0.92),
        static_cast<Float>(0.93), static_cast<Float>(0.94), static_cast<Float>(0.95),
        static_cast<Float>(0.96), static_cast<Float>(0.97), static_cast<Float>(0.98),
        static_cast<Float>(0.99), static_cast<Float>(1.00),
    };
    constexpr auto dx         = xs[1] - xs[0];  // Assumes equal spacing
    constexpr std::array cols = {
        RGB{.r = 0x44, .g = 0x01, .b = 0x54}, RGB{.r = 0x45, .g = 0x04, .b = 0x57},
        RGB{.r = 0x46, .g = 0x08, .b = 0x5C}, RGB{.r = 0x46, .g = 0x0B, .b = 0x5E},
        RGB{.r = 0x47, .g = 0x10, .b = 0x63}, RGB{.r = 0x47, .g = 0x13, .b = 0x65},
        RGB{.r = 0x48, .g = 0x17, .b = 0x69}, RGB{.r = 0x48, .g = 0x1A, .b = 0x6C},
        RGB{.r = 0x48, .g = 0x1D, .b = 0x6F}, RGB{.r = 0x48, .g = 0x21, .b = 0x73},
        RGB{.r = 0x48, .g = 0x24, .b = 0x75}, RGB{.r = 0x48, .g = 0x28, .b = 0x78},
        RGB{.r = 0x47, .g = 0x2A, .b = 0x7A}, RGB{.r = 0x47, .g = 0x2E, .b = 0x7C},
        RGB{.r = 0x46, .g = 0x30, .b = 0x7E}, RGB{.r = 0x46, .g = 0x34, .b = 0x80},
        RGB{.r = 0x45, .g = 0x37, .b = 0x81}, RGB{.r = 0x44, .g = 0x3A, .b = 0x83},
        RGB{.r = 0x43, .g = 0x3E, .b = 0x85}, RGB{.r = 0x42, .g = 0x40, .b = 0x86},
        RGB{.r = 0x41, .g = 0x44, .b = 0x87}, RGB{.r = 0x40, .g = 0x46, .b = 0x88},
        RGB{.r = 0x3E, .g = 0x49, .b = 0x89}, RGB{.r = 0x3E, .g = 0x4C, .b = 0x8A},
        RGB{.r = 0x3C, .g = 0x4F, .b = 0x8A}, RGB{.r = 0x3B, .g = 0x52, .b = 0x8B},
        RGB{.r = 0x3A, .g = 0x54, .b = 0x8C}, RGB{.r = 0x38, .g = 0x58, .b = 0x8C},
        RGB{.r = 0x37, .g = 0x5A, .b = 0x8C}, RGB{.r = 0x36, .g = 0x5D, .b = 0x8D},
        RGB{.r = 0x35, .g = 0x5F, .b = 0x8D}, RGB{.r = 0x33, .g = 0x62, .b = 0x8D},
        RGB{.r = 0x32, .g = 0x64, .b = 0x8E}, RGB{.r = 0x31, .g = 0x67, .b = 0x8E},
        RGB{.r = 0x30, .g = 0x6A, .b = 0x8E}, RGB{.r = 0x2F, .g = 0x6C, .b = 0x8E},
        RGB{.r = 0x2E, .g = 0x6F, .b = 0x8E}, RGB{.r = 0x2D, .g = 0x71, .b = 0x8E},
        RGB{.r = 0x2C, .g = 0x73, .b = 0x8E}, RGB{.r = 0x2B, .g = 0x75, .b = 0x8E},
        RGB{.r = 0x2A, .g = 0x78, .b = 0x8E}, RGB{.r = 0x29, .g = 0x7A, .b = 0x8E},
        RGB{.r = 0x28, .g = 0x7D, .b = 0x8E}, RGB{.r = 0x27, .g = 0x80, .b = 0x8E},
        RGB{.r = 0x26, .g = 0x82, .b = 0x8E}, RGB{.r = 0x25, .g = 0x84, .b = 0x8E},
        RGB{.r = 0x24, .g = 0x86, .b = 0x8E}, RGB{.r = 0x23, .g = 0x89, .b = 0x8E},
        RGB{.r = 0x22, .g = 0x8B, .b = 0x8D}, RGB{.r = 0x21, .g = 0x8E, .b = 0x8D},
        RGB{.r = 0x21, .g = 0x91, .b = 0x8C}, RGB{.r = 0x20, .g = 0x92, .b = 0x8C},
        RGB{.r = 0x1F, .g = 0x95, .b = 0x8B}, RGB{.r = 0x1F, .g = 0x97, .b = 0x8B},
        RGB{.r = 0x1F, .g = 0x9A, .b = 0x8A}, RGB{.r = 0x1E, .g = 0x9C, .b = 0x89},
        RGB{.r = 0x1F, .g = 0x9F, .b = 0x88}, RGB{.r = 0x1F, .g = 0xA1, .b = 0x88},
        RGB{.r = 0x20, .g = 0xA3, .b = 0x86}, RGB{.r = 0x21, .g = 0xA6, .b = 0x85},
        RGB{.r = 0x22, .g = 0xA8, .b = 0x84}, RGB{.r = 0x25, .g = 0xAB, .b = 0x82},
        RGB{.r = 0x26, .g = 0xAD, .b = 0x81}, RGB{.r = 0x29, .g = 0xAF, .b = 0x7F},
        RGB{.r = 0x2C, .g = 0xB1, .b = 0x7E}, RGB{.r = 0x2F, .g = 0xB4, .b = 0x7C},
        RGB{.r = 0x32, .g = 0xB6, .b = 0x7A}, RGB{.r = 0x37, .g = 0xB8, .b = 0x78},
        RGB{.r = 0x3B, .g = 0xBB, .b = 0x75}, RGB{.r = 0x3F, .g = 0xBC, .b = 0x73},
        RGB{.r = 0x44, .g = 0xBF, .b = 0x70}, RGB{.r = 0x48, .g = 0xC1, .b = 0x6E},
        RGB{.r = 0x4E, .g = 0xC3, .b = 0x6B}, RGB{.r = 0x52, .g = 0xC5, .b = 0x69},
        RGB{.r = 0x58, .g = 0xC7, .b = 0x65}, RGB{.r = 0x5E, .g = 0xC9, .b = 0x62},
        RGB{.r = 0x63, .g = 0xCB, .b = 0x5F}, RGB{.r = 0x69, .g = 0xCD, .b = 0x5B},
        RGB{.r = 0x6E, .g = 0xCE, .b = 0x58}, RGB{.r = 0x75, .g = 0xD0, .b = 0x54},
        RGB{.r = 0x7A, .g = 0xD1, .b = 0x51}, RGB{.r = 0x81, .g = 0xD3, .b = 0x4D},
        RGB{.r = 0x86, .g = 0xD5, .b = 0x49}, RGB{.r = 0x8E, .g = 0xD6, .b = 0x45},
        RGB{.r = 0x95, .g = 0xD8, .b = 0x40}, RGB{.r = 0x9B, .g = 0xD9, .b = 0x3C},
        RGB{.r = 0xA2, .g = 0xDA, .b = 0x37}, RGB{.r = 0xA8, .g = 0xDB, .b = 0x34},
        RGB{.r = 0xB0, .g = 0xDD, .b = 0x2F}, RGB{.r = 0xB5, .g = 0xDE, .b = 0x2B},
        RGB{.r = 0xBD, .g = 0xDF, .b = 0x26}, RGB{.r = 0xC2, .g = 0xDF, .b = 0x23},
        RGB{.r = 0xCA, .g = 0xE1, .b = 0x1F}, RGB{.r = 0xD2, .g = 0xE2, .b = 0x1B},
        RGB{.r = 0xD8, .g = 0xE2, .b = 0x19}, RGB{.r = 0xDF, .g = 0xE3, .b = 0x18},
        RGB{.r = 0xE5, .g = 0xE4, .b = 0x19}, RGB{.r = 0xEC, .g = 0xE5, .b = 0x1B},
        RGB{.r = 0xF1, .g = 0xE5, .b = 0x1D}, RGB{.r = 0xF8, .g = 0xE6, .b = 0x21},
        RGB{.r = 0xFD, .g = 0xE7, .b = 0x25},
    };
    static_assert(xs.size() == cols.size());
    return get_pixel(xs, dx, cols);
  }
}

}  // namespace Zap::Renderer

template <>
struct std::formatter<Zap::Renderer::RGB, char> {
  template <typename ParseContext>
  static constexpr auto parse(ParseContext& ctx) noexcept {
    return ctx.begin();
  }
  template <typename FormatContext>
  static constexpr auto format(const Zap::Renderer::RGB& col, FormatContext& ctx) noexcept {
    return std::format_to(
        ctx.out(), "{{ .r = {:#04x}, .g = {:#04x}, .b = {:#04x} }}", col.r, col.g, col.b);
  }
};

#endif  // ZAP_RENDERER_COLOR_HPP_
