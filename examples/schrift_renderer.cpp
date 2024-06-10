#include <algorithm>
#include <cassert>
#include <fstream>

#include <libschrift/schrift.h>

#include "Igor.hpp"

// -------------------------------------------------------------------------------------------------
class ImageWrapper {
 public:
  std::vector<char> pixels;
  SFT_Image image{};

  constexpr ImageWrapper(int width, int height) noexcept {
    assert(width > 0);
    assert(height > 0);

    pixels.resize(static_cast<size_t>(width) * static_cast<size_t>(height));
    image = SFT_Image{.pixels = pixels.data(), .width = width, .height = height};
  }
};

// -------------------------------------------------------------------------------------------------
auto get_image(SFT* sft, SFT_UChar c) noexcept -> std::optional<ImageWrapper> {
  assert(sft != nullptr);

  SFT_Glyph glyph{};
  if (sft_lookup(sft, c, &glyph) < 0) {
    Igor::Warn("Could not lookup glyph for letter '{}'.", static_cast<char>(c));
    return std::nullopt;
  }

  SFT_GMetrics metrics{};
  if (sft_gmetrics(sft, glyph, &metrics) < 0) {
    Igor::Warn("Could not get the glyph metrics for letter '{}'", static_cast<char>(c));
    return std::nullopt;
  }
  assert(metrics.minWidth > 0 && metrics.minHeight > 0);

  ImageWrapper image(metrics.minWidth, metrics.minHeight);

  sft_render(sft, glyph, image.image);

  return image;
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto save_image(const ImageWrapper& img, char letter) noexcept -> bool {
  std::string output_file = "letters/$.pgm";
  output_file[8]          = letter;

  std::ofstream out(output_file);
  if (!out) {
    Igor::Warn("Could not open file `{}` for writing: {}", output_file, std::strerror(errno));
    return false;
  }

  out << "P5" << ' ' << img.image.width << ' ' << img.image.height << '\n' << 255 << '\n';
  if (!out.write(img.pixels.data(),
                 static_cast<std::streamsize>(img.pixels.size() * sizeof(img.pixels[0])))) {
    Igor::Warn("Could not write data to file `{}`: {}", output_file, std::strerror(errno));
    return false;
  }

  // Igor::Info("Saved letter '{}' to file `{}`", letter, output_file);
  return true;
}

// -------------------------------------------------------------------------------------------------
auto main() -> int {
  // TODO: Look into https://freetype.org/

  auto font_free = [](SFT* sft) { sft_freefont(sft->font); };
  SFT stack_sft{};
  std::unique_ptr<SFT, decltype(font_free)> sft{&stack_sft};

  // constexpr auto* font_filename = "../assets/LinuxLibertine-RdWo.ttf";
  constexpr auto* font_filename = "../assets/LinLibertine_R.ttf";
  // constexpr auto* font_filename = "../assets/RobotoMonoNerdFont-Regular.ttf";
  sft->font = sft_loadfile(font_filename);
  if (sft->font == nullptr) {
    Igor::Warn("Could not load font from file `{}`", font_filename);
    return 1;
  }

  sft->xScale = 128;
  sft->yScale = 128;
  sft->flags  = SFT_DOWNWARD_Y;

  for (char letter = '0'; letter <= '9'; ++letter) {
    const auto opt_img = get_image(sft.get(), static_cast<SFT_UChar>(letter));
    if (!opt_img.has_value()) {
      Igor::Warn("Could not get image for letter '{}'.", letter);
      return 1;
    }
    Igor::Info("'{}' =>  {}x{}", letter, opt_img->image.width, opt_img->image.height);
    if (!save_image(*opt_img, letter)) {
      return 1;
    }
  }

  for (char letter = 'a'; letter <= 'z'; ++letter) {
    const auto opt_img = get_image(sft.get(), static_cast<SFT_UChar>(letter));
    if (!opt_img.has_value()) {
      Igor::Warn("Could not get image for letter '{}'.", letter);
      return 1;
    }
    Igor::Info("'{}' =>  {}x{}", letter, opt_img->image.width, opt_img->image.height);
    if (!save_image(*opt_img, letter)) {
      return 1;
    }
  }

  // for (char letter = 'A'; letter <= 'Z'; ++letter) {
  //   const auto opt_img = get_image(sft.get(), static_cast<SFT_UChar>(letter));
  //   if (!opt_img.has_value()) {
  //     Igor::Warn("Could not get image for letter '{}'.", letter);
  //     return 1;
  //   }
  //   Igor::Info("'{}' =>  {}x{}", letter, opt_img->image.width, opt_img->image.height);
  //   if (!save_image(*opt_img, letter)) {
  //     return 1;
  //   }
  // }

  std::cout
      << "\n================================================================================\n\n";

  const char left_letter = 'A';
  SFT_Glyph left_glyph{};
  if (sft_lookup(sft.get(), left_letter, &left_glyph) < 0) {
    Igor::Warn("Could not lookup glyph for letter '{}'.", left_letter);
    return 1;
  }

  const char right_letter = 'j';
  SFT_Glyph right_glyph{};
  if (sft_lookup(sft.get(), right_letter, &right_glyph) < 0) {
    Igor::Warn("Could not lookup glyph for letter '{}'.", right_letter);
    return 1;
  }

  SFT_Kerning kerning{};
  if (sft_kerning(sft.get(), left_glyph, right_glyph, &kerning) < 0) {
    Igor::Warn("Could not get kerning for glyph '{}' and '{}'.", left_letter, right_letter);
    return 1;
  }

  Igor::Info("Kerning for '{}' and '{}':", left_letter, right_letter);
  Igor::Info("kerning.xShift = {}", kerning.xShift);
  Igor::Info("kerning.yShift = {}", kerning.yShift);

  SFT_GMetrics left_metrics{};
  if (sft_gmetrics(sft.get(), left_glyph, &left_metrics) < 0) {
    Igor::Warn("Could not get metrics for glyph of letter '{}'", left_letter);
  }

  SFT_GMetrics right_metrics{};
  if (sft_gmetrics(sft.get(), right_glyph, &right_metrics) < 0) {
    Igor::Warn("Could not get metrics for glyph of letter '{}'", right_letter);
  }

  Igor::Info("left_metrics.advanceWidth = {}", left_metrics.advanceWidth);
  Igor::Info("left_metrics.yOffset = {}", left_metrics.yOffset);
  Igor::Info("right_metrics.advanceWidth = {}", right_metrics.advanceWidth);
  Igor::Info("right_metrics.yOffset = {}", right_metrics.yOffset);
}
