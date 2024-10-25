#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <algorithm>
#include <array>
#include <charconv>
#include <optional>
#include <string_view>

#include "Igor/Logging.hpp"

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto parse_size_t(std::string_view str) noexcept -> size_t {
  size_t num;
  const auto res = std::from_chars(str.data(), str.data() + str.size(), num);

  if (res.ec != std::errc{}) {
    Igor::Panic(
        "Could not parse string `{}` to size_t: {}", str, std::make_error_code(res.ec).message());
  }
  return num;
};

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto parse_double(std::string_view str) noexcept -> double {
  constexpr size_t LOCAL_BUFFER_SIZE = 1024;
  std::array<char, LOCAL_BUFFER_SIZE> local_buffer{};
  IGOR_ASSERT(str.size() < LOCAL_BUFFER_SIZE,
              "Local buffer must have enough space for the string and the null termiator, string "
              "with size {} is too large for buffer with capacity {}.",
              str.size(),
              LOCAL_BUFFER_SIZE);
  std::copy(std::cbegin(str), std::cend(str), std::begin(local_buffer));

  char* end        = nullptr;
  const double val = std::strtod(local_buffer.data(), &end);
  if (val == HUGE_VAL) { Igor::Panic("Could not parse string `{}` to double: Out of range.", str); }
  if (local_buffer.data() == end) { Igor::Panic("Could not parse string `{}` to double.", str); }
  return val;
};

// - Get next command line argument ----------------------------------------------------------------
[[nodiscard]] auto next_arg(int& argc, char**& argv) noexcept -> std::optional<std::string_view> {
  if (argc <= 0) { return std::nullopt; }
  std::string_view res = *argv;
  argc -= 1;
  argv += 1;  // NOLINT
  return res;
}

// - Simpson's rule to integrate a function in 1D --------------------------------------------------
template <typename Float, typename FUNC>
[[nodiscard]] constexpr auto simpsons_rule_1d(FUNC f, Float x_min, Float x_max, size_t n) noexcept
    -> Float {
  IGOR_ASSERT(n > 0, "n must be larger than zero but is {}", n);
  if (n % 2 == 1) { n += 1; }
  Float res      = 0;
  const Float dx = (x_max - x_min) / static_cast<Float>(n);
  for (size_t i = 1; i <= n / 2; ++i) {
    res += f(static_cast<Float>(2 * i - 2) * dx + x_min) +
           4 * f(static_cast<Float>(2 * i - 1) * dx + x_min) +
           f(static_cast<Float>(2 * i) * dx + x_min);
  }
  return res * dx / 3;
}

// - Simpson's rule to integrate a function in 2D --------------------------------------------------
template <typename Float, typename FUNC>
[[nodiscard]] auto simpsons_rule_2d(
    FUNC f, Float x_min, Float x_max, Float y_min, Float y_max, size_t nx, size_t ny) noexcept
    -> Float {
  IGOR_ASSERT(nx > 0, "n must be larger than zero but is {}", nx);
  IGOR_ASSERT(ny > 0, "n must be larger than zero but is {}", ny);
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

  Float res = 0;
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
  return res * dx * dy / 9;
}

#endif  // COMMON_HPP_
