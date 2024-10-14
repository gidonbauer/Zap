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

#endif  // COMMON_HPP_
