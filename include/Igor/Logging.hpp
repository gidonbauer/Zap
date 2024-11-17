// Copyright 2024 Gidon Bauer <gidon.bauer@rwth-aachen.de>

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:

// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef IGOR_LOGGING_HPP_
#define IGOR_LOGGING_HPP_

#ifndef IGOR_USE_FMT
#include <format>
#else
#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#endif  // IGOR_USE_FMT

#include <cstdint>
#include <iostream>
#include <ranges>
#include <source_location>
#include <string>
#include <utility>

namespace Igor {

enum class ExitCode : int {  // NOLINT(performance-enum-size)
  NO_EARLY_EXIT,
  PANIC,
  TODO,
  ASSERT,
};

namespace detail {

#ifndef IGOR_USE_FMT
using std::format;
using std::format_string;
#else
using fmt::format;
using fmt::format_string;
#endif  // IGOR_USE_FMT

[[nodiscard]] constexpr auto strip_path(std::string_view full_path) noexcept -> std::string_view {
#if defined(WIN32) || defined(_WIN32)
#error "Not implemented yet: Requires different path separator ('\\') and potentially uses wchar"
#else
  constexpr char separator = '/';
#endif

  size_t counter = 0;
  for (char c : std::ranges::reverse_view(full_path)) {
    if (c == separator) { break; }
    ++counter;
  }

  return full_path.substr(full_path.size() - counter, counter);
}

[[nodiscard]] constexpr auto
error_loc(const std::source_location loc = std::source_location::current()) noexcept
    -> std::string {
  try {
    return detail::format("`{}` (\033[95m{}:{}:{}\033[0m)",
                          loc.function_name(),
#ifdef IGOR_ERROR_LOC_FULL_PATH
                          loc.file_name(),
#else
                          strip_path(loc.file_name()),
#endif  // IGOR_ERROR_LOC_FULL_PATH
                          loc.line(),
                          loc.column());
  } catch (const std::exception& e) {
    std::cerr << "Could not format the error location: " << e.what() << '\n';
    std::exit(static_cast<int>(ExitCode::PANIC));
  }
}

enum class Level : std::uint8_t {
  INFO,
  WARN,
  TODO,
  PANIC,
  DEBUG,
  TIME,
  ASSERT,
};

consteval auto level_stream(Level level) noexcept -> std::ostream& {
  switch (level) {
    case Level::INFO:
    case Level::TIME:   return std::cout;
    case Level::WARN:
    case Level::TODO:
    case Level::PANIC:
    case Level::DEBUG:
    case Level::ASSERT: return std::cerr;
  }
}

consteval auto level_repr(Level level) noexcept {
  switch (level) {
    case Level::INFO:   return "\033[32m[INFO]\033[0m ";
    case Level::WARN:   return "\033[33m[WARN]\033[0m ";
    case Level::TODO:   return "\033[35m[TODO]\033[0m ";
    case Level::PANIC:  return "\033[31m[ERROR]\033[0m ";
    case Level::DEBUG:  return "\033[94m[DEBUG]\033[0m ";
    case Level::TIME:   return "\033[94m[TIME]\033[0m ";
    case Level::ASSERT: return "\033[31m[ASSERT]\033[0m ";
  }
}

template <Level level, ExitCode exit_code, typename... Args>
class Print {
 protected:
  constexpr Print(detail::format_string<Args...> fmt, Args&&... args) noexcept {
    auto& out = level_stream(level);
    out << level_repr(level) << detail::format(fmt, std::forward<Args>(args)...) << '\n';
    if constexpr (exit_code != ExitCode::NO_EARLY_EXIT) { std::exit(static_cast<int>(exit_code)); }
  }

  constexpr Print(const std::source_location loc,
                  detail::format_string<Args...> fmt,
                  Args&&... args) noexcept {
    auto& out = level_stream(level);
    out << level_repr(level) << error_loc(loc) << ": "
        << detail::format(fmt, std::forward<Args>(args)...) << '\n';
    if constexpr (exit_code != ExitCode::NO_EARLY_EXIT) { std::exit(static_cast<int>(exit_code)); }
  }
};

// -------------------------------------------------------------------------------------------------
template <typename... Args>
class [[maybe_unused]] Time final
    : detail::Print<detail::Level::TIME, ExitCode::NO_EARLY_EXIT, Args...> {
  using P = detail::Print<detail::Level::TIME, ExitCode::NO_EARLY_EXIT, Args...>;

 public:
  constexpr Time(detail::format_string<Args...> fmt, Args&&... args) noexcept
      : P{fmt, std::forward<Args>(args)...} {}
};
template <typename... Args>
Time(detail::format_string<Args...>, Args&&...) -> Time<Args...>;

}  // namespace detail

// -------------------------------------------------------------------------------------------------
template <typename... Args>
class [[maybe_unused]] Info final
    : detail::Print<detail::Level::INFO, ExitCode::NO_EARLY_EXIT, Args...> {
  using P = detail::Print<detail::Level::INFO, ExitCode::NO_EARLY_EXIT, Args...>;

 public:
  constexpr Info(detail::format_string<Args...> fmt, Args&&... args) noexcept
      : P{fmt, std::forward<Args>(args)...} {}
};
template <typename... Args>
Info(detail::format_string<Args...>, Args&&...) -> Info<Args...>;

// -------------------------------------------------------------------------------------------------
template <typename... Args>
class [[maybe_unused]] Warn final
    : detail::Print<detail::Level::WARN, ExitCode::NO_EARLY_EXIT, Args...> {
  using P = detail::Print<detail::Level::WARN, ExitCode::NO_EARLY_EXIT, Args...>;

 public:
  constexpr Warn(detail::format_string<Args...> fmt,
                 Args&&... args,
                 const std::source_location loc = std::source_location::current()) noexcept
      : P{loc, fmt, std::forward<Args>(args)...} {}
};
template <typename... Args>
Warn(detail::format_string<Args...>, Args&&...) -> Warn<Args...>;

// -------------------------------------------------------------------------------------------------
template <typename... Args>
class [[maybe_unused]] Todo final : detail::Print<detail::Level::TODO, ExitCode::TODO, Args...> {
  using P = detail::Print<detail::Level::TODO, ExitCode::TODO, Args...>;

 public:
  [[noreturn]] constexpr Todo(
      const std::source_location loc = std::source_location::current()) noexcept
      : P{loc, "Not implemented yet."} {
    std::unreachable();
  }
  [[noreturn]] constexpr Todo(
      detail::format_string<Args...> fmt,
      Args&&... args,
      const std::source_location loc = std::source_location::current()) noexcept
      : P{loc, fmt, std::forward<Args>(args)...} {
    std::unreachable();
  }
};
template <typename... Args>
Todo(detail::format_string<Args...>, Args&&...) -> Todo<Args...>;

// -------------------------------------------------------------------------------------------------
template <typename... Args>
class [[maybe_unused]] Panic final : detail::Print<detail::Level::PANIC, ExitCode::PANIC, Args...> {
  using P = detail::Print<detail::Level::PANIC, ExitCode::PANIC, Args...>;

 public:
  [[noreturn]] constexpr Panic(
      detail::format_string<Args...> fmt,
      Args&&... args,
      const std::source_location loc = std::source_location::current()) noexcept
      : P{loc, fmt, std::forward<Args>(args)...} {
    std::unreachable();
  }
};
template <typename... Args>
Panic(detail::format_string<Args...>, Args&&...) -> Panic<Args...>;

// -------------------------------------------------------------------------------------------------
template <typename... Args>
class [[maybe_unused]] Assert final
    : detail::Print<detail::Level::ASSERT, ExitCode::ASSERT, Args...> {
  using P = detail::Print<detail::Level::ASSERT, ExitCode::ASSERT, Args...>;

 public:
  [[noreturn]] constexpr Assert(
      detail::format_string<Args...> fmt,
      Args&&... args,
      const std::source_location loc = std::source_location::current()) noexcept
      : P{loc, fmt, std::forward<Args>(args)...} {
    std::unreachable();
  }
};
template <typename... Args>
Assert(detail::format_string<Args...>, Args&&...) -> Assert<Args...>;

#ifndef IGOR_NDEBUG
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define IGOR_ASSERT(cond, ...)                                                                     \
  do {                                                                                             \
    /* NOLINTNEXTLINE(readability-simplify-boolean-expr) */                                        \
    if (!(cond)) [[unlikely]] {                                                                    \
      using va_args_tuple__ = decltype(std::make_tuple(__VA_ARGS__));                              \
      static_assert(std::tuple_size_v<va_args_tuple__> > 0,                                        \
                    "`IGOR_ASSERT` requires an error message, please provide a format string and " \
                    "the corresponding values");                                                   \
      static_assert(std::is_convertible_v<std::tuple_element_t<0, va_args_tuple__>, std::string>,  \
                    "First argument for error message must be a format string.");                  \
      Igor::Assert("Assertion `{}` failed: {}", #cond, Igor::detail::format(__VA_ARGS__));         \
    }                                                                                              \
  } while (false)
#else
#define IGOR_ASSERT(cond, ...) ((void)0)
#endif  // IGOR_NDEBUG

// -------------------------------------------------------------------------------------------------
template <typename... Args>
class [[maybe_unused]] Debug final
    : detail::Print<detail::Level::DEBUG, ExitCode::NO_EARLY_EXIT, Args...> {
  using P = detail::Print<detail::Level::DEBUG, ExitCode::NO_EARLY_EXIT, Args...>;

 public:
  constexpr Debug(detail::format_string<Args...> fmt, Args&&... args) noexcept
      : P{fmt, std::forward<Args>(args)...} {}
};
template <typename... Args>
Debug(detail::format_string<Args...>, Args&&...) -> Debug<Args...>;

#define IGOR_DEBUG_PRINT(x) Igor::Debug("{} = {}", #x, x)  // NOLINT(cppcoreguidelines-macro-usage)

}  // namespace Igor

#endif  // IGOR_LOGGING_HPP_
