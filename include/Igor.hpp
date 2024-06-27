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

#ifndef IGOR_HPP_
#define IGOR_HPP_

#include <cassert>
#include <chrono>
#include <format>
#include <iomanip>
#include <iostream>
#include <memory>
#include <source_location>
#include <string>
#include <type_traits>
#include <utility>
#include <version>

#include <cxxabi.h>

static_assert(__cpp_lib_source_location >= 201907L, "Requires source_location.");
// static_assert(__cpp_lib_format >= 201907L, "Requires std::format.");

namespace Igor {

enum class ExitCode : int {  // NOLINT(performance-enum-size)
  NO_EARLY_EXIT,
  PANIC,
  TODO,
};

namespace detail {

[[nodiscard]] constexpr auto error_loc(
    const std::source_location loc = std::source_location::current()) noexcept -> std::string {
  try {
    return std::format("`{}` (\033[95m{}:{}:{}\033[0m)",
                       loc.function_name(),
                       loc.file_name(),
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
};

consteval auto level_stream(Level level) noexcept -> std::ostream& {
  switch (level) {
    case Level::INFO:
    case Level::TIME:
      return std::cout;
    case Level::WARN:
    case Level::TODO:
    case Level::PANIC:
    case Level::DEBUG:
      return std::cerr;
  }
}

consteval auto level_repr(Level level) noexcept {
  switch (level) {
    case Level::INFO:
      return "\033[32m[INFO]\033[0m ";
    case Level::WARN:
      return "\033[33m[WARN]\033[0m ";
    case Level::TODO:
      return "\033[35m[TODO]\033[0m ";
    case Level::PANIC:
      return "\033[31m[ERROR]\033[0m ";
    case Level::DEBUG:
      return "\033[94m[DEBUG]\033[0m ";
    case Level::TIME:
      return "\033[94m[TIME]\033[0m ";
  }
}

template <Level level, ExitCode exit_code, typename... Args>
class Print {
 protected:
  constexpr Print(std::format_string<Args...> fmt, Args&&... args) noexcept {
    auto& out = level_stream(level);
    out << level_repr(level) << std::format(fmt, std::forward<Args>(args)...) << '\n';
    if constexpr (exit_code != ExitCode::NO_EARLY_EXIT) {
      std::exit(static_cast<int>(exit_code));
    }
  }

  constexpr Print(const std::source_location loc,
                  std::format_string<Args...> fmt,
                  Args&&... args) noexcept {
    auto& out = level_stream(level);
    out << level_repr(level) << error_loc(loc) << ": "
        << std::format(fmt, std::forward<Args>(args)...) << '\n';
    if constexpr (exit_code != ExitCode::NO_EARLY_EXIT) {
      std::exit(static_cast<int>(exit_code));
    }
  }
};

// -------------------------------------------------------------------------------------------------
template <typename... Args>
class [[maybe_unused]] Time final
    : detail::Print<detail::Level::TIME, ExitCode::NO_EARLY_EXIT, Args...> {
  using P = detail::Print<detail::Level::TIME, ExitCode::NO_EARLY_EXIT, Args...>;

 public:
  constexpr Time(std::format_string<Args...> fmt, Args&&... args) noexcept
      : P{fmt, std::forward<Args>(args)...} {}
};
template <typename... Args>
Time(std::format_string<Args...>, Args&&...) -> Time<Args...>;

}  // namespace detail

// -------------------------------------------------------------------------------------------------
template <typename... Args>
class [[maybe_unused]] Info final
    : detail::Print<detail::Level::INFO, ExitCode::NO_EARLY_EXIT, Args...> {
  using P = detail::Print<detail::Level::INFO, ExitCode::NO_EARLY_EXIT, Args...>;

 public:
  constexpr Info(std::format_string<Args...> fmt, Args&&... args) noexcept
      : P{fmt, std::forward<Args>(args)...} {}
};
template <typename... Args>
Info(std::format_string<Args...>, Args&&...) -> Info<Args...>;

// -------------------------------------------------------------------------------------------------
template <typename... Args>
class [[maybe_unused]] Warn final
    : detail::Print<detail::Level::WARN, ExitCode::NO_EARLY_EXIT, Args...> {
  using P = detail::Print<detail::Level::WARN, ExitCode::NO_EARLY_EXIT, Args...>;

 public:
  constexpr Warn(std::format_string<Args...> fmt,
                 Args&&... args,
                 const std::source_location loc = std::source_location::current()) noexcept
      : P{loc, fmt, std::forward<Args>(args)...} {}
};
template <typename... Args>
Warn(std::format_string<Args...>, Args&&...) -> Warn<Args...>;

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
      std::format_string<Args...> fmt,
      Args&&... args,
      const std::source_location loc = std::source_location::current()) noexcept
      : P{loc, fmt, std::forward<Args>(args)...} {
    std::unreachable();
  }
};
template <typename... Args>
Todo(std::format_string<Args...>, Args&&...) -> Todo<Args...>;

// -------------------------------------------------------------------------------------------------
template <typename... Args>
class [[maybe_unused]] Panic final : detail::Print<detail::Level::PANIC, ExitCode::PANIC, Args...> {
  using P = detail::Print<detail::Level::PANIC, ExitCode::PANIC, Args...>;

 public:
  [[noreturn]] constexpr Panic(
      std::format_string<Args...> fmt,
      Args&&... args,
      const std::source_location loc = std::source_location::current()) noexcept
      : P{loc, fmt, std::forward<Args>(args)...} {
    std::unreachable();
  }
};
template <typename... Args>
Panic(std::format_string<Args...>, Args&&...) -> Panic<Args...>;

// -------------------------------------------------------------------------------------------------
template <typename... Args>
class [[maybe_unused]] Debug final
    : detail::Print<detail::Level::DEBUG, ExitCode::NO_EARLY_EXIT, Args...> {
  using P = detail::Print<detail::Level::DEBUG, ExitCode::NO_EARLY_EXIT, Args...>;

 public:
  constexpr Debug(std::format_string<Args...> fmt, Args&&... args) noexcept
      : P{fmt, std::forward<Args>(args)...} {}
};
template <typename... Args>
Debug(std::format_string<Args...>, Args&&...) -> Debug<Args...>;

#define IGOR_DEBUG_PRINT(x) Igor::Debug("{} = {}", #x, x)  // NOLINT(cppcoreguidelines-macro-usage)

// -------------------------------------------------------------------------------------------------
#define IGOR_STRINGIFY(s) IGOR_XSTRINGIFY(s)  // NOLINT(cppcoreguidelines-macro-usage)
#define IGOR_XSTRINGIFY(s) #s                 // NOLINT(cppcoreguidelines-macro-usage)

// -------------------------------------------------------------------------------------------------
template <typename T>
[[nodiscard]] constexpr auto type_name() -> std::string {
  using namespace std::string_literals;

  int status;
  constexpr auto free_deleter = [](void* p) constexpr noexcept {
    std::free(p);  // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc)
  };
  std::unique_ptr<char, decltype(free_deleter)> name_cstr{
      abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status), free_deleter};

  if (status != 0 || name_cstr == nullptr) {
    switch (status) {
      case -1:
        throw std::runtime_error(
            "Demagleing failed with status -1: A memory allocation failure occurred.");
      case -2:
        throw std::runtime_error("Demagleing failed with status -2: mangled_name is not a valid "
                                 "name under the C++ ABI mangling rules.");
      case -3:
        throw std::runtime_error(
            "Demagleing failed with status -3: One of the arguments is invalid.");
      default:
        throw std::runtime_error("Demagleing failed with unknown status "s +
                                 std::to_string(status) + "."s);
    }
  }

  std::string name{name_cstr.get()};
  if (std::is_volatile_v<std::remove_reference_t<T>>) {
    name += " volatile"s;
  }
  if (std::is_const_v<std::remove_reference_t<T>>) {
    name += " const"s;
  }
  if (std::is_lvalue_reference_v<T>) {
    name += "&"s;
  }
  if (std::is_rvalue_reference_v<T>) {
    name += "&&"s;
  }

  return name;
}

template <typename T>
[[nodiscard]] constexpr auto type_name(T /*ignored*/) -> std::string {
  return type_name<T>();
}

// - Times the duration of the scope in wall-clock time --------------------------------------------
class ScopeTimer {
  std::string m_scope_name;
  std::chrono::high_resolution_clock::time_point m_t_begin;

 public:
  [[nodiscard]] ScopeTimer(std::string scope_name = "Scope") noexcept
      : m_scope_name(std::move(scope_name)),
        m_t_begin(std::chrono::high_resolution_clock::now()) {}

  ScopeTimer(const ScopeTimer& other) noexcept                    = delete;
  ScopeTimer(ScopeTimer&& other) noexcept                         = delete;
  auto operator=(const ScopeTimer& other) noexcept -> ScopeTimer& = delete;
  auto operator=(ScopeTimer&& other) noexcept -> ScopeTimer&      = delete;
  ~ScopeTimer() noexcept {
    const auto t_duration =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - m_t_begin);
    detail::Time("{} took {}.", m_scope_name, t_duration);
  }
};

#define IGOR_COMBINE1(X, Y) X##Y                // NOLINT
#define IGOR_COMBINE(X, Y) IGOR_COMBINE1(X, Y)  // NOLINT

// NOLINTNEXTLINE
#define IGOR_TIME_SCOPE(...)                                                                       \
  if (const auto IGOR_COMBINE(IGOR__SCOPE__TIMER__NAME__, __LINE__) =                              \
          Igor::ScopeTimer{__VA_ARGS__};                                                           \
      true)

// - Simple command line progress bar --------------------------------------------------------------
class ProgressBar {
  std::size_t m_max_progress;
  std::size_t m_length;
  std::size_t m_progress = 0UZ;

  enum : char {
    DONE_CHAR     = '#',
    NOT_DONE_CHAR = '.',
  };

 public:
  constexpr ProgressBar(std::size_t max_progress, std::size_t length) noexcept
      : m_max_progress(max_progress),
        m_length(length - 5UZ) {
    assert(length < std::numeric_limits<int>::max());
  }

  constexpr void update() noexcept {
    m_progress = std::min(m_progress + 1, m_max_progress);
    show();
  }

  void show() const noexcept {
    const auto done_length = (m_length * m_progress) / m_max_progress;
    const auto done_prct   = (100UZ * m_progress) / m_max_progress;
    std::cout << "\r[";
    std::cout << std::string(done_length, DONE_CHAR);
    std::cout << std::string(m_length - done_length, NOT_DONE_CHAR);
    std::cout << "] ";
    std::cout << std::setw(3) << done_prct << "%" << std::flush;
  }
};

// - Transforms memory in bytes to a human readable string -----------------------------------------
[[nodiscard]] auto memory_to_string(uint64_t mem_in_bytes) noexcept -> std::string {
  using namespace std::string_literals;
  constexpr uint64_t step_factor = 1024;

  if (mem_in_bytes < step_factor) {
    return std::to_string(mem_in_bytes) + " B"s;
  }
  if (mem_in_bytes < step_factor * step_factor) {
    return std::to_string(mem_in_bytes / step_factor) + " kB"s;
  }
  if (mem_in_bytes < step_factor * step_factor * step_factor) {
    return std::to_string(mem_in_bytes / (step_factor * step_factor)) + " MB"s;
  }
  return std::to_string(mem_in_bytes / (step_factor * step_factor * step_factor)) + " GB"s;
}

}  // namespace Igor

#endif  // IGOR_HPP_
