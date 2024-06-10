#ifndef IGOR_HPP_
#define IGOR_HPP_

#include <cxxabi.h>
#include <format>
#include <iostream>
#include <memory>
#include <source_location>
#include <string>
#include <type_traits>
#include <utility>
#include <version>

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
};

consteval auto level_stream(Level level) noexcept -> std::ostream& {
  switch (level) {
    case Level::INFO:
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

}  // namespace Igor

#endif  // IGOR_HPP_
