#ifndef ZAP_IO_COMMON_HPP_
#define ZAP_IO_COMMON_HPP_

#include <string_view>
#include <type_traits>

#include "Igor.hpp"

namespace Zap::IO {

// -------------------------------------------------------------------------------------------------
template <typename Scalar>
[[nodiscard]] static consteval auto dtype_short() noexcept -> std::string_view {
  using namespace std::string_view_literals;
  static_assert(!std::is_reference_v<Scalar>);
  using S = std::remove_cvref_t<Scalar>;
  if constexpr (std::is_same_v<S, std::uint8_t>) {
    return " u8"sv;
  }
  if constexpr (std::is_same_v<S, std::uint16_t>) {
    return "u16"sv;
  }
  if constexpr (std::is_same_v<S, std::uint32_t>) {
    return "u32"sv;
  }
  if constexpr (std::is_same_v<S, std::uint64_t>) {
    return "u64"sv;
  }
  if constexpr (std::is_same_v<S, std::int8_t>) {
    return " i8"sv;
  }
  if constexpr (std::is_same_v<S, std::int16_t>) {
    return "i16"sv;
  }
  if constexpr (std::is_same_v<S, std::int32_t>) {
    return "i32"sv;
  }
  if constexpr (std::is_same_v<S, std::int64_t>) {
    return "i64"sv;
  }
  if constexpr (std::is_same_v<S, float>) {
    return "f32"sv;
  }
  if constexpr (std::is_same_v<S, double>) {
    return "f64"sv;
  }
  Igor::Panic("Unknown scalar type `{}`", Igor::type_name<Scalar>());
  std::unreachable();
}

}  // namespace Zap::IO

#endif  // ZAP_IO_COMMON_HPP_
