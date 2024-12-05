#ifndef ZAP_CELL_BASED_DEFINITIONS_HPP_
#define ZAP_CELL_BASED_DEFINITIONS_HPP_

#include <bit>
#include <concepts>
#include <cstddef>

#include <AD/ad.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/fwd.hpp>
namespace mp = boost::multiprecision;

#include "Igor/Logging.hpp"

namespace Zap::CellBased {

// - Type traits for passive and active types ------------------------------------------------------
template <typename T>
inline constexpr bool is_passive_real = std::is_floating_point_v<T> || mp::is_number<T>::value;

template <typename T>
inline constexpr bool is_active_real =
    ad::mode<T>::is_ad_type && is_passive_real<typename ad::mode<T>::passive_t>;

template <typename T>
inline constexpr bool is_real = is_passive_real<T> || is_active_real<T>;

enum class Orientation : uint8_t { X, Y, FREE };
using enum Orientation;

// - Tolerances for given floating point accuracy --------------------------------------------------
template <typename Float>
requires(is_real<Float> && !mp::is_number<Float>::value)
[[nodiscard]] constexpr auto EPS() noexcept;

template <>
[[nodiscard]] constexpr auto EPS<float>() noexcept {
  return 1e-6F;
}
template <>
[[nodiscard]] constexpr auto EPS<double>() noexcept {
  return 1e-8;
}

template <>
[[nodiscard]] constexpr auto EPS<typename ad::gt1s<float>::type>() noexcept {
  return EPS<float>();
}
template <>
[[nodiscard]] constexpr auto EPS<typename ad::gt1s<double>::type>() noexcept {
  return EPS<double>();
}

template <typename Float>
requires mp::is_number<Float>::value
[[nodiscard]] constexpr auto EPS() noexcept {
  return mp::sqrt(std::numeric_limits<Float>::epsilon());
}

// - Check if two floating point values are approximately equal ------------------------------------
template <typename Float>
[[nodiscard]] constexpr auto approx_eq(Float a, Float b, Float tol = EPS<Float>()) noexcept
    -> bool {
  return a >= b ? (a - b <= tol) : (b - a <= tol);
};

// - Sign function ---------------------------------------------------------------------------------
template <typename T>
[[nodiscard]] constexpr auto sign(const T& x) noexcept -> T {
  // return static_cast<T>(x > EPS<T>()) - static_cast<T>(x < -EPS<T>());
  return static_cast<T>(x > T{0}) - static_cast<T>(x < T{0});
}

// - Points in grid and simulation coordinates -----------------------------------------------------
template <typename PointType>
concept Point2D_c = requires(PointType t) {
  t.x;
  t.y;

  { PointType::Zero() } -> std::same_as<PointType>;
  { t.norm() } -> std::same_as<decltype(t.x)>;
  { t.normalized() } -> std::same_as<PointType>;

  { std::declval<PointType>().dot(std::declval<PointType>()) } -> std::same_as<decltype(t.x)>;

  { std::declval<PointType>() + std::declval<PointType>() } -> std::same_as<PointType>;
  { std::declval<PointType>() - std::declval<PointType>() } -> std::same_as<PointType>;
  { std::declval<decltype(t.x)>() * std::declval<PointType>() } -> std::same_as<PointType>;
  { std::declval<PointType>() * std::declval<decltype(t.x)>() } -> std::same_as<PointType>;
  { std::declval<PointType>() / std::declval<decltype(t.x)>() } -> std::same_as<PointType>;

  { std::declval<PointType>() += std::declval<PointType>() } -> std::same_as<PointType&>;
  { std::declval<PointType>() -= std::declval<PointType>() } -> std::same_as<PointType&>;
  { std::declval<PointType>() *= std::declval<decltype(t.x)>() } -> std::same_as<PointType&>;
  { std::declval<PointType>() /= std::declval<decltype(t.x)>() } -> std::same_as<PointType&>;
};

// NOTE: GridCoord and SimCoord are functionally the same type, but we want the C++ type checker to
//       make sure that we never mix simulation- and grid-coordinates. Therefore they are defined to
//       be two different types that both work exactly the same by using the Preprocessor.
// NOLINTBEGIN(cppcoreguidelines-macro-usage,bugprone-macro-parentheses)
#define DEF_POINT2D_TYPE(name)                                                                     \
  template <typename Scalar>                                                                       \
  struct name {                                                                                    \
    Scalar x;                                                                                      \
    Scalar y;                                                                                      \
                                                                                                   \
    constexpr name() noexcept = default;                                                           \
    constexpr name(Scalar x_, Scalar y_) noexcept                                                  \
        : x(std::move(x_)),                                                                        \
          y(std::move(y_)) {}                                                                      \
    constexpr name(const name& other) noexcept                    = default;                       \
    constexpr name(name&& other) noexcept                         = default;                       \
    constexpr auto operator=(const name& other) noexcept -> name& = default;                       \
    constexpr auto operator=(name&& other) noexcept -> name&      = default;                       \
    constexpr ~name() noexcept                                    = default;                       \
    template <typename PassiveScalar>                                                              \
    requires(ad::mode<Scalar>::is_ad_type &&                                                       \
             std::is_same_v<PassiveScalar, typename ad::mode<Scalar>::passive_t>)                  \
    constexpr name(const name<PassiveScalar>& other) noexcept                                      \
        : x(other.x),                                                                              \
          y(other.y) {}                                                                            \
                                                                                                   \
    [[nodiscard]] static constexpr auto Zero() noexcept -> name { return {0, 0}; }                 \
                                                                                                   \
    [[nodiscard]] constexpr auto norm() const noexcept -> Scalar {                                 \
      return std::sqrt(x * x + y * y);                                                             \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] constexpr auto normalized() const noexcept -> name {                             \
      const auto n = norm();                                                                       \
      IGOR_ASSERT(std::abs(n) >= EPS<Scalar>(), "Cannot normalize null vector.");                  \
      return *this / norm();                                                                       \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] constexpr auto dot(const name& other) const noexcept -> Scalar {                 \
      return x * other.x + y * other.y;                                                            \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] friend constexpr auto operator+(const name& lhs, const name& rhs) noexcept       \
        -> name {                                                                                  \
      return {lhs.x + rhs.x, lhs.y + rhs.y};                                                       \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] friend constexpr auto operator-(const name& lhs, const name& rhs) noexcept       \
        -> name {                                                                                  \
      return {lhs.x - rhs.x, lhs.y - rhs.y};                                                       \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] friend constexpr auto operator*(const name& lhs, Scalar scalar) noexcept         \
        -> name {                                                                                  \
      return {lhs.x * scalar, lhs.y * scalar};                                                     \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] friend constexpr auto operator*(Scalar scalar, const name& rhs) noexcept         \
        -> name {                                                                                  \
      return {scalar * rhs.x, scalar * rhs.y};                                                     \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] friend constexpr auto operator/(const name& lhs, Scalar scalar) noexcept         \
        -> name {                                                                                  \
      return {lhs.x / scalar, lhs.y / scalar};                                                     \
    }                                                                                              \
                                                                                                   \
    constexpr auto operator+=(const name& other) noexcept -> name& {                               \
      x += other.x;                                                                                \
      y += other.y;                                                                                \
      return *this;                                                                                \
    }                                                                                              \
                                                                                                   \
    constexpr auto operator-=(const name& other) noexcept -> name& {                               \
      x -= other.x;                                                                                \
      y -= other.y;                                                                                \
      return *this;                                                                                \
    }                                                                                              \
                                                                                                   \
    constexpr auto operator*=(Scalar scalar) noexcept -> name& {                                   \
      x *= scalar;                                                                                 \
      y *= scalar;                                                                                 \
      return *this;                                                                                \
    }                                                                                              \
                                                                                                   \
    constexpr auto operator/=(Scalar scalar) noexcept -> name& {                                   \
      x /= scalar;                                                                                 \
      y /= scalar;                                                                                 \
      return *this;                                                                                \
    }                                                                                              \
  }

#define DEFINE_TEMPLATE_CHECK(name)                                                                \
  template <typename T>                                                                            \
  struct is_##name : std::false_type {};                                                           \
  template <typename T>                                                                            \
  struct is_##name<name<T>> : std::true_type {};                                                   \
  template <typename T>                                                                            \
  constexpr bool is_##name##_v = is_##name<T>::value

#define DEF_POINT2D_FORMATTER(name)                                                                \
  template <typename T, typename CharT>                                                            \
  struct fmt::formatter<name<T>, CharT> {                                                          \
    template <typename ParseContext>                                                               \
    static constexpr auto parse(ParseContext& ctx) noexcept {                                      \
      return ctx.begin();                                                                          \
    }                                                                                              \
    template <typename FormatContext>                                                              \
    static constexpr auto format(const name<T>& p, FormatContext& ctx) noexcept {                  \
      return fmt::format_to(ctx.out(), "[{}, {}]", p.x, p.y);                                      \
    }                                                                                              \
  }
// NOLINTEND(cppcoreguidelines-macro-usage,bugprone-macro-parentheses)

DEF_POINT2D_TYPE(GenCoord);

DEF_POINT2D_TYPE(GridCoord);
DEFINE_TEMPLATE_CHECK(GridCoord);

DEF_POINT2D_TYPE(SimCoord);
DEFINE_TEMPLATE_CHECK(SimCoord);

static_assert(std::is_trivial_v<GenCoord<double>>,
              "GenCoord type should be trivial for trivial Scalar.");
static_assert(std::is_trivial_v<GridCoord<double>>,
              "GridCoord type should be trivial for trivial Scalar.");
static_assert(std::is_trivial_v<SimCoord<double>>,
              "SimCoord type should be trivial for trivial Scalar.");

static_assert(is_GridCoord_v<GridCoord<double>>,
              "is_GridCoord_v must evaluate to true for GridCoord.");
static_assert(!is_SimCoord_v<GridCoord<double>>,
              "is_SimCoord_v must evaluate to false for GridCoord.");
static_assert(Point2D_c<GridCoord<double>>,
              "GridCoord must fulfill the Point2D concepts requirements.");

static_assert(!is_GridCoord_v<SimCoord<double>>,
              "is_GridCoord_v must evaluate to false for SimCoord.");
static_assert(is_SimCoord_v<SimCoord<double>>, "is_SimCoord_v must evaluate to true for SimCoord.");
static_assert(Point2D_c<SimCoord<double>>,
              "GridCoord must fulfill the Point2D concepts requirements.");

enum class CoordType : uint8_t { GRID_C, SIM_C };
using enum CoordType;

template <typename Float, CoordType COORD_TYPE>
using CoordType2PointType =
    std::conditional_t<COORD_TYPE == GRID_C, GridCoord<Float>, SimCoord<Float>>;

template <typename PointType>
requires(is_GridCoord_v<PointType> || is_SimCoord_v<PointType>)
constexpr CoordType PointType2CoordType = is_SimCoord_v<PointType> ? SIM_C : GRID_C;

// -------------------------------------------------------------------------------------------------
enum Side : uint8_t {
  BOTTOM = 0b0001,
  RIGHT  = 0b0010,
  TOP    = 0b0100,
  LEFT   = 0b1000,
  ALL    = LEFT | RIGHT | BOTTOM | TOP,
};

[[nodiscard]] constexpr auto make_cut_type(Side entry, Side exit) -> uint8_t {
  return static_cast<uint8_t>((exit << 4) | entry);
}

[[nodiscard]] constexpr auto count_sides(Side side) noexcept {
  return std::popcount(std::to_underlying(side));
}

// -------------------------------------------------------------------------------------------------
enum class ExtendType : uint8_t { NONE, NEAREST, MAX };

}  // namespace Zap::CellBased

template <>
struct fmt::formatter<Zap::CellBased::Side> {
  template <typename ParseContext>
  static constexpr auto parse(ParseContext& ctx) noexcept {
    return ctx.begin();
  }
  template <typename FormatContext>
  static constexpr auto format(Zap::CellBased::Side side, FormatContext& ctx) noexcept {
    if (side == Zap::CellBased::ALL) { return fmt::format_to(ctx.out(), "ALL"); }
    if (side == 0) { return fmt::format_to(ctx.out(), "NONE"); }

    std::string s;
    if ((side & Zap::CellBased::BOTTOM) > 0) { s += "BOTTOM"; }
    if ((side & Zap::CellBased::RIGHT) > 0) {
      if (!s.empty()) { s += " | "; }
      s += "RIGHT";
    }
    if ((side & Zap::CellBased::TOP) > 0) {
      if (!s.empty()) { s += " | "; }
      s += "TOP";
    }
    if ((side & Zap::CellBased::LEFT) > 0) {
      if (!s.empty()) { s += " | "; }
      s += "LEFT";
    }
    return fmt::format_to(ctx.out(), "{}", s);
  }
};

DEF_POINT2D_FORMATTER(Zap::CellBased::GenCoord);
DEF_POINT2D_FORMATTER(Zap::CellBased::GridCoord);
DEF_POINT2D_FORMATTER(Zap::CellBased::SimCoord);

#undef DEF_POINT2D_TYPE
#undef DEFINE_TEMPLATE_CHECK
#undef DEF_POINT2D_FORMATTER

#endif  // ZAP_CELL_BASED_DEFINITIONS_HPP_
