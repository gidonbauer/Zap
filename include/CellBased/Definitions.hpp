#ifndef ZAP_CELL_BASED_DEFINITIONS_HPP_
#define ZAP_CELL_BASED_DEFINITIONS_HPP_

#include <concepts>
#include <cstddef>
#include <format>

#include <Eigen/Core>

namespace Zap::CellBased {

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
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define DEF_POINT2D_TYPE(name)                                                                     \
  template <typename Scalar>                                                                       \
  struct name {                                                                                    \
    Scalar x;                                                                                      \
    Scalar y;                                                                                      \
                                                                                                   \
    [[nodiscard]] static constexpr auto Zero() noexcept -> name { return {.x = 0, .y = 0}; }       \
                                                                                                   \
    [[nodiscard]] constexpr auto norm() const noexcept -> Scalar {                                 \
      return std::sqrt(x * x + y * y);                                                             \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] constexpr auto normalized() const noexcept -> name { return *this / norm(); }    \
                                                                                                   \
    [[nodiscard]] constexpr auto dot(const name& other) const noexcept -> Scalar {                 \
      return x * other.x + y * other.y;                                                            \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] friend constexpr auto operator+(const name& lhs, const name& rhs) noexcept       \
        -> name {                                                                                  \
      return {.x = lhs.x + rhs.x, .y = lhs.y + rhs.y};                                             \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] friend constexpr auto operator-(const name& lhs, const name& rhs) noexcept       \
        -> name {                                                                                  \
      return {.x = lhs.x - rhs.x, .y = lhs.y - rhs.y};                                             \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] friend constexpr auto operator*(const name& lhs, Scalar scalar) noexcept         \
        -> name {                                                                                  \
      return {.x = lhs.x * scalar, .y = lhs.y * scalar};                                           \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] friend constexpr auto operator*(Scalar scalar, const name& rhs) noexcept         \
        -> name {                                                                                  \
      return {.x = scalar * rhs.x, .y = scalar * rhs.y};                                           \
    }                                                                                              \
                                                                                                   \
    [[nodiscard]] friend constexpr auto operator/(const name& lhs, Scalar scalar) noexcept         \
        -> name {                                                                                  \
      return {.x = lhs.x / scalar, .y = lhs.y / scalar};                                           \
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

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define DEF_POINT2D_FORMATTER(name)                                                                \
  template <typename T, typename CharT>                                                            \
  struct formatter<name<T>, CharT> {                                                               \
    template <typename ParseContext>                                                               \
    static constexpr auto parse(ParseContext& ctx) noexcept {                                      \
      return ctx.begin();                                                                          \
    }                                                                                              \
    template <typename FormatContext>                                                              \
    static constexpr auto format(const name<T>& p, FormatContext& ctx) noexcept {                  \
      return std::format_to(ctx.out(), "[{}, {}]", p.x, p.y);                                      \
    }                                                                                              \
  }

DEF_POINT2D_TYPE(GridCoord);
static_assert(Point2D_c<GridCoord<double>>,
              "GridCoord must fulfill the Point2D concepts requirements.");

DEF_POINT2D_TYPE(SimCoord);
static_assert(Point2D_c<SimCoord<double>>,
              "GridCoord must fulfill the Point2D concepts requirements.");
}  // namespace Zap::CellBased

namespace std {

DEF_POINT2D_FORMATTER(Zap::CellBased::GridCoord);
DEF_POINT2D_FORMATTER(Zap::CellBased::SimCoord);

}  // namespace std

namespace Zap::CellBased {

enum class CoordType : uint8_t { GRID, SIM };

template <typename Float, CoordType COORD_TYPE>
using Coord_t =
    std::conditional_t<COORD_TYPE == CoordType::GRID, GridCoord<Float>, SimCoord<Float>>;

#undef DEF_POINT2D_TYPE
#undef DEF_POINT2D_FORMATTER

// - Generic 2D Point based on Eigen Vector --------------------------------------------------------
enum : size_t { X, Y, POINT_SIZE };  // NOLINT

template <typename T>
using Point = Eigen::Vector<T, POINT_SIZE>;

// -------------------------------------------------------------------------------------------------
enum Side : int {  // NOLINT
  BOTTOM = 0b0001,
  RIGHT  = 0b0010,
  TOP    = 0b0100,
  LEFT   = 0b1000,
  ALL    = LEFT | RIGHT | BOTTOM | TOP,
};

// - Tolerances for given floating point accuracy --------------------------------------------------
template <std::floating_point Float>
inline constexpr Float EPS;
template <>
inline constexpr float EPS<float> = 1e-6f;
template <>
inline constexpr double EPS<double> = 1e-8;

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_DEFINITIONS_HPP_
