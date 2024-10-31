#ifndef ZAP_IO_POINTS_TO_NPY_HPP_
#define ZAP_IO_POINTS_TO_NPY_HPP_

#include <fstream>

#include <AD/ad.hpp>

#include "CellBased/Definitions.hpp"

#include "Igor/Logging.hpp"

namespace Zap::IO {

template <CellBased::Point2D_c PointType>
[[nodiscard]] auto points_to_npy(const std::string& filename,
                                 const std::vector<PointType>& points) noexcept -> bool {
  using namespace std::string_literals;
  static_assert(
      std::is_same_v<decltype(std::declval<PointType>().x), decltype(std::declval<PointType>().y)>,
      "Expected x- and y-component to have the same type.");
  using ActiveFloat  = std::remove_cvref_t<decltype(std::declval<PointType>().x)>;
  using PassiveFloat = std::conditional_t<ad::mode<ActiveFloat>::is_ad_type,
                                          typename ad::mode<ActiveFloat>::passive_t,
                                          ActiveFloat>;
  static_assert(std::is_floating_point_v<PassiveFloat>,
                "Expected PassiveFloat to be a C++ floating point type.");

  std::ofstream out(filename, std::ios::out | std::ios::binary);
  if (!out) {
    Igor::Warn(
        "Could not open file `{}` for writing in binary: {}.", filename, std::strerror(errno));
    return false;
  }

  // Write magic string
  constexpr std::streamsize magic_string_len = 6;
  if (!out.write("\x93NUMPY", magic_string_len)) {
    Igor::Warn("Could not write magic string to  {}: {}", filename, std::strerror(errno));
    return false;
  }

  // Write format version, use 1.0
  constexpr std::streamsize version_len = 2;
  if (!out.write("\x01\x00", version_len)) {
    Igor::Warn("Could not write version number to  {}: {}", filename, std::strerror(errno));
    return false;
  }

  // Length of length entry
  constexpr std::streamsize header_len_len = 2;

  // Create header
  // clang-format off
  std::string header =
    "{"s +
        "'descr': '<f"s + std::to_string(sizeof(PassiveFloat)) + "', "s +
        "'fortran_order': False, "s +
        "'shape': ("s + std::to_string(points.size()) + ", 2), "s +
    "}"s;
  // clang-format on

  // Pad header with spaces s.t. magic string, version, header length and header together are
  // divisible by 64
  for (auto preamble_len = magic_string_len + version_len + header_len_len + header.size() + 1;
       preamble_len % 64 != 0;
       preamble_len = magic_string_len + version_len + header_len_len + header.size() + 1) {
    header.push_back('\x20');
  }
  header.push_back('\n');

  // Write header length
  assert(header.size() <= std::numeric_limits<uint16_t>::max());
  const auto header_len = static_cast<uint16_t>(header.size());
  if (!out.write(reinterpret_cast<const char*>(&header_len), header_len_len)) {  // NOLINT
    Igor::Warn("Could not write header length to  {}: {}", filename, std::strerror(errno));
    return false;
  }

  // Write header
  if (!out.write(header.data(), static_cast<std::streamsize>(header.size()))) {
    Igor::Warn("Could not write header to  {}: {}", filename, std::strerror(errno));
    return false;
  }

  // Write data
  for (const auto& p : points) {
    const std::array<PassiveFloat, 2> p_value{ad::value(p.x), ad::value(p.y)};
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    if (!out.write(reinterpret_cast<const char*>(p_value.data()), sizeof(p_value))) {
      Igor::Warn("Could not write data to  {}: {}", filename, std::strerror(errno));
    }
  }

  return true;
}

}  // namespace Zap::IO

#endif  // ZAP_IO_POINTS_TO_NPY_HPP_
