#ifndef ZAP_IO_POINTS_TO_NPY_HPP_
#define ZAP_IO_POINTS_TO_NPY_HPP_

#include <concepts>
#include <fstream>

#include <AD/ad.hpp>

#include "CellBased/Definitions.hpp"

#include "Igor/Logging.hpp"

namespace Zap::IO {

namespace detail {

// -------------------------------------------------------------------------------------------------
template <std::floating_point Float>
[[nodiscard]] auto write_npy_header(std::ostream& out,
                                    const std::string& filename,
                                    bool row_major,
                                    const std::vector<size_t>& shape) noexcept -> bool {
  using namespace std::string_literals;

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
  std::string header = "{"s;

  // Data type
  header += "'descr': '<f"s + std::to_string(sizeof(Float)) + "', "s;
  // Data order, Fortran order (column major) or C order (row major)
  header += "'fortran_order': "s + (row_major ? "False"s : "True"s) + ", "s;
  // Data shape
  IGOR_ASSERT(shape.size() >= 1UZ, "Shape must have at least one entry but has {}", shape.size());
  header += "'shape': ("s;
  for (size_t s : shape) {
    header += std::to_string(s) + ", "s;
  }
  header += "), "s;

  header += "}"s;

  // Pad header with spaces s.t. magic string, version, header length and header together are
  // divisible by 64
  for (auto preamble_len = magic_string_len + version_len + header_len_len + header.size() + 1;
       preamble_len % 64 != 0;
       preamble_len = magic_string_len + version_len + header_len_len + header.size() + 1) {
    header.push_back('\x20');
  }
  header.push_back('\n');

  // Write header length
  IGOR_ASSERT(
      header.size() <= std::numeric_limits<uint16_t>::max(),
      "Size cannot be larger than the max for an unsigned 16-bit integer as it is stored in one.");
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

  return true;
}

}  // namespace detail

// -------------------------------------------------------------------------------------------------
template <std::floating_point Float>
[[nodiscard]] auto vector_to_npy(const std::string& filename,
                                 const std::vector<Float>& vec) noexcept -> bool {
  std::ofstream out(filename, std::ios::out | std::ios::binary);
  if (!out) {
    Igor::Warn(
        "Could not open file `{}` for writing in binary: {}.", filename, std::strerror(errno));
    return false;
  }

  if (!detail::write_npy_header<Float>(out, filename, true, {vec.size()})) {
    Igor::Warn("Could not write header.");
    return false;
  }

  // Write data
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  if (!out.write(reinterpret_cast<const char*>(vec.data()),
                 static_cast<std::streamsize>(sizeof(Float) * vec.size()))) {
    Igor::Warn("Could not write data to {}: {}", filename, std::strerror(errno));
  }

  return true;
}

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat>
requires ad::mode<ActiveFloat>::is_ad_type
[[nodiscard]] auto vector_to_npy(const std::string& primal_filename,
                                 const std::string& derivative_filename,
                                 const std::vector<ActiveFloat>& vec) noexcept -> bool {
  using PassiveFloat = ad::mode<ActiveFloat>::passive_t;
  static_assert(std::is_floating_point_v<PassiveFloat>,
                "Passive type has to be a C++ buildin floating point type.");

  std::vector<PassiveFloat> passive_vec(vec.size());
  std::transform(std::cbegin(vec), std::cend(vec), std::begin(passive_vec), [](auto& value) {
    return ad::value(value);
  });
  const auto wrote_primal = vector_to_npy(primal_filename, passive_vec);

  std::transform(std::cbegin(vec), std::cend(vec), std::begin(passive_vec), [](auto& value) {
    return ad::derivative(value);
  });
  const auto wrote_derivative = vector_to_npy(derivative_filename, passive_vec);

  return wrote_primal && wrote_derivative;
}

// -------------------------------------------------------------------------------------------------
template <std::floating_point Float>
[[nodiscard]] auto matrix_to_npy(const std::string& filename,
                                 const std::vector<std::vector<Float>>& mat) noexcept -> bool {
  IGOR_ASSERT(mat.size() >= 1UZ, "Size must be at least one but is {}", mat.size());
  const auto n_cols = mat[0].size();
  IGOR_ASSERT(std::all_of(std::cbegin(mat),
                          std::cend(mat),
                          [n_cols](auto& vec) { return vec.size() == n_cols; }),
              "All rows must have the same number of columns.");

  std::ofstream out(filename, std::ios::out | std::ios::binary);
  if (!out) {
    Igor::Warn(
        "Could not open file `{}` for writing in binary: {}.", filename, std::strerror(errno));
    return false;
  }

  if (!detail::write_npy_header<Float>(out, filename, true, {mat.size(), mat[0].size()})) {
    Igor::Warn("Could not write header.");
    return false;
  }

  // Write data
  for (const auto& vec : mat) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    if (!out.write(reinterpret_cast<const char*>(vec.data()),
                   static_cast<std::streamsize>(sizeof(Float) * vec.size()))) {
      Igor::Warn("Could not write data to {}: {}", filename, std::strerror(errno));
    }
  }

  return true;
}

// -------------------------------------------------------------------------------------------------
template <typename ActiveFloat>
requires ad::mode<ActiveFloat>::is_ad_type
[[nodiscard]] auto matrix_to_npy(const std::string& primal_filename,
                                 const std::string& derivative_filename,
                                 const std::vector<std::vector<ActiveFloat>>& mat) noexcept
    -> bool {
  IGOR_ASSERT(mat.size() >= 1UZ, "Size must be at least one but is {}", mat.size());
  const auto n_cols = mat[0].size();
  IGOR_ASSERT(std::all_of(std::cbegin(mat),
                          std::cend(mat),
                          [n_cols](auto& vec) { return vec.size() == n_cols; }),
              "All rows must have the same number of columns.");

  using PassiveFloat = ad::mode<ActiveFloat>::passive_t;
  static_assert(std::is_floating_point_v<PassiveFloat>,
                "Passive type has to be a C++ buildin floating point type.");

  std::vector<std::vector<PassiveFloat>> passive_mat(mat.size(), std::vector<PassiveFloat>(n_cols));

  for (size_t i = 0; i < mat.size(); ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      passive_mat[i][j] = ad::value(mat[i][j]);
    }
  }
  const auto wrote_primal = matrix_to_npy(primal_filename, passive_mat);

  for (size_t i = 0; i < mat.size(); ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      passive_mat[i][j] = ad::derivative(mat[i][j]);
    }
  }
  const auto wrote_derivative = matrix_to_npy(derivative_filename, passive_mat);

  return wrote_primal && wrote_derivative;
}

// -------------------------------------------------------------------------------------------------
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

  if (!detail::write_npy_header<PassiveFloat>(out, filename, true, {points.size(), 2UZ})) {
    Igor::Warn("Could not write header.");
    return false;
  }

  // Write data
  for (const auto& p : points) {
    const std::array<PassiveFloat, 2> p_value{ad::value(p.x), ad::value(p.y)};
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    if (!out.write(reinterpret_cast<const char*>(p_value.data()), sizeof(p_value))) {
      Igor::Warn("Could not write data to {}: {}", filename, std::strerror(errno));
    }
  }

  return true;
}

}  // namespace Zap::IO

#endif  // ZAP_IO_POINTS_TO_NPY_HPP_
