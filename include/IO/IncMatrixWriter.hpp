#ifndef ZAP_IO_INC_MATRIX_WRITER_HPP_
#define ZAP_IO_INC_MATRIX_WRITER_HPP_

#include <cstdint>
#include <fstream>
#include <string_view>
#include <type_traits>

#include <Eigen/Dense>

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

namespace HeaderLayout {

using namespace std::string_view_literals;

constexpr auto MAGIC_STRING      = "\x24SHOCK"sv;  // NOLINT(modernize-raw-string-literal)
constexpr auto DTYPE_SIZE        = 3UZ;
constexpr auto ROWS_SIZE         = 8UZ;
constexpr auto COLS_SIZE         = 8UZ;
constexpr auto IS_ROW_MAJOR_SIZE = 1UZ;
constexpr auto HEADER_SIZE =
    MAGIC_STRING.size() + DTYPE_SIZE + ROWS_SIZE + COLS_SIZE + IS_ROW_MAJOR_SIZE;

}  // namespace HeaderLayout

// -------------------------------------------------------------------------------------------------
template <typename Scalar, int ROWS, int COLS, int OPTIONS>
class IncMatrixWriter {
  std::string m_filename;
  std::ofstream m_out;

  // -----------------------------------------------------------------------------------------------
  IncMatrixWriter(std::string filename)
      : m_filename(std::move(filename)),
        m_out(std::ofstream(m_filename, std::ios::out | std::ios::binary)) {}

 public:
  IncMatrixWriter(std::string filename, int64_t rows, int64_t cols, int8_t is_row_major)
      : IncMatrixWriter(std::move(filename)) {
    if (!m_out) {
      Igor::Warn("Could not open file `{}` for writing: ", m_filename, std::strerror(errno));
      throw std::runtime_error(
          std::format("Opening file `{}` failed: {}", m_filename, std::strerror(errno)));
    }

    if (!write_header(rows, cols, is_row_major)) {
      throw std::runtime_error("Writing the header failed.");
    }
  }

  // -----------------------------------------------------------------------------------------------
  IncMatrixWriter(std::string filename, const Eigen::Matrix<Scalar, ROWS, COLS, OPTIONS>& mat)
      : IncMatrixWriter(std::move(filename), mat.rows(), mat.cols(), mat.IsRowMajor) {}

  // -----------------------------------------------------------------------------------------------
 private:
  [[nodiscard]] auto
  write_header(int64_t rows, int64_t cols, int8_t is_row_major) noexcept -> bool {
    using namespace std::string_view_literals;

    if (!m_out.write(HeaderLayout::MAGIC_STRING.data(), HeaderLayout::MAGIC_STRING.size())) {
      Igor::Warn("Could not write magic string to `{}`.", m_filename);
      return false;
    }

    {
      constexpr auto dtype = dtype_short<Scalar>();
      assert(dtype.size() == HeaderLayout::DTYPE_SIZE);
      if (!m_out.write(dtype.data(), dtype.size())) {
        Igor::Warn("Could not write the data type to `{}`.", m_filename);
        return false;
      }
    }

    // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
    static_assert(sizeof(rows) == HeaderLayout::ROWS_SIZE);
    if (!m_out.write(reinterpret_cast<const char*>(&rows), sizeof(rows))) {
      Igor::Warn("Could not write number of rows to `{}`.", m_filename);
      return false;
    }
    static_assert(sizeof(cols) == HeaderLayout::COLS_SIZE);
    if (!m_out.write(reinterpret_cast<const char*>(&cols), sizeof(cols))) {
      Igor::Warn("Could not write number of cols to `{}`.", m_filename);
      return false;
    }
    static_assert(sizeof(is_row_major) == HeaderLayout::IS_ROW_MAJOR_SIZE);
    if (!m_out.write(reinterpret_cast<const char*>(&is_row_major), sizeof(is_row_major))) {
      Igor::Warn("Could not write storage order to `{}`.", m_filename);
      return false;
    }
    // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)

    return true;
  }

  // -----------------------------------------------------------------------------------------------
 public:
  [[nodiscard]] auto
  write_data(const Eigen::Matrix<Scalar, ROWS, COLS, OPTIONS>& mat) noexcept -> bool {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    if (!m_out.write(reinterpret_cast<const char*>(mat.data()),
                     static_cast<std::streamsize>(sizeof(Scalar)) * mat.rows() * mat.cols())) {
      Igor::Warn("Could not write data to file to `{}`.", m_filename);
      return false;
    }

    return true;
  }

  [[nodiscard]] auto write_data(const Scalar& s) noexcept -> bool
  requires(ROWS == 1 && COLS == 1)
  {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    if (!m_out.write(reinterpret_cast<const char*>(&s),
                     static_cast<std::streamsize>(sizeof(Scalar)))) {
      Igor::Warn("Could not write data to file to `{}`.", m_filename);
      return false;
    }

    return true;
  }
};

// -------------------------------------------------------------------------------------------------
class NoopWriter {
 public:
  [[nodiscard]] constexpr auto write_data(const auto& /*ignored*/) noexcept -> bool { return true; }
};

}  // namespace Zap::IO

#endif  // ZAP_IO_INC_MATRIX_WRITER_HPP_
