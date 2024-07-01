#ifndef ZAP_IO_INC_MATRIX_READER_HPP_
#define ZAP_IO_INC_MATRIX_READER_HPP_

#include "IO/IncMatrixWriter.hpp"

namespace Zap::IO {

template <typename Scalar>
class IncMatrixReader {
  std::string m_filename;
  std::ifstream m_in;
  int64_t m_rows{};
  int64_t m_cols{};
  int8_t m_is_row_major{};
  int64_t m_elem_idx = -1;
  int64_t m_num_elem{};
  std::vector<Scalar> m_element;

  // -----------------------------------------------------------------------------------------------
  IncMatrixReader(std::string filename, std::ifstream in)
      : m_filename(std::move(filename)),
        m_in(std::move(in)) {}

  // -----------------------------------------------------------------------------------------------
 public:
  IncMatrixReader(const std::string& filename)
      : IncMatrixReader(filename, std::ifstream(filename, std::ios::in | std::ios::binary)) {
    if (!m_in) {
      Igor::Warn("Could not open `{}` for reading: {}", m_filename, std::strerror(errno));
      throw std::runtime_error(
          std::format("Opening file `{}` failed: {}", m_filename, std::strerror(errno)));
    }

    if (!read_header()) {
      Igor::Warn("Could not read header from file `{}`", m_filename);
      throw std::runtime_error(std::format("Reading header from file `{}` failed.", m_filename));
    }

    assert(m_rows >= 0 && m_cols >= 0);
    m_element.resize(static_cast<size_t>(m_rows * m_cols));

    if (!read_num_elements()) {
      Igor::Warn("Could not read number of elements in file `{}`", m_filename);
      throw std::runtime_error(
          std::format("Reading number of elements in file `{}` failed.", m_filename));
    }

    // if (!read_next()) {
    //   Igor::Warn("Could not read next element in file `{}`", m_filename);
    //   throw std::runtime_error(
    //       std::format("Reading next element in file `{}` failed.", m_filename));
    // }
  }

  // -----------------------------------------------------------------------------------------------
 private:
  [[nodiscard]] auto read_header() noexcept -> bool {
    using namespace std::string_view_literals;

    // Read magic string
    {
      std::array<char, IncMatrixHeaderLayout::MAGIC_STRING.size()> magic_string{};
      if (!m_in.read(magic_string.data(), magic_string.size())) {
        Igor::Warn("Could not read magic string from `{}`: {}", m_filename, std::strerror(errno));
        return false;
      }
      if (std::string_view magic_string_sv{magic_string.data(), magic_string.size()};
          IncMatrixHeaderLayout::MAGIC_STRING != magic_string_sv) {
        Igor::Warn("Invalid magic string `{}` from `{}`, expected `{}`.",
                   magic_string_sv,
                   m_filename,
                   IncMatrixHeaderLayout::MAGIC_STRING);
        return false;
      }
    }

    // Read dtype
    {
      std::array<char, IncMatrixHeaderLayout::DTYPE_SIZE> dtype{};
      if (!m_in.read(dtype.data(), dtype.size())) {
        Igor::Warn("Could not read dtype from `{}`: {}", m_filename, std::strerror(errno));
        return false;
      }
      std::string_view dtype_sv{dtype.data(), dtype.size()};
      if (dtype_short<Scalar>() != dtype_sv) {
        Igor::Warn(
            "Incompatible dtype, got `{}` but expected `{}`", dtype_sv, dtype_short<Scalar>());
        return false;
      }
    }

    // Read rows
    {
      std::array<char, IncMatrixHeaderLayout::ROWS_SIZE> rows{};
      if (!m_in.read(rows.data(), rows.size())) {
        Igor::Warn("Could not read number of rows from `{}`: {}", m_filename, std::strerror(errno));
        return false;
      }
      std::memcpy(&m_rows, rows.data(), rows.size());
    }

    // Read cols
    {
      std::array<char, IncMatrixHeaderLayout::COLS_SIZE> cols{};
      if (!m_in.read(cols.data(), cols.size())) {
        Igor::Warn(
            "Could not read number of columns from `{}`: {}", m_filename, std::strerror(errno));
        return false;
      }
      std::memcpy(&m_cols, cols.data(), cols.size());
    }

    // Read is_row_major
    {
      char byte = 0;
      static_assert(IncMatrixHeaderLayout::IS_ROW_MAJOR_SIZE == sizeof(byte));
      if (!m_in.read(&byte, sizeof(byte))) {
        Igor::Warn("Could not read `is_row_major` from `{}`: {}", m_filename, std::strerror(errno));
        return false;
      }
      m_is_row_major = byte > 0 ? 1 : 0;
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto matrix_size_bytes() -> int64_t {
    return m_rows * m_cols * static_cast<int64_t>(sizeof(Scalar));
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto read_num_elements() -> bool {
    assert(m_in.tellg() == IncMatrixHeaderLayout::HEADER_SIZE);

    const auto begin = m_in.tellg();
    if (begin == -1) {
      Igor::Warn("Could not get position in steam.");
      return false;
    }

    if (!m_in.seekg(0, std::ios::end)) {
      Igor::Warn("Could not move input position indicator.");
    }

    const auto end = m_in.tellg();
    if (end == -1) {
      Igor::Warn("Could not get position in steam.");
      return false;
    }

    if (!m_in.seekg(begin)) {
      Igor::Warn("Could not move input position indicator.");
    }

    if ((end - begin) % matrix_size_bytes() != 0) {
      Igor::Warn("Data size ({} bytes) is not a multiple of the matrix size ({} bytes)",
                 end - begin,
                 matrix_size_bytes());
      return false;
    }
    m_num_elem = (end - begin) / matrix_size_bytes();

    return true;
  }

  // -----------------------------------------------------------------------------------------------
 public:
  template <bool WARN_END = true>
  [[nodiscard]] auto read_next() -> bool {
    m_elem_idx += 1;
    if (m_elem_idx >= m_num_elem) {
      if constexpr (WARN_END) {
        Igor::Warn("Reached end of file.");
      }
      return false;
    }

    if (!m_in.read(reinterpret_cast<char*>(m_element.data()),  // NOLINT
                   matrix_size_bytes())) {
      Igor::Warn("Could not read element at index {} from file `{}`", m_elem_idx, m_filename);
      return false;
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto operator()(int64_t row, int64_t col) noexcept -> Scalar& {
    assert(row >= 0 && row < m_rows);
    assert(col >= 0 && col < m_cols);
    if (m_elem_idx < 0) {
      Igor::Panic("Need to call `read_next` at least once before indexing into the matrix.");
    }

    if (m_is_row_major != 0) {
      return m_element[static_cast<size_t>(row * m_cols + col)];
    }
    return m_element[static_cast<size_t>(col * m_rows + row)];
  }

  [[nodiscard]] auto operator()(int64_t row, int64_t col) const noexcept -> const Scalar& {
    assert(row >= 0 && row < m_rows);
    assert(col >= 0 && col < m_cols);
    if (m_elem_idx < 0) {
      Igor::Panic("Need to call `read_next` at least once before indexing into the matrix.");
    }

    if (m_is_row_major != 0) {
      return m_element[static_cast<size_t>(row * m_cols + col)];
    }
    return m_element[static_cast<size_t>(col * m_rows + row)];
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto is_row_major() const noexcept -> bool { return m_is_row_major > 0; }
  [[nodiscard]] auto rows() const noexcept -> int64_t { return m_rows; }
  [[nodiscard]] auto cols() const noexcept -> int64_t { return m_cols; }
  [[nodiscard]] auto num_elem() const noexcept -> int64_t { return m_num_elem; }
  [[nodiscard]] auto data() noexcept -> std::vector<Scalar>& {
    if (m_elem_idx < 0) {
      Igor::Panic("Need to call `read_next` at least once before accessing the data.");
    }
    return m_element;
  }
  [[nodiscard]] auto data() const noexcept -> const std::vector<Scalar>& {
    if (m_elem_idx < 0) {
      Igor::Panic("Need to call `read_next` at least once before accessing the data.");
    }
    return m_element;
  }
};

}  // namespace Zap::IO

#endif  // ZAP_IO_INC_MATRIX_READER_HPP_
