#ifndef ZAP_IO_INC_CELL_READER_HPP_
#define ZAP_IO_INC_CELL_READER_HPP_

#include "IO/IncCellWriter.hpp"

namespace Zap::IO {

template <typename Float, size_t DIM>
class IncCellReader {
  struct ReducedCartesianCell {
    Float x_min;
    Float dx;
    Float y_min;
    Float dy;
    Eigen::Vector<Float, DIM> value;
  };

  std::string m_filename;
  std::ifstream m_in;

  Float m_x_min = std::numeric_limits<Float>::quiet_NaN();
  Float m_x_max = std::numeric_limits<Float>::quiet_NaN();
  Float m_y_min = std::numeric_limits<Float>::quiet_NaN();
  Float m_y_max = std::numeric_limits<Float>::quiet_NaN();
  size_t m_nx   = std::numeric_limits<size_t>::max();
  size_t m_ny   = std::numeric_limits<size_t>::max();

  std::vector<ReducedCartesianCell> m_cells{};
  bool m_read_once = false;

  // -----------------------------------------------------------------------------------------------
  IncCellReader(std::string filename, std::ifstream in) noexcept
      : m_filename(std::move(filename)),
        m_in(std::move(in)) {}

 public:
  // -----------------------------------------------------------------------------------------------
  IncCellReader(const std::string& filename)
      : IncCellReader(filename, std::ifstream(filename, std::ios::binary)) {
    if (!m_in) {
      Igor::Warn("Could not open file `{}` for writing: {}", filename, std::strerror(errno));
      throw std::runtime_error(
          std::format("Could not open file `{}` for writing: {}", filename, std::strerror(errno)));
    }

    if (!read_header()) {
      Igor::Warn("Could not read header.");
      throw std::runtime_error("Could not read header.");
    }
    m_cells.resize(m_nx * m_ny);
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto read_next_float(Float& f) noexcept -> bool {
    if constexpr (std::is_same_v<float, std::remove_cvref_t<Float>>) {
      static_assert(sizeof(Float) == 4, "Invalid floating point type, does not have correct size.");
    } else if constexpr (std::is_same_v<double, std::remove_cvref_t<Float>>) {
      static_assert(sizeof(Float) == 8, "Invalid floating point type, does not have correct size.");
    } else {
      static_assert(false, "Unknown floating point type.");
    }

    if (!m_in.read(reinterpret_cast<char*>(&f), sizeof(f))) {  // NOLINT
      Igor::Warn("Could not read the next Float from `{}`: {}", m_filename, std::strerror(errno));
      return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto read_next_size_t(size_t& s) noexcept -> bool {
    static_assert(sizeof(size_t) == IncCellHeaderLayout::DIM_SIZE,
                  "size_t has wrong size, expected 64 bit");

    if (!m_in.read(reinterpret_cast<char*>(&s), sizeof(s))) {  // NOLINT
      Igor::Warn("Could not read the next size_t from `{}`: {}", m_filename, std::strerror(errno));
      return false;
    }
    return true;
  }

 private:
  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto read_header() noexcept -> bool {
    // Read magic string
    {
      std::array<char, IncCellHeaderLayout::MAGIC_STRING.size()> magic_string{};
      if (!m_in.read(magic_string.data(), magic_string.size())) {
        Igor::Warn("Could not read magic string from `{}`: {}", m_filename, std::strerror(errno));
        return false;
      }
      if (std::string_view magic_string_sv{magic_string.data(), magic_string.size()};
          IncCellHeaderLayout::MAGIC_STRING != magic_string_sv) {
        Igor::Warn("Invalid magic string `{}` from `{}`, expected `{}`.",
                   magic_string_sv,
                   m_filename,
                   IncCellHeaderLayout::MAGIC_STRING);
        return false;
      }
    }

    {
      std::array<char, IncCellHeaderLayout::DTYPE_SIZE> dtype_arr{};
      if (!m_in.read(dtype_arr.data(), dtype_arr.size())) {
        Igor::Warn("Could not read dtype from `{}`: {}", m_filename, std::strerror(errno));
        return false;
      }
      std::string_view dtype_sv{dtype_arr.data(), dtype_arr.size()};
      if (dtype_short<Float>() != dtype_sv) {
        Igor::Warn(
            "Incompatible dtype, got `{}` but expected `{}`", dtype_sv, dtype_short<Float>());
        return false;
      }
    }

    if (!read_next_float(m_x_min)) {
      return false;
    }
    if (!read_next_float(m_x_max)) {
      return false;
    }
    if (!read_next_size_t(m_nx)) {
      return false;
    }

    if (!read_next_float(m_y_min)) {
      return false;
    }
    if (!read_next_float(m_y_max)) {
      return false;
    }
    if (!read_next_size_t(m_ny)) {
      return false;
    }

    {
      size_t dim{};
      if (!read_next_size_t(dim)) {
        return false;
      }
      if (dim != DIM) {
        Igor::Warn("Incompatible dimension {}, expected {}", dim, DIM);
        return false;
      }
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
 public:
  template <bool WARN_END = true>
  [[nodiscard]] auto read_next() -> bool {
    for (size_t i = 0; i < m_nx * m_ny; ++i) {
      char cell_type{};
      if (!m_in.get(cell_type)) {
        if constexpr (WARN_END) {
          Igor::Warn("Could not read cell type from `{}`.", m_filename);
          if (m_in.eof()) {
            Igor::Warn("Reached end of file.");
          }
        }
        return false;
      }

      if (static_cast<IncCellHeaderLayout::CellType>(cell_type) ==
          IncCellHeaderLayout::CellType::Cartesian) {
        auto& cell = m_cells.at(i);
        if (!read_next_float(cell.x_min)) {
          if constexpr (WARN_END) {
            Igor::Warn("Could not read x_min from `{}`", m_filename);
            if (m_in.eof()) {
              Igor::Warn("Reached end of file.");
            }
          }
          return false;
        }
        if (!read_next_float(cell.dx)) {
          if constexpr (WARN_END) {
            Igor::Warn("Could not read dx from `{}`", m_filename);
            if (m_in.eof()) {
              Igor::Warn("Reached end of file.");
            }
          }
          return false;
        }
        if (!read_next_float(cell.y_min)) {
          if constexpr (WARN_END) {
            Igor::Warn("Could not read y_min from `{}`", m_filename);
            if (m_in.eof()) {
              Igor::Warn("Reached end of file.");
            }
          }
          return false;
        }
        if (!read_next_float(cell.dy)) {
          if constexpr (WARN_END) {
            Igor::Warn("Could not read dy from `{}`", m_filename);
            if (m_in.eof()) {
              Igor::Warn("Reached end of file.");
            }
          }
          return false;
        }
        if (!m_in.read(reinterpret_cast<char*>(cell.value.data()),  // NOLINT
                       DIM * sizeof(Float))) {
          if constexpr (WARN_END) {
            Igor::Warn("Could not read value from `{}`", m_filename);
            if (m_in.eof()) {
              Igor::Warn("Reached end of file.");
            }
          }
          return false;
        }
      } else {
        Igor::Panic("Unknown cell type {}", static_cast<int>(cell_type));
        return false;
      }
    }

    m_read_once = true;
    return true;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto
  operator[](size_t idx) const noexcept -> const ReducedCartesianCell& {
    if (!m_read_once) {
      Igor::Panic("No data available, call `read_next`.");
    }

    assert(idx < m_nx * m_ny);
    return m_cells[idx];
  }

  [[nodiscard]] constexpr auto cells() const noexcept -> const std::vector<ReducedCartesianCell>& {
    return m_cells;
  }
  [[nodiscard]] constexpr auto x_min() const noexcept -> Float { return m_x_min; }
  [[nodiscard]] constexpr auto x_max() const noexcept -> Float { return m_x_max; }
  [[nodiscard]] constexpr auto y_min() const noexcept -> Float { return m_y_min; }
  [[nodiscard]] constexpr auto y_max() const noexcept -> Float { return m_y_max; }
  [[nodiscard]] constexpr auto nx() const noexcept -> size_t { return m_nx; }
  [[nodiscard]] constexpr auto ny() const noexcept -> size_t { return m_ny; }
};

}  // namespace Zap::IO

#endif  // ZAP_IO_INC_CELL_READER_HPP_
