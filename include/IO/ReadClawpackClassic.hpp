#ifndef ZAP_IO_READ_CLAWPACK_CLASSIC_HPP_
#define ZAP_IO_READ_CLAWPACK_CLASSIC_HPP_

#include <algorithm>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "Igor/Logging.hpp"

namespace Zap::IO {

template <typename Float>
class ClawpackSolution {
  Float m_time  = std::numeric_limits<Float>::quiet_NaN();
  size_t m_nx   = 0;
  size_t m_ny   = 0;
  Float m_x_min = std::numeric_limits<Float>::quiet_NaN();
  Float m_y_min = std::numeric_limits<Float>::quiet_NaN();
  Float m_dx    = std::numeric_limits<Float>::quiet_NaN();
  Float m_dy    = std::numeric_limits<Float>::quiet_NaN();
  std::vector<Float> m_data{};

  // -----------------------------------------------------------------------------------------------
  auto read_field(std::ifstream& in,
                  const std::string& filename,
                  auto& fieldvalue,
                  std::string& fieldname) noexcept -> bool {
    if (!(in >> fieldvalue)) {
      Igor::Warn("Could not read field value from `{}`", filename);
      return false;
    }
    if (!(in >> fieldname)) {
      Igor::Warn("Could not read field name from `{}`", filename);
      return false;
    }
    return true;
  }

  // -----------------------------------------------------------------------------------------------
  auto read_t_file(const std::string& filename) noexcept -> bool {
    using namespace std::string_literals;

    std::ifstream t_in(filename);
    if (!t_in) {
      Igor::Warn("Could not open file `{}`: {}", filename, std::strerror(errno));
      return false;
    }
    std::string fieldname{};
    int value = -1;

    // - Read time -----------------------------------------
    if (!read_field(t_in, filename, m_time, fieldname)) {
      Igor::Warn("Could not read time from `{}`", filename);
      return false;
    }
    if (fieldname != "time"s) {
      Igor::Warn("Expected first field to be `time`, but is `{}`", fieldname);
      return false;
    }

    // - Read number of equations --------------------------
    if (!read_field(t_in, filename, value, fieldname)) {
      Igor::Warn("Could not read number equations from `{}`", filename);
      return false;
    }
    if (value != 1 || fieldname != "num_eqn"s) {
      Igor::Warn(
          "Expected `num_eqn` with value 1, but got field `{}` with value {}", fieldname, value);
      return false;
    }

    // - Read number of states -----------------------------
    if (!read_field(t_in, filename, value, fieldname)) {
      Igor::Warn("Could not read number of states from `{}`", filename);
      return false;
    }
    if (value != 1 || fieldname != "nstates"s) {
      Igor::Warn(
          "Expected `nstates` with value 1, but got field `{}` with value {}", fieldname, value);
      return false;
    }

    // - Read number of auxiliary variables ----------------
    if (!read_field(t_in, filename, value, fieldname)) {
      Igor::Warn("Could not read number of auxiliary variables from `{}`", filename);
      return false;
    }
    if (value != 0 || fieldname != "num_aux"s) {
      Igor::Warn(
          "Expected `num_aux` with value 0, but got field `{}` with value {}", fieldname, value);
      return false;
    }

    // - Read number of dimensions -------------------------
    if (!read_field(t_in, filename, value, fieldname)) {
      Igor::Warn("Could not read number of dimensions from `{}`", filename);
      return false;
    }
    if (value != 2 || fieldname != "num_dim"s) {
      Igor::Warn(
          "Expected `num_dim` with value 2, but got field `{}` with value {}", fieldname, value);
      return false;
    }

    // - Read number of ghost cells ------------------------
    if (!read_field(t_in, filename, value, fieldname)) {
      Igor::Warn("Could not read number of ghost cells from `{}`", filename);
      return false;
    }
    if (value != 0 || fieldname != "num_ghost"s) {
      Igor::Warn(
          "Expected `num_ghost` with value 0, but got field `{}` with value {}", fieldname, value);
      return false;
    }

    // - Read file format ----------------------------------
    std::string format{};
    if (!read_field(t_in, filename, format, fieldname)) {
      Igor::Warn("Could not read file-format from `{}`", filename);
      return false;
    }
    if (format != "ascii"s || fieldname != "file_format"s) {
      Igor::Warn("Expected `file_format` with value `ascii`, but got field `{}` with value `{}`",
                 fieldname,
                 format);
      return false;
    }

    return true;
  }

  // -----------------------------------------------------------------------------------------------
  auto read_q_file(const std::string& filename) noexcept -> bool {
    using namespace std::string_literals;

    std::ifstream q_in(filename);
    if (!q_in) {
      Igor::Warn("Could not open file `{}`: {}", filename, std::strerror(errno));
      return false;
    }
    std::string fieldname{};
    int value = -1;

    // - Read patch value ----------------------------------
    if (!read_field(q_in, filename, value, fieldname)) {
      Igor::Warn("Could not read patch_number from `{}`", filename);
      return false;
    }
    if (value != 1 || fieldname != "patch_number"s) {
      Igor::Warn(
          "Expected `patch_number` with value 1, but got `{}` with value {}", fieldname, value);
      return false;
    }

    // - Read AMR level ------------------------------------
    if (!read_field(q_in, filename, value, fieldname)) {
      Igor::Warn("Could not read AMR level from `{}`", filename);
      return false;
    }
    if (value != 1 || fieldname != "AMR_level"s) {
      Igor::Warn("Expected `AMR_level` with value 1, but got `{}` with value {}", fieldname, value);
      return false;
    }

    // - Read number of cells in x-direction ---------------
    if (!read_field(q_in, filename, m_nx, fieldname)) {
      Igor::Warn("Could not read number of cells in x-direction from `{}`", filename);
      return false;
    }
    if (fieldname != "mx"s) {
      Igor::Warn("Expected `mx`, but got `{}`", fieldname);
      return false;
    }

    // - Read number of cells in y-direction ---------------
    if (!read_field(q_in, filename, m_ny, fieldname)) {
      Igor::Warn("Could not read number of cells in y-direction from `{}`", filename);
      return false;
    }
    if (fieldname != "my"s) {
      Igor::Warn("Expected `my`, but got `{}`", fieldname);
      return false;
    }

    // - Read number x_min ---------------------------------
    if (!read_field(q_in, filename, m_x_min, fieldname)) {
      Igor::Warn("Could not read number of cells in x_min from `{}`", filename);
      return false;
    }
    if (fieldname != "xlow"s) {
      Igor::Warn("Expected `xlow`, but got `{}`", fieldname);
      return false;
    }

    // - Read number y_min ---------------------------------
    if (!read_field(q_in, filename, m_y_min, fieldname)) {
      Igor::Warn("Could not read number of cells in y_min from `{}`", filename);
      return false;
    }
    if (fieldname != "ylow"s) {
      Igor::Warn("Expected `ylow`, but got `{}`", fieldname);
      return false;
    }

    // - Read number dx ---------------------------------
    if (!read_field(q_in, filename, m_dx, fieldname)) {
      Igor::Warn("Could not read number of cells in dx from `{}`", filename);
      return false;
    }
    if (fieldname != "dx"s) {
      Igor::Warn("Expected `dx`, but got `{}`", fieldname);
      return false;
    }

    // - Read number dy ---------------------------------
    if (!read_field(q_in, filename, m_dy, fieldname)) {
      Igor::Warn("Could not read number of cells in dy from `{}`", filename);
      return false;
    }
    if (fieldname != "dy"s) {
      Igor::Warn("Expected `dy`, but got `{}`", fieldname);
      return false;
    }

    return read_data(filename, q_in);
  }

  // -------------------------------------------------------------------------------------------------
  [[nodiscard]] auto vec_idx(size_t xi, size_t yi) const noexcept -> size_t {
    IGOR_ASSERT(xi < m_nx && yi < m_ny,
                "Index ({}, {}) is out of range for grid with size ({}, {})",
                xi,
                yi,
                m_nx,
                m_ny);
    return xi + yi * m_nx;
  }

  // -------------------------------------------------------------------------------------------------
  auto read_data(const std::string& q_filename, std::ifstream& q_in) noexcept -> bool {
    m_data.resize(m_nx * m_ny);

    std::string value{};
    for (size_t yi = 0; yi < m_ny; ++yi) {
      for (size_t xi = 0; xi < m_nx; ++xi) {
        if (!(q_in >> value)) {
          Igor::Warn("Could not read data for index ({}, {}) from `{}`.", xi, yi, q_filename);
        }

        char* end{};
        errno                   = 0;
        m_data[vec_idx(xi, yi)] = std::strtod(value.c_str(), &end);
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        if (end != value.c_str() + value.size()) {
          Igor::Warn("Invalid value `{}`.", value);
          return false;
        } else if (errno == ERANGE) {
          // Igor::Warn("Value `{}` is out of range for double, use 0 instead.", value);
          m_data[vec_idx(xi, yi)] = 0;
        }
      }
    }

    return true;
  }

 public:
  // -----------------------------------------------------------------------------------------------
  ClawpackSolution(std::string q_filename, std::string t_filename) {
    if (!read_t_file(t_filename)) {
      Igor::Warn("Could not read t file.");
      throw std::runtime_error("Could not read t file.");
    }

    if (!read_q_file(q_filename)) {
      Igor::Warn("Could not read q file.");
      throw std::runtime_error("Could not read q file.");
    }
  }

  // -------------------------------------------------------------------------------------------------
  [[nodiscard]] auto x_min() const noexcept -> Float { return m_x_min; }
  [[nodiscard]] auto x_max() const noexcept -> Float {
    return m_x_min + m_dx * static_cast<Float>(m_nx);
  }
  [[nodiscard]] auto y_min() const noexcept -> Float { return m_y_min; }
  [[nodiscard]] auto y_max() const noexcept -> Float {
    return m_y_min + m_dy * static_cast<Float>(m_ny);
  }

  // -------------------------------------------------------------------------------------------------
  [[nodiscard]] auto operator()(size_t xi, size_t yi) const noexcept -> Float {
    return m_data[vec_idx(xi, yi)];
  }

  // -------------------------------------------------------------------------------------------------
  [[nodiscard]] auto operator()(Float x, Float y) const noexcept -> Float {
    IGOR_ASSERT(x_min() <= x && x <= x_max(),
                "x={} is out of range for grid with x-range [{}, {}]",
                x,
                x_min(),
                x_max());
    IGOR_ASSERT(y_min() <= y && y <= y_max(),
                "y={} is out of range for grid with y-range [{}, {}]",
                y,
                y_min(),
                y_max());
    const auto xi =
        std::clamp(static_cast<size_t>(std::round((x - x_min()) / m_dx)), 0UZ, m_nx - 1UZ);
    const auto yi =
        std::clamp(static_cast<size_t>(std::round((y - y_min()) / m_dy)), 0UZ, m_ny - 1UZ);
    return operator()(xi, yi);
  }
};

}  // namespace Zap::IO

#endif  // ZAP_IO_READ_CLAWPACK_CLASSIC_HPP_
