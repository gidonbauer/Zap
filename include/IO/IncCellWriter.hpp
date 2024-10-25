#ifndef ZAP_IO_INC_CELL_WRITER_HPP_
#define ZAP_IO_INC_CELL_WRITER_HPP_

#include <cassert>
#include <fstream>
#include <string>

#include "CellBased/Grid.hpp"
#include "IO/Common.hpp"

namespace Zap::IO {

namespace IncCellHeaderLayout {

constexpr std::string_view MAGIC_STRING = "%SHOCK";
constexpr size_t DTYPE_SIZE             = 3;
constexpr size_t DIM_SIZE               = 8;
constexpr size_t HEADER_SIZE_F32 =
    MAGIC_STRING.size() + DTYPE_SIZE + 4 + 4 + DIM_SIZE + 4 + 4 + DIM_SIZE + DIM_SIZE;
constexpr size_t HEADER_SIZE_F64 =
    MAGIC_STRING.size() + DTYPE_SIZE + 8 + 8 + DIM_SIZE + 8 + 8 + DIM_SIZE + DIM_SIZE;

enum struct CellType : char { Cartesian, Cut };

}  // namespace IncCellHeaderLayout

template <typename ActiveFloat, typename PassiveFloat, size_t DIM>
class IncCellWriter {
  std::string m_filename;
  std::ofstream m_out;

  // -----------------------------------------------------------------------------------------------
  IncCellWriter(const std::string& filename) noexcept
      : m_filename(filename),
        m_out(filename, std::ios::binary) {}

 public:
  // -----------------------------------------------------------------------------------------------
  IncCellWriter(const std::string& filename,
                const CellBased::UniformGrid<ActiveFloat, PassiveFloat, DIM>& grid)
      : IncCellWriter(filename) {
    if (!m_out) {
      Igor::Warn("Could not open file `{}` for writing: {}", filename, std::strerror(errno));
      throw std::runtime_error(
          std::format("Could not open file `{}` for writing: {}", filename, std::strerror(errno)));
    }

    if (!write_header(
            grid.x_min(), grid.x_max(), grid.nx(), grid.y_min(), grid.y_max(), grid.ny())) {
      Igor::Warn("Could not write header.");
      throw std::runtime_error("Could not write header.");
    }
  }

 private:
  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto write_header(PassiveFloat x_min,
                                  PassiveFloat x_max,
                                  size_t nx,
                                  PassiveFloat y_min,
                                  PassiveFloat y_max,
                                  size_t ny) noexcept -> bool {
    if (!m_out.write(IncCellHeaderLayout::MAGIC_STRING.data(),
                     IncCellHeaderLayout::MAGIC_STRING.size())) {
      Igor::Warn("Could not write magic string to `{}`.", m_filename);
      return false;
    }

    {
      constexpr auto dtype = dtype_short<PassiveFloat>();
      static_assert(dtype.size() == IncCellHeaderLayout::DTYPE_SIZE);
      if (!m_out.write(dtype.data(), dtype.size())) {
        Igor::Warn("Could not write the data type to `{}`.", m_filename);
        return false;
      }
    }

    // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
    if (!m_out.write(reinterpret_cast<const char*>(&x_min), sizeof(x_min))) {
      Igor::Warn("Could not write x_min to `{}`", m_filename);
      return false;
    }
    if (!m_out.write(reinterpret_cast<const char*>(&x_max), sizeof(x_max))) {
      Igor::Warn("Could not write x_max to `{}`", m_filename);
      return false;
    }
    static_assert(sizeof(nx) == IncCellHeaderLayout::DIM_SIZE);
    if (!m_out.write(reinterpret_cast<const char*>(&nx), sizeof(nx))) {
      Igor::Warn("Could not write nx to `{}`", m_filename);
      return false;
    }
    if (!m_out.write(reinterpret_cast<const char*>(&y_min), sizeof(y_min))) {
      Igor::Warn("Could not write y_min to `{}`", m_filename);
      return false;
    }
    if (!m_out.write(reinterpret_cast<const char*>(&y_max), sizeof(y_max))) {
      Igor::Warn("Could not write y_max to `{}`", m_filename);
      return false;
    }
    static_assert(sizeof(ny) == IncCellHeaderLayout::DIM_SIZE);
    if (!m_out.write(reinterpret_cast<const char*>(&ny), sizeof(ny))) {
      Igor::Warn("Could not write ny to `{}`", m_filename);
      return false;
    }
    static_assert(sizeof(DIM) == IncCellHeaderLayout::DIM_SIZE);
    if (size_t dim = DIM; !m_out.write(reinterpret_cast<const char*>(&dim), sizeof(DIM))) {
      Igor::Warn("Could not write DIM to `{}`", m_filename);
      return false;
    }
    // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)

    const std::streamoff stream_pos = m_out.tellp();
    assert(stream_pos > 0);
    if constexpr (std::is_same_v<float, std::remove_cvref_t<PassiveFloat>>) {
      if (!(static_cast<size_t>(stream_pos) == IncCellHeaderLayout::HEADER_SIZE_F32)) {
        Igor::Warn("Wrote header of size {}, but expected it to have size {}.",
                   stream_pos,
                   IncCellHeaderLayout::HEADER_SIZE_F32);
        return false;
      }
    } else if constexpr (std::is_same_v<double, std::remove_cvref_t<PassiveFloat>>) {
      if (!(static_cast<size_t>(stream_pos) == IncCellHeaderLayout::HEADER_SIZE_F64)) {
        Igor::Warn("Wrote header of size {}, but expected it to have size {}.",
                   stream_pos,
                   IncCellHeaderLayout::HEADER_SIZE_F64);
        return false;
      }
    } else {
      Igor::Panic("Incompatible passive float-type `{}`.", Igor::type_name<PassiveFloat>());
    }
    return true;
  }

 public:
  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto
  write_data(const CellBased::UniformGrid<ActiveFloat, PassiveFloat, DIM>& grid) noexcept -> bool {
    for (const auto& cell : grid) {
      // Extend of cell
      const PassiveFloat x_min = cell.template x_min<CellBased::SIM_C>();
      const PassiveFloat y_min = cell.template y_min<CellBased::SIM_C>();
      const PassiveFloat dx    = cell.template dx<CellBased::SIM_C>();
      const PassiveFloat dy    = cell.template dy<CellBased::SIM_C>();

      if (cell.is_cartesian()) {
        // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
        m_out.put(static_cast<char>(IncCellHeaderLayout::CellType::Cartesian));

        m_out.write(reinterpret_cast<const char*>(&x_min), sizeof(x_min));
        m_out.write(reinterpret_cast<const char*>(&dx), sizeof(dx));
        m_out.write(reinterpret_cast<const char*>(&y_min), sizeof(y_min));
        m_out.write(reinterpret_cast<const char*>(&dy), sizeof(dy));

        if constexpr (ad::mode<ActiveFloat>::is_ad_type) {
          // TODO: Should we save the derivative values as well?
          Eigen::Vector<PassiveFloat, DIM> passive_value;
          std::transform(std::cbegin(cell.get_cartesian().value),
                         std::cend(cell.get_cartesian().value),
                         std::begin(passive_value),
                         [](const ActiveFloat& v) { return ad::value(v); });

          m_out.write(reinterpret_cast<const char*>(passive_value.data()),
                      DIM * sizeof(PassiveFloat));
        } else {
          static_assert(std::is_floating_point_v<ActiveFloat>,
                        "We assume ActiveFloat is a floating point type.");
          m_out.write(reinterpret_cast<const char*>(cell.get_cartesian().value.data()),
                      DIM * sizeof(ActiveFloat));
        }
        // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
      } else if (cell.is_cut()) {
        // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
        m_out.put(static_cast<char>(IncCellHeaderLayout::CellType::Cut));

        m_out.write(reinterpret_cast<const char*>(&x_min), sizeof(x_min));
        m_out.write(reinterpret_cast<const char*>(&dx), sizeof(dx));
        m_out.write(reinterpret_cast<const char*>(&y_min), sizeof(y_min));
        m_out.write(reinterpret_cast<const char*>(&dy), sizeof(dy));

        // Cut
        const auto& cell_value = cell.get_cut();
        m_out.put(static_cast<char>(cell_value.type));

        const auto cut1 = cell.template cut1<CellBased::SIM_C>();
        const auto cut2 = cell.template cut2<CellBased::SIM_C>();
        m_out.write(reinterpret_cast<const char*>(&cut1.x), sizeof(cut1.x));
        m_out.write(reinterpret_cast<const char*>(&cut1.y), sizeof(cut1.y));
        m_out.write(reinterpret_cast<const char*>(&cut2.x), sizeof(cut2.x));
        m_out.write(reinterpret_cast<const char*>(&cut2.y), sizeof(cut2.y));

        // Value
        if constexpr (ad::mode<ActiveFloat>::is_ad_type) {
          // TODO: Should we save the derivative values as well?
          Eigen::Vector<PassiveFloat, DIM> passive_value;

          std::transform(std::cbegin(cell_value.left_value),
                         std::cend(cell_value.left_value),
                         std::begin(passive_value),
                         [](const ActiveFloat& v) { return ad::value(v); });
          m_out.write(reinterpret_cast<const char*>(passive_value.data()),
                      DIM * sizeof(PassiveFloat));

          std::transform(std::cbegin(cell_value.right_value),
                         std::cend(cell_value.right_value),
                         std::begin(passive_value),
                         [](const ActiveFloat& v) { return ad::value(v); });
          m_out.write(reinterpret_cast<const char*>(passive_value.data()),
                      DIM * sizeof(PassiveFloat));
        } else {
          static_assert(std::is_floating_point_v<ActiveFloat>,
                        "We assume ActiveFloat is a floating point type.");
          m_out.write(reinterpret_cast<const char*>(cell_value.left_value.data()),
                      DIM * sizeof(ActiveFloat));
          m_out.write(reinterpret_cast<const char*>(cell_value.right_value.data()),
                      DIM * sizeof(ActiveFloat));
        }
        // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
      } else {
        Igor::Panic("Unknown cell type with variant index {}", cell.value.index());
      }
    }
    if (!m_out) {
      Igor::Warn("Could not write data to `{}`.", m_filename);
      return false;
    }

    return true;
  }
};

}  // namespace Zap::IO

#endif  // ZAP_IO_INC_CELL_WRITER_HPP_
