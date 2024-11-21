#ifndef ZAP_IO_INC_CELL_WRITER_HPP_
#define ZAP_IO_INC_CELL_WRITER_HPP_

#include <cassert>
#include <fstream>
#include <optional>
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

enum struct CutType : char {
  BOTTOM_LEFT,
  BOTTOM_RIGHT,
  TOP_RIGHT,
  TOP_LEFT,
  MIDDLE_HORI,
  MIDDLE_VERT,
};

}  // namespace IncCellHeaderLayout

template <typename ActiveFloat, typename PassiveFloat>
class IncCellWriter {
  std::string m_primal_filename;
  std::ofstream m_primal_out;

  std::optional<std::string> m_tangent_filename = std::nullopt;
  std::optional<std::ofstream> m_tangent_out    = std::nullopt;

  // -----------------------------------------------------------------------------------------------
  IncCellWriter(const std::string& filename) noexcept
      : m_primal_filename(filename),
        m_primal_out(filename, std::ios::binary) {}

  // -----------------------------------------------------------------------------------------------
  IncCellWriter(const std::string& primal_filename, const std::string& tangent_filename) noexcept
      : m_primal_filename(primal_filename),
        m_primal_out(primal_filename, std::ios::binary),
        m_tangent_filename(tangent_filename),
        m_tangent_out(std::ofstream(tangent_filename, std::ios::binary)) {}

 public:
  // -----------------------------------------------------------------------------------------------
  IncCellWriter(const std::string& filename,
                const CellBased::UniformGrid<ActiveFloat, PassiveFloat>& grid)
      : IncCellWriter(filename) {
    if (!m_primal_out) {
      Igor::Warn("Could not open file `{}` for writing: {}", filename, std::strerror(errno));
      throw std::runtime_error(
          std::format("Could not open file `{}` for writing: {}", filename, std::strerror(errno)));
    }

    if (!write_header(m_primal_out,
                      m_primal_filename,
                      grid.x_min(),
                      grid.x_max(),
                      grid.nx(),
                      grid.y_min(),
                      grid.y_max(),
                      grid.ny())) {
      Igor::Warn("Could not write header.");
      throw std::runtime_error("Could not write header.");
    }
  }

  // -----------------------------------------------------------------------------------------------
  IncCellWriter(const std::string& primal_filename,
                const std::string& tangent_filename,
                const CellBased::UniformGrid<ActiveFloat, PassiveFloat>& grid)
      : IncCellWriter(primal_filename, tangent_filename) {
    if constexpr (!ad::mode<ActiveFloat>::is_ad_type) {
      Igor::Warn("Cannot write derivative data when no AD type is present, ActiveFloat is {}",
                 Igor::type_name<ActiveFloat>());
      throw std::runtime_error(
          fmt::format("Cannot write derivative data when no AD type is present, ActiveFloat is {}",
                      Igor::type_name<ActiveFloat>()));
    }

    // - Check output file for primal solution -----------------------------------------------------
    if (!m_primal_out) {
      Igor::Warn(
          "Could not open file `{}` for writing: {}", m_primal_filename, std::strerror(errno));
      throw std::runtime_error(fmt::format(
          "Could not open file `{}` for writing: {}", m_primal_filename, std::strerror(errno)));
    }

    if (!write_header(m_primal_out,
                      m_primal_filename,
                      grid.x_min(),
                      grid.x_max(),
                      grid.nx(),
                      grid.y_min(),
                      grid.y_max(),
                      grid.ny())) {
      Igor::Warn("Could not write header in {}.", m_primal_filename);
      throw std::runtime_error(fmt::format("Could not write header in {}.", m_primal_filename));
    }

    // - Check output file for derivative ----------------------------------------------------------
    IGOR_ASSERT(m_tangent_out.has_value() && m_tangent_filename.has_value(),
                "Both tangent out-stream and filename must have a value.");
    if (!m_tangent_out) {
      Igor::Warn(
          "Could not open file `{}` for writing: {}", *m_tangent_filename, std::strerror(errno));
      throw std::runtime_error(fmt::format(
          "Could not open file `{}` for writing: {}", *m_tangent_filename, std::strerror(errno)));
    }

    if (!write_header(*m_tangent_out,
                      *m_tangent_filename,
                      grid.x_min(),
                      grid.x_max(),
                      grid.nx(),
                      grid.y_min(),
                      grid.y_max(),
                      grid.ny())) {
      Igor::Warn("Could not write header in {}.", *m_tangent_filename);
      throw std::runtime_error(fmt::format("Could not write header in {}.", *m_tangent_filename));
    }
  }

 private:
  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto write_header(std::ofstream& out,
                                  const std::string& filename,
                                  PassiveFloat x_min,
                                  PassiveFloat x_max,
                                  size_t nx,
                                  PassiveFloat y_min,
                                  PassiveFloat y_max,
                                  size_t ny) noexcept -> bool {
    if (!out.write(IncCellHeaderLayout::MAGIC_STRING.data(),
                   IncCellHeaderLayout::MAGIC_STRING.size())) {
      Igor::Warn("Could not write magic string to `{}`.", filename);
      return false;
    }

    {
      constexpr auto dtype = dtype_short<PassiveFloat>();
      static_assert(dtype.size() == IncCellHeaderLayout::DTYPE_SIZE);
      if (!out.write(dtype.data(), dtype.size())) {
        Igor::Warn("Could not write the data type to `{}`.", filename);
        return false;
      }
    }

    // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
    if (!out.write(reinterpret_cast<const char*>(&x_min), sizeof(x_min))) {
      Igor::Warn("Could not write x_min to `{}`", filename);
      return false;
    }
    if (!out.write(reinterpret_cast<const char*>(&x_max), sizeof(x_max))) {
      Igor::Warn("Could not write x_max to `{}`", filename);
      return false;
    }
    static_assert(sizeof(nx) == IncCellHeaderLayout::DIM_SIZE);
    if (!out.write(reinterpret_cast<const char*>(&nx), sizeof(nx))) {
      Igor::Warn("Could not write nx to `{}`", filename);
      return false;
    }
    if (!out.write(reinterpret_cast<const char*>(&y_min), sizeof(y_min))) {
      Igor::Warn("Could not write y_min to `{}`", filename);
      return false;
    }
    if (!out.write(reinterpret_cast<const char*>(&y_max), sizeof(y_max))) {
      Igor::Warn("Could not write y_max to `{}`", filename);
      return false;
    }
    static_assert(sizeof(ny) == IncCellHeaderLayout::DIM_SIZE);
    if (!out.write(reinterpret_cast<const char*>(&ny), sizeof(ny))) {
      Igor::Warn("Could not write ny to `{}`", filename);
      return false;
    }
    constexpr size_t DIM = 1;  // DIM is included for legacy reasons, we set it to 1
    static_assert(sizeof(DIM) == IncCellHeaderLayout::DIM_SIZE);
    if (!out.write(reinterpret_cast<const char*>(&DIM), sizeof(DIM))) {
      Igor::Warn("Could not write DIM to `{}`", filename);
      return false;
    }
    // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)

    const std::streamoff stream_pos = out.tellp();
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

  template <typename Accessor>
  [[nodiscard]] auto write_data_impl(std::ofstream& out,
                                     const std::string& filename,
                                     const CellBased::UniformGrid<ActiveFloat, PassiveFloat>& grid,
                                     Accessor accessor) noexcept -> bool {
    for (const auto& cell : grid) {
      // Extend of cell
      const PassiveFloat x_min = cell.template x_min<CellBased::SIM_C>();
      const PassiveFloat y_min = cell.template y_min<CellBased::SIM_C>();
      const PassiveFloat dx    = cell.template dx<CellBased::SIM_C>();
      const PassiveFloat dy    = cell.template dy<CellBased::SIM_C>();

      if (cell.is_cartesian()) {
        // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
        out.put(static_cast<char>(IncCellHeaderLayout::CellType::Cartesian));

        out.write(reinterpret_cast<const char*>(&x_min), sizeof(x_min));
        out.write(reinterpret_cast<const char*>(&dx), sizeof(dx));
        out.write(reinterpret_cast<const char*>(&y_min), sizeof(y_min));
        out.write(reinterpret_cast<const char*>(&dy), sizeof(dy));

        const PassiveFloat passive_value = accessor(cell.get_cartesian().value);
        out.write(reinterpret_cast<const char*>(&passive_value), sizeof(PassiveFloat));
        // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
      } else if (cell.is_cut()) {
        // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
        out.put(static_cast<char>(IncCellHeaderLayout::CellType::Cut));

        out.write(reinterpret_cast<const char*>(&x_min), sizeof(x_min));
        out.write(reinterpret_cast<const char*>(&dx), sizeof(dx));
        out.write(reinterpret_cast<const char*>(&y_min), sizeof(y_min));
        out.write(reinterpret_cast<const char*>(&dy), sizeof(dy));

        // Cut
        const auto& cell_value = cell.get_cut();
        switch (cell_value.entry_loc | cell_value.exit_loc) {
          case CellBased::LEFT | CellBased::BOTTOM:
            out.put(static_cast<char>(IncCellHeaderLayout::CutType::BOTTOM_LEFT));
            break;
          case CellBased::BOTTOM | CellBased::RIGHT:
            out.put(static_cast<char>(IncCellHeaderLayout::CutType::BOTTOM_RIGHT));
            break;
          case CellBased::RIGHT | CellBased::TOP:
            out.put(static_cast<char>(IncCellHeaderLayout::CutType::TOP_RIGHT));
            break;
          case CellBased::TOP | CellBased::LEFT:
            out.put(static_cast<char>(IncCellHeaderLayout::CutType::TOP_LEFT));
            break;
          case CellBased::LEFT | CellBased::RIGHT:
            out.put(static_cast<char>(IncCellHeaderLayout::CutType::MIDDLE_HORI));
            break;
          case CellBased::BOTTOM | CellBased::TOP:
            out.put(static_cast<char>(IncCellHeaderLayout::CutType::MIDDLE_VERT));
            break;
          default:
            Igor::Panic("Invalid combination entry_loc = {} and exit_loc = {}",
                        cell_value.entry_loc,
                        cell_value.exit_loc);
            std::unreachable();
        }

        CellBased::SimCoord<ActiveFloat> cut1;
        CellBased::SimCoord<ActiveFloat> cut2;
        ActiveFloat left_value;
        ActiveFloat right_value;
        if (cell_value.entry_loc < cell_value.exit_loc) {
          cut1        = cell.template cut_entry<CellBased::SIM_C>();
          cut2        = cell.template cut_exit<CellBased::SIM_C>();
          left_value  = cell_value.left_value;
          right_value = cell_value.right_value;
        } else {
          cut1        = cell.template cut_exit<CellBased::SIM_C>();
          cut2        = cell.template cut_entry<CellBased::SIM_C>();
          left_value  = cell_value.right_value;
          right_value = cell_value.left_value;
        }

        // TODO: Save derivative of cut1 and cut2
        CellBased::SimCoord<PassiveFloat> passive_cut = {ad::value(cut1.x), ad::value(cut1.y)};
        out.write(reinterpret_cast<const char*>(&passive_cut.x), sizeof(passive_cut.x));
        out.write(reinterpret_cast<const char*>(&passive_cut.y), sizeof(passive_cut.y));

        passive_cut = {ad::value(cut2.x), ad::value(cut2.y)};
        out.write(reinterpret_cast<const char*>(&passive_cut.x), sizeof(passive_cut.x));
        out.write(reinterpret_cast<const char*>(&passive_cut.y), sizeof(passive_cut.y));

        // Value
        PassiveFloat passive_value;

        passive_value = accessor(left_value);
        out.write(reinterpret_cast<const char*>(&passive_value), sizeof(PassiveFloat));

        passive_value = accessor(right_value);
        out.write(reinterpret_cast<const char*>(&passive_value), sizeof(PassiveFloat));
        // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
      } else {
        Igor::Panic("Unknown cell type with variant index {}", cell.cell_type.index());
      }
    }
    if (!out) {
      Igor::Warn("Could not write data to `{}`.", filename);
      return false;
    }

    return true;
  }

 public:
  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] auto
  write_data(const CellBased::UniformGrid<ActiveFloat, PassiveFloat>& grid) noexcept -> bool {
    const auto success_primal = write_data_impl(
        m_primal_out, m_primal_filename, grid, [](const ActiveFloat& v) -> PassiveFloat {
          return ad::value(v);
        });

    if (m_tangent_out.has_value() && m_tangent_filename.has_value()) {
      const auto success_tangent = write_data_impl(
          *m_tangent_out, *m_tangent_filename, grid, [](const ActiveFloat& v) -> PassiveFloat {
            return ad::derivative(v);
          });

      return success_primal && success_tangent;
    }

    return success_primal;
  }
};  // namespace Zap::IO

}  // namespace Zap::IO

#endif  // ZAP_IO_INC_CELL_WRITER_HPP_
