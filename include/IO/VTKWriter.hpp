#ifndef ZAP_IO_BASED_VTK_WRITER_HPP_
#define ZAP_IO_BASED_VTK_WRITER_HPP_

#include <fstream>
#include <string_view>

#include "CellBased/Grid.hpp"

namespace Zap::IO {

template <VTKFormat FORMAT = VTKFormat::UNSTRUCTURED_GRID>
class VTKWriter {
  std::string_view m_output_prefix;
  std::size_t m_counter = 0;

 public:
  constexpr VTKWriter(std::string_view output_prefix) noexcept
      : m_output_prefix(output_prefix) {}

  // -----------------------------------------------------------------------------------------------
  template <typename Float, std::size_t DIM>
  [[nodiscard]] constexpr auto
  write_data(const CellBased::UniformGrid<Float, DIM>& grid) noexcept -> bool {
    const std::string filename = std::format("{}_{}.vtk", m_output_prefix, m_counter);
    m_counter += 1;
    return to_vtk(filename, grid);
  }

  // -----------------------------------------------------------------------------------------------
  template <typename Float, std::size_t DIM>
  [[nodiscard]] constexpr auto
  to_vtk(const std::string& filename,
         const CellBased::UniformGrid<Float, DIM>& grid) const noexcept -> bool {
    static_assert(std::is_floating_point_v<std::remove_cvref_t<Float>>);
    std::string datatype{};
    if constexpr (std::is_same_v<double, std::remove_cvref_t<Float>>) {
      datatype = "double";
    } else if constexpr (std::is_same_v<float, std::remove_cvref_t<Float>>) {
      datatype = "float";
    }

    std::ofstream out(filename);
    if (!out) {
      Igor::Warn("Could not open file `{}`: {}", filename, std::strerror(errno));
      return false;
    }

    if constexpr (FORMAT == VTKFormat::STRUCTURED_POINTS) {
      out << "# vtk DataFile Version 2.0\n";
      out << "Grid data for 2D Burgers equation\n";
      out << "ASCII\n";  // TODO: Can we use binary format?
      out << "DATASET STRUCTURED_GRID \n";
      out << "DIMENSIONS " << grid.m_nx << ' ' << grid.m_ny << " 1\n";
      out << "POINTS " << grid.m_nx * grid.m_ny << ' ' << datatype << '\n';
      for (const auto& cell : grid.m_cells) {
        const auto x_mid     = cell.x_min + static_cast<Float>(0.5) * cell.dx;
        const auto y_mid     = cell.y_min + static_cast<Float>(0.5) * cell.dy;
        constexpr auto z_mid = Float{0};
        out << x_mid << ' ' << y_mid << ' ' << z_mid << '\n';
      }
      out << '\n';

      out << "POINT_DATA " << grid.m_nx * grid.m_ny << '\n';
      for (int d = 0; d < static_cast<int>(DIM); ++d) {
        out << "SCALARS u" << d + 1 << ' ' << datatype << ' ' << 1 << '\n';
        out << "LOOKUP_TABLE default\n";
        for (const auto& cell : grid.m_cells) {
          out << (cell.value(d) > 1e-8 ? cell.value(d) : 0) << '\n';
        }
        out << '\n';
      }
    } else if constexpr (FORMAT == VTKFormat::STRUCTURED_GRID) {
      out << "# vtk DataFile Version 2.0\n";
      out << "Grid data for 2D Burgers equation\n";
      out << "ASCII\n";  // TODO: Can we use binary format?
      out << "DATASET STRUCTURED_GRID \n";
      out << "DIMENSIONS " << grid.m_nx + 1 << ' ' << grid.m_ny + 1 << " 1\n";
      out << "POINTS " << (grid.m_nx + 1) * (grid.m_ny + 1) << ' ' << datatype << '\n';

      assert(grid.m_nx < std::numeric_limits<int>::max());
      assert(grid.m_ny < std::numeric_limits<int>::max());
      for (int yi = -1; yi < static_cast<int>(grid.m_ny); ++yi) {
        for (int xi = -1; xi < static_cast<int>(grid.m_nx); ++xi) {
          if (xi == -1 && yi == -1) {
            const auto& cell = grid.m_cells[grid.to_vec_idx(0, 0)];
            out << cell.x_min << ' ' << cell.y_min << ' ' << 0 << '\n';
          } else if (xi == -1) {
            const auto& cell = grid.m_cells[grid.to_vec_idx(0, static_cast<std::size_t>(yi))];
            out << cell.x_min << ' ' << cell.y_min + cell.dy << ' ' << 0 << '\n';
          } else if (yi == -1) {
            const auto& cell = grid.m_cells[grid.to_vec_idx(static_cast<std::size_t>(xi), 0)];
            out << cell.x_min + cell.dx << ' ' << cell.y_min << ' ' << 0 << '\n';
          } else {
            const auto& cell = grid.m_cells[grid.to_vec_idx(static_cast<std::size_t>(xi),
                                                            static_cast<std::size_t>(yi))];
            out << cell.x_min + cell.dx << ' ' << cell.y_min + cell.dy << ' ' << 0 << '\n';
          }
        }
      }
      out << '\n';

      out << "CELL_DATA " << grid.m_nx * grid.m_ny << '\n';
      for (int d = 0; d < static_cast<int>(DIM); ++d) {
        out << "SCALARS u" << d + 1 << ' ' << datatype << ' ' << 1 << '\n';
        out << "LOOKUP_TABLE default\n";
        for (const auto& cell : grid.m_cells) {
          out << (cell.value(d) > 1e-8 ? cell.value(d) : 0) << '\n';
        }
        out << '\n';
      }
    } else if constexpr (FORMAT == VTKFormat::UNSTRUCTURED_GRID) {
      out << "# vtk DataFile Version 2.0\n";
      out << "Grid data for 2D Burgers equation\n";
      out << "ASCII\n";  // TODO: Can we use binary format?
      out << "DATASET UNSTRUCTURED_GRID \n";
      out << '\n';

      out << "POINTS " << (grid.m_nx + 1) * (grid.m_ny + 1) << ' ' << datatype << '\n';
      assert(grid.m_nx < std::numeric_limits<int>::max());
      assert(grid.m_ny < std::numeric_limits<int>::max());
      for (int yi = -1; yi < static_cast<int>(grid.m_ny); ++yi) {
        for (int xi = -1; xi < static_cast<int>(grid.m_nx); ++xi) {
          if (xi == -1 && yi == -1) {
            const auto& cell = grid.m_cells[grid.to_vec_idx(0, 0)];
            out << cell.x_min << ' ' << cell.y_min << ' ' << 0 << '\n';
          } else if (xi == -1) {
            const auto& cell = grid.m_cells[grid.to_vec_idx(0, static_cast<std::size_t>(yi))];
            out << cell.x_min << ' ' << cell.y_min + cell.dy << ' ' << 0 << '\n';
          } else if (yi == -1) {
            const auto& cell = grid.m_cells[grid.to_vec_idx(static_cast<std::size_t>(xi), 0)];
            out << cell.x_min + cell.dx << ' ' << cell.y_min << ' ' << 0 << '\n';
          } else {
            const auto& cell = grid.m_cells[grid.to_vec_idx(static_cast<std::size_t>(xi),
                                                            static_cast<std::size_t>(yi))];
            out << cell.x_min + cell.dx << ' ' << cell.y_min + cell.dy << ' ' << 0 << '\n';
          }
        }
      }
      out << '\n';

      out << "CELLS " << grid.m_nx * grid.m_ny << ' ' << grid.m_nx * grid.m_ny * 5 << '\n';
      for (std::size_t yi = 0; yi < grid.m_ny; ++yi) {
        for (std::size_t xi = 0; xi < grid.m_nx; ++xi) {
          constexpr auto num_corners = 4;
          const auto corners         = grid.to_corner_idx(xi, yi);
          out << num_corners << ' '           //
              << corners.top_left << ' '      //
              << corners.top_right << ' '     //
              << corners.bottom_right << ' '  //
              << corners.bottom_left << '\n';
        }
      }
      out << '\n';

      out << "CELL_TYPES " << grid.m_nx * grid.m_ny << '\n';
      for (std::size_t i = 0; i < grid.m_cells.size(); ++i) {
        constexpr auto vtk_pixel_cell_type = 9;
        out << vtk_pixel_cell_type << '\n';
      }
      out << '\n';

      out << "CELL_DATA " << grid.m_nx * grid.m_ny << '\n';
      for (int d = 0; d < static_cast<int>(DIM); ++d) {
        out << "SCALARS u" << d + 1 << ' ' << datatype << ' ' << 1 << '\n';
        out << "LOOKUP_TABLE default\n";
        for (std::size_t yi = 0; yi < grid.m_ny; ++yi) {
          for (std::size_t xi = 0; xi < grid.m_nx; ++xi) {
            const auto& cell = grid.m_cells[grid.to_vec_idx(xi, yi)];
            out << (cell.value(d) > 1e-8 ? cell.value(d) : 0) << '\n';
          }
        }
        out << '\n';
      }
    } else {
      Igor::Panic("Unknown VTKFormat with value: {}",
                  static_cast<std::underlying_type_t<VTKFormat>>(FORMAT));
    }

    return true;
  }
};

}  // namespace Zap::IO

#endif  // ZAP_CELL_BASED_VTK_WRITER_HPP_
