#ifndef ZAP_IO_FWD_HPP_
#define ZAP_IO_FWD_HPP_

#include <cstddef>
#include <cstdint>

namespace Zap::IO {

enum struct VTKFormat : std::uint8_t { STRUCTURED_POINTS, STRUCTURED_GRID, UNSTRUCTURED_GRID };

template <VTKFormat>
class VTKWriter;

template <typename, size_t>
class IncCellWriter;

template <typename, size_t>
class IncCellReader;

}  // namespace Zap::IO

#endif  // ZAP_IO_FWD_HPP_
