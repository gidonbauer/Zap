#ifndef ZAP_CELL_BASED_CHUNK_HPP_
#define ZAP_CELL_BASED_CHUNK_HPP_

#include <cstddef>
#include <utility>
#include <vector>

#include <omp.h>

#include "Igor/Logging.hpp"

namespace Zap::CellBased {

// -------------------------------------------------------------------------------------------------
struct Chunk {
  size_t xi_min;
  size_t xi_max;
  size_t yi_min;
  size_t yi_max;
};

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto optimal_chunk_ratio(size_t nt_target, size_t nx, size_t ny) noexcept
    -> std::pair<size_t, size_t> {
  IGOR_ASSERT(nt_target > 0, "nt  must be greater than zero but is {}", nt_target);

  bool swap = true;
  if (ny > nx) {
    std::swap(nx, ny);
    swap = false;
  }

  auto calc_val = [nx, ny](const std::pair<size_t, size_t>& d) {
    // nx / ny = d2 / d1 <=> nx * d1 = ny * d2 <=> nx * d1 - ny * d2 = 0
    const auto a = nx * d.first;
    const auto b = ny * d.second;
    return a >= b ? a - b : b - a;
  };

  std::pair<size_t, size_t> opt_pair = {1, nt_target};
  auto opt_val                       = calc_val(opt_pair);

  const auto upper_bound =
      static_cast<size_t>(std::floor(std::sqrt(static_cast<double>(nt_target))));
  for (size_t lower = 2; lower <= upper_bound; ++lower) {
    if (nt_target % lower == 0) {
      const std::pair<size_t, size_t> p = {lower, nt_target / lower};
      const auto val                    = calc_val(p);
      if (val < opt_val) {
        opt_pair = p;
        opt_val  = val;
      }
    }
  }

  if (swap) { std::swap(opt_pair.first, opt_pair.second); }

  // Handle case where nx and ny are smaller than num_chunks_x and num_chunks_y
  opt_pair.first  = std::min(opt_pair.first, nx / 2UZ);
  opt_pair.second = std::min(opt_pair.second, ny / 2UZ);

  return opt_pair;
};

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto get_num_threads() -> size_t {
  int n = -1;
  // clang-format off
  #pragma omp parallel
  {
    #pragma omp single
    { n = omp_get_num_threads(); }
  }
  // clang-format on
  IGOR_ASSERT(n > 0, "Number of threads must be larger than zero but is {}", n);
  return static_cast<size_t>(n);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto
generate_chunks(size_t num_chunks_x, size_t num_chunks_y, size_t nx, size_t ny) noexcept
    -> std::vector<Chunk> {
  const auto chunk_width     = nx / num_chunks_x;
  const auto chunk_height    = ny / num_chunks_y;
  const auto chunk_missing_x = nx - chunk_width * num_chunks_x;
  const auto chunk_missing_y = ny - chunk_height * num_chunks_y;

  IGOR_ASSERT(
      chunk_width >= 2UZ, "Chunk width must be larger or equal to 2 but is {}", chunk_width);
  IGOR_ASSERT(
      chunk_height >= 2UZ, "Chunk height must be larger or equal to 2 but is {}", chunk_height);

  std::vector<Chunk> chunks(num_chunks_x * num_chunks_y);
  for (size_t cyi = 0; cyi < num_chunks_y; ++cyi) {
    for (size_t cxi = 0; cxi < num_chunks_x; ++cxi) {
      const auto i = cxi + cyi * num_chunks_x;
      chunks[i]    = Chunk{
             .xi_min = cxi * chunk_width,
             .xi_max = (cxi + 1) * chunk_width + ((cxi + 1 == num_chunks_x) ? chunk_missing_x : 0),
             .yi_min = cyi * chunk_height,
             .yi_max = (cyi + 1) * chunk_height + ((cyi + 1 == num_chunks_y) ? chunk_missing_y : 0),
      };
    }
  }

  return chunks;
}

// -------------------------------------------------------------------------------------------------
template <typename UPDATE>
constexpr void update_chunk(const Chunk& chunk, const UPDATE& do_update) {
  for (size_t yi = chunk.yi_min + 1; yi < chunk.yi_max - 1; ++yi) {
    for (size_t xi = chunk.xi_min + 1; xi < chunk.xi_max - 1; ++xi) {
      do_update(xi, yi);
    }
  }
}

// -------------------------------------------------------------------------------------------------
template <typename UPDATE>
constexpr void stich_chunk(const Chunk& chunk, const UPDATE& do_update) {
  for (size_t xi = chunk.xi_min; xi < chunk.xi_max; ++xi) {
    do_update(xi, chunk.yi_min);
    do_update(xi, chunk.yi_max - 1);
  }
  for (size_t yi = chunk.yi_min + 1; yi < chunk.yi_max - 1; ++yi) {
    do_update(chunk.xi_min, yi);
    do_update(chunk.xi_max - 1, yi);
  }
}

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_CHUNK_HPP_
