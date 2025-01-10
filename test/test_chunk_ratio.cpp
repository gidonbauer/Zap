#include <gtest/gtest.h>

#include "CellBased/Chunk.hpp"

using namespace Zap::CellBased;

auto objective_function(size_t NX, size_t NY, size_t CX, size_t CY) noexcept -> size_t {
  const size_t a = NX * CY > NY * CX ? NX * CY - NY * CX : NY * CX - NX * CY;
  return a;
}

template <size_t NT, size_t NX, size_t NY>
void test_optimality(size_t CX, size_t CY) noexcept {
  EXPECT_LE(2 * CX, NX);
  EXPECT_LE(2 * CY, NY);
  EXPECT_LE(CX * CY, NT);

  auto constraints = [&](size_t CX_test, size_t CY_test) {
    if (NX * NY > 4 * NT) {
      return 2 * CX_test <= NX && 2 * CY_test <= NY && CX_test * CY_test == NT;
    }
    return 2 * CX_test <= NX && 2 * CY_test <= NY && CX_test * CY_test <= NT;
  };

  const auto objective_value = objective_function(NX, NY, CX, CY);
  for (size_t CX_test = 1; CX_test <= NT; ++CX_test) {
    for (size_t CY_test = 1; CY_test <= NT; ++CY_test) {
      if (constraints(CX_test, CY_test)) {
        EXPECT_LE(objective_value, objective_function(NX, NY, CX_test, CY_test))
            << "CX=" << CX_test << " and CY=" << CY_test
            << " is a better fit for problem with NT=" << NT << ", NX=" << NX << ", and NY=" << NY
            << " than CX=" << CX << " and CY=" << CY;
      }
    }
  }
}

TEST(ChunkRatio, SquareGrid_SquareThreads) {
  constexpr size_t NT = 16;
  constexpr size_t NX = 100;
  constexpr size_t NY = 100;
  const auto [CX, CY] = optimal_chunk_ratio(NT, NX, NY);
  test_optimality<NT, NX, NY>(CX, CY);

  EXPECT_EQ(CX * CY, NT);
  EXPECT_EQ(NX * CY, NY * CX);
}

TEST(ChunkRatio, SquareGrid_RectThreads) {
  constexpr size_t NT = 12;
  constexpr size_t NX = 100;
  constexpr size_t NY = 100;
  const auto [CX, CY] = optimal_chunk_ratio(NT, NX, NY);

  test_optimality<NT, NX, NY>(CX, CY);

  EXPECT_EQ(CX * CY, NT);
}

TEST(ChunkRatio, ManyThreads) {
  constexpr size_t NT = 1 << 16;
  constexpr size_t NX = 100;
  constexpr size_t NY = 100;
  const auto [CX, CY] = optimal_chunk_ratio(NT, NX, NY);

  test_optimality<NT, NX, NY>(CX, CY);

  EXPECT_LT(CX * CY, NT);
}

TEST(ChunkRatio, SingleThread) {
  constexpr size_t NT = 1;
  constexpr size_t NX = 100;
  constexpr size_t NY = 100;
  const auto [CX, CY] = optimal_chunk_ratio(NT, NX, NY);

  test_optimality<NT, NX, NY>(CX, CY);
}

TEST(ChunkRatio, RectGrid) {
  constexpr size_t NT = 32;
  constexpr size_t NX = 2000;
  constexpr size_t NY = 100;
  const auto [CX, CY] = optimal_chunk_ratio(NT, NX, NY);

  test_optimality<NT, NX, NY>(CX, CY);
  EXPECT_EQ(CX * CY, NT);
}

TEST(ChunkRatio, PrimeNumberThreads) {
  constexpr size_t NT = 47;
  constexpr size_t NX = 1000;
  constexpr size_t NY = 1000;
  const auto [CX, CY] = optimal_chunk_ratio(NT, NX, NY);

  test_optimality<NT, NX, NY>(CX, CY);
  EXPECT_EQ(CX * CY, NT);
}
