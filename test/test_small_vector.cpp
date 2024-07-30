#include <gtest/gtest.h>

#define ZAP_SMALL_VECTOR_CAPACITY 14
#include "CellBased/SmallVector.hpp"

TEST(SmallVector, Initialize) {
  {
    Zap::CellBased::SmallVector<int> vec{};
    EXPECT_EQ(vec.size(), 0);
  }

  {
    Zap::CellBased::SmallVector<int> vec(10);
    EXPECT_EQ(vec.size(), 10);
  }

  {
    Zap::CellBased::SmallVector<int> vec{1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_EQ(vec.size(), 8);
    for (size_t i = 0; i < vec.size(); ++i) {
      EXPECT_EQ(i + 1, vec[i]);
    }
  }

  {
    Zap::CellBased::SmallVector<std::array<float, 2>> vec{{1, 2}, {3, 4}, {5, 6}, {7, 8}};
    ASSERT_EQ(vec.size(), 4);
    for (size_t i = 0; i < vec.size(); ++i) {
      EXPECT_FLOAT_EQ(static_cast<float>(2 * i + 1), vec[i][0]);
      EXPECT_FLOAT_EQ(static_cast<float>(2 * i + 2), vec[i][1]);
    }
  }

  {
    EXPECT_DEATH(Zap::CellBased::SmallVector<int> vec(Zap::CellBased::SMALL_VECTOR_CAPACITY + 1);
                 , "");
  }
}

TEST(SmallVector, IndexOperator) {
  Zap::CellBased::SmallVector<int> vec{1, 2, 3, 4, 5, 6, 7, 8};
  ASSERT_EQ(vec.size(), 8);
  EXPECT_DEATH(const auto _ = vec[9], "");
}

TEST(SmallVector, RangeBasedLoops) {
  {
    const Zap::CellBased::SmallVector<int> vec{};
    for ([[maybe_unused]] const auto e : vec) {
      EXPECT_TRUE(false) << "The loop should never be executed.";
    }
  }

  {
    const Zap::CellBased::SmallVector<int> vec{1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_EQ(vec.size(), 8);
    int i = 1;
    for (const auto e : vec) {
      EXPECT_EQ(i, e);
      i += 1;
    }
  }

  {
    static_assert(Zap::CellBased::SMALL_VECTOR_CAPACITY == 14,
                  "This test assumes that the small vector capacity is 14.");
    const Zap::CellBased::SmallVector<int> vec{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
    ASSERT_EQ(vec.size(), 14);
    int i = 1;
    for (const auto e : vec) {
      EXPECT_EQ(i, e);
      i += 1;
    }
  }

  {
    Zap::CellBased::SmallVector<int> vec{1, 2, 3, 4, 5, 6, 7, 8};
    for (auto& e : vec) {
      e *= 2;
    }
    int i = 1;
    for (const auto e : vec) {
      EXPECT_EQ(2 * i, e);
      i += 1;
    }
  }
}

TEST(SmallVector, PushBack) {
  Zap::CellBased::SmallVector<size_t> vec;
  for (size_t i = 0; i < Zap::CellBased::SMALL_VECTOR_CAPACITY; ++i) {
    vec.push_back(i);
  }

  for (size_t i = 0; i < Zap::CellBased::SMALL_VECTOR_CAPACITY; ++i) {
    EXPECT_EQ(i, vec[i]);
  }

  EXPECT_DEATH(vec.push_back(42), "");
}

TEST(SmallVector, Erase) {
  {
    Zap::CellBased::SmallVector<size_t> vec{0, 1, 2, 3, 4, 5};
    ASSERT_EQ(vec.size(), 6);
    const auto it = vec.erase(std::next(std::begin(vec), 3));  // NOLINT
    EXPECT_EQ(*it, 4);

    ASSERT_EQ(vec.size(), 5);
    for (size_t i = 0; i < 3; ++i) {
      EXPECT_EQ(i, vec[i]);
    }
    for (size_t i = 3; i < 5; ++i) {
      EXPECT_EQ(i + 1, vec[i]);
    }
  }

  {
    Zap::CellBased::SmallVector<size_t> vec{0, 1, 2, 3, 4, 5};
    ASSERT_EQ(vec.size(), 6);
    const auto it =  // NOLINT
        vec.erase(std::next(std::begin(vec), 3), std::next(std::begin(vec), 3));
    EXPECT_EQ(*it, 3);

    ASSERT_EQ(vec.size(), 6);
    for (size_t i = 0; i < 6; ++i) {
      EXPECT_EQ(i, vec[i]);
    }
  }

  {
    Zap::CellBased::SmallVector<size_t> vec{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    const auto remove =  // NOLINT
        std::find_if(std::begin(vec), std::end(vec), [](size_t e) { return e > 5; });
    const auto it = vec.erase(remove, std::end(vec));  // NOLINT
    EXPECT_EQ(it, std::end(vec));

    ASSERT_EQ(vec.size(), 6);
    for (size_t i = 0; i < 6; ++i) {
      EXPECT_EQ(i, vec[i]);
    }
  }

  {
    Zap::CellBased::SmallVector<size_t> vec{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    const auto it =  // NOLINT
        vec.erase(std::next(std::begin(vec), 3), std::next(std::begin(vec), 6));
    EXPECT_EQ(*it, 6);

    ASSERT_EQ(vec.size(), 11);
    for (size_t i = 0; i < 3; ++i) {
      EXPECT_EQ(i, vec[i]);
    }

    for (size_t i = 3; i < 11; ++i) {
      EXPECT_EQ(i + 3, vec[i]);
    }
  }

  // -----------------------------------------------------------------------------------------------

  {
    Zap::CellBased::SmallVector<size_t> vec{0, 1, 2, 3, 4, 5};
    ASSERT_EQ(vec.size(), 6);
    const auto it = vec.erase(std::next(std::cbegin(vec), 3));  // NOLINT
    EXPECT_EQ(*it, 4);

    ASSERT_EQ(vec.size(), 5);
    for (size_t i = 0; i < 3; ++i) {
      EXPECT_EQ(i, vec[i]);
    }
    for (size_t i = 3; i < 5; ++i) {
      EXPECT_EQ(i + 1, vec[i]);
    }
  }

  {
    Zap::CellBased::SmallVector<size_t> vec{0, 1, 2, 3, 4, 5};
    ASSERT_EQ(vec.size(), 6);
    const auto it =  // NOLINT
        vec.erase(std::next(std::cbegin(vec), 3), std::next(std::cbegin(vec), 3));
    EXPECT_EQ(*it, 3);

    ASSERT_EQ(vec.size(), 6);
    for (size_t i = 0; i < 6; ++i) {
      EXPECT_EQ(i, vec[i]);
    }
  }

  {
    Zap::CellBased::SmallVector<size_t> vec{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    const auto remove =  // NOLINT
        std::find_if(std::cbegin(vec), std::cend(vec), [](size_t e) { return e > 5; });
    const auto it = vec.erase(remove, std::cend(vec));  // NOLINT
    EXPECT_EQ(it, std::cend(vec));

    ASSERT_EQ(vec.size(), 6);
    for (size_t i = 0; i < 6; ++i) {
      EXPECT_EQ(i, vec[i]);
    }
  }

  {
    Zap::CellBased::SmallVector<size_t> vec{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    const auto it =  // NOLINT
        vec.erase(std::next(std::cbegin(vec), 3), std::next(std::cbegin(vec), 6));
    EXPECT_EQ(*it, 6);

    ASSERT_EQ(vec.size(), 11);
    for (size_t i = 0; i < 3; ++i) {
      EXPECT_EQ(i, vec[i]);
    }

    for (size_t i = 3; i < 11; ++i) {
      EXPECT_EQ(i + 3, vec[i]);
    }
  }
}
