#ifndef ZAP_CELL_BASED_SMALL_VECTOR_HPP_
#define ZAP_CELL_BASED_SMALL_VECTOR_HPP_

// TODO: Replace by Igor::StaticVector

#include <array>
#include <initializer_list>
#ifdef ZAP_SMALL_VECTOR_CAPACITY
#include <type_traits>
#endif  // ZAP_SMALL_VECTOR_CAPACITY

#include "Igor/Logging.hpp"

namespace Zap::CellBased {

#ifndef ZAP_SMALL_VECTOR_CAPACITY
constexpr size_t SMALL_VECTOR_CAPACITY = 16;
#else
static_assert(std::is_convertible_v<decltype(ZAP_SMALL_VECTOR_CAPACITY), size_t>,
              "ZAP_SMALL_VECTOR_CAPACITY must be a value that can be convertible to size_t.");
constexpr size_t SMALL_VECTOR_CAPACITY = ZAP_SMALL_VECTOR_CAPACITY;
#endif  // ZAP_SMALL_VECTOR_CAPACITY

template <typename Element>
class SmallVector {
  std::array<Element, SMALL_VECTOR_CAPACITY> m_data{};
  size_t m_size = 0;

 public:
  // -----------------------------------------------------------------------------------------------
  using iterator        = typename std::array<Element, SMALL_VECTOR_CAPACITY>::iterator;
  using const_iterator  = typename std::array<Element, SMALL_VECTOR_CAPACITY>::const_iterator;
  using difference_type = typename std::array<Element, SMALL_VECTOR_CAPACITY>::difference_type;

  // -----------------------------------------------------------------------------------------------
  constexpr SmallVector() noexcept = default;

  constexpr SmallVector(size_t size) noexcept {
    if (size > SMALL_VECTOR_CAPACITY) {
      Igor::Panic("Size {} is larger than maximal capacity {}. Consider increasing the capacity "
                  "with the `ZAP_SMALL_VECTOR_CAPACITY` macro.",
                  size,
                  SMALL_VECTOR_CAPACITY);
    }
    m_size = size;
  }

  constexpr SmallVector(std::initializer_list<Element> init) noexcept {
    if (init.size() > SMALL_VECTOR_CAPACITY) {
      Igor::Panic("Size {} is larger than maximal capacity {}. Consider increasing the capacity "
                  "with the `ZAP_SMALL_VECTOR_CAPACITY` macro.",
                  init.size(),
                  SMALL_VECTOR_CAPACITY);
    }
    m_size = init.size();
    std::move(std::begin(init), std::end(init), std::begin(m_data));
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto empty() const noexcept -> bool { return m_size == 0; }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto size() const noexcept -> size_t { return m_size; }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto operator[](size_t idx) noexcept -> Element& {
    assert(idx < m_size);
    return m_data[idx];
  }
  [[nodiscard]] constexpr auto operator[](size_t idx) const noexcept -> const Element& {
    assert(idx < m_size);
    return m_data[idx];
  }

  // -----------------------------------------------------------------------------------------------
  constexpr void clear() noexcept { m_size = 0; }

  // -----------------------------------------------------------------------------------------------
  constexpr void push_back(const Element& e) noexcept {
    if (m_size + 1 > SMALL_VECTOR_CAPACITY) {
      Igor::Panic("Try to push an element into a full SmallVector (size={}, capacity={}). Consider "
                  "increasing the capacity with the `ZAP_SMALL_VECTOR_CAPACITY` macro.",
                  m_size,
                  SMALL_VECTOR_CAPACITY);
    }
    m_data[m_size++] = e;
  }
  constexpr void push_back(Element&& e) noexcept {
    if (m_size + 1 > SMALL_VECTOR_CAPACITY) {
      Igor::Panic("Try to push an element into a full SmallVector (size={}, capacity={})",
                  m_size,
                  SMALL_VECTOR_CAPACITY);
    }
    m_data[m_size++] = e;
  }

  // -----------------------------------------------------------------------------------------------
  constexpr auto erase(iterator pos) noexcept -> iterator { return erase(pos, std::next(pos)); }
  constexpr auto erase(iterator first, iterator last) noexcept -> iterator {
    std::move(last, end(), first);

    const auto num_elems_removed = std::distance(first, last);
    assert(num_elems_removed >= 0 && m_size >= static_cast<size_t>(num_elems_removed));
    m_size -= static_cast<size_t>(num_elems_removed);

    return first;
  }

  constexpr auto erase(const_iterator pos) noexcept -> const_iterator {
    return erase(pos, std::next(pos));
  }
  constexpr auto erase(const_iterator first, const_iterator last) noexcept -> const_iterator {
    iterator non_const_first = std::next(begin(), std::distance(cbegin(), first));
    std::move(last, cend(), non_const_first);

    const auto num_elems_removed = std::distance(first, last);
    assert(num_elems_removed >= 0 && m_size >= static_cast<size_t>(num_elems_removed));
    m_size -= static_cast<size_t>(num_elems_removed);

    return first;
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto begin() noexcept { return m_data.begin(); }
  [[nodiscard]] constexpr auto begin() const noexcept { return m_data.begin(); }
  [[nodiscard]] constexpr auto cbegin() const noexcept { return m_data.cbegin(); }
  [[nodiscard]] constexpr auto end() noexcept {
    return std::next(m_data.begin(), static_cast<difference_type>(m_size));
  }
  [[nodiscard]] constexpr auto end() const noexcept {
    return std::next(m_data.begin(), static_cast<difference_type>(m_size));
  }
  [[nodiscard]] constexpr auto cend() const noexcept {
    return std::next(m_data.cbegin(), static_cast<difference_type>(m_size));
  }
};

}  // namespace Zap::CellBased

#endif  // ZAP_CELL_BASED_SMALL_VECTOR_HPP_
