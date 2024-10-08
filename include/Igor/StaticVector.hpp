// Copyright 2024 Gidon Bauer <gidon.bauer@rwth-aachen.de>

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:

// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef IGOR_STATIC_VECTOR_HPP_
#define IGOR_STATIC_VECTOR_HPP_

#include <cstddef>
#include <initializer_list>

#include "./Logging.hpp"

namespace Igor {

template <typename Element, size_t CAPACITY>
class StaticVector {
  std::array<Element, CAPACITY> m_data{};
  size_t m_size = 0;

 public:
  using value_type             = Element;
  using size_type              = std::array<Element, CAPACITY>::size_type;
  using difference_type        = std::array<Element, CAPACITY>::difference_type;
  using reference              = value_type&;
  using const_reference        = const value_type&;
  using pointer                = std::array<Element, CAPACITY>::pointer;
  using const_pointer          = std::array<Element, CAPACITY>::const_pointer;
  using iterator               = std::array<Element, CAPACITY>::iterator;
  using const_iterator         = std::array<Element, CAPACITY>::const_iterator;
  using reverse_iterator       = std::array<Element, CAPACITY>::reverse_iterator;
  using const_reverse_iterator = std::array<Element, CAPACITY>::const_reverse_iterator;

  // -----------------------------------------------------------------------------------------------
  constexpr StaticVector() noexcept = default;

  constexpr StaticVector(size_t size) noexcept
      : m_size(size) {
    if (size > CAPACITY) {
      Igor::Panic("Size {} is greater than capacity {}.", size, CAPACITY);
    }
  }

  constexpr StaticVector(size_t size, Element value) noexcept
      : m_size(size) {
    if (size > CAPACITY) {
      Igor::Panic("Size {} is greater than capacity {}.", size, CAPACITY);
    }
    std::fill(begin(), end(), value);
  }

  constexpr StaticVector(std::initializer_list<Element> values) noexcept {
    if (values.size() > CAPACITY) {
      Igor::Panic("Size {} is greater than capacity {}.", values.size(), CAPACITY);
    }
    m_size = values.size();
    std::move(values.begin(), values.end(), m_data.begin());
  }

  template <size_t OTHER_CAPACITY>
  constexpr StaticVector(const StaticVector<Element, OTHER_CAPACITY>& other) noexcept
      : m_size(other.size()) {
    if constexpr (CAPACITY < OTHER_CAPACITY) {
      if (other.size() > CAPACITY) {
        Igor::Panic("Size {} is greater than capacity {}.", other.size(), CAPACITY);
      }
    }
    std::copy(other.begin(), other.end(), begin());
  }

  template <size_t OTHER_CAPACITY>
  constexpr StaticVector(StaticVector<Element, OTHER_CAPACITY>&& other) noexcept
      : m_size(other.size()) {
    if constexpr (CAPACITY < OTHER_CAPACITY) {
      if (other.size() > CAPACITY) {
        Igor::Panic("Size {} is greater than capacity {}.", other.size(), CAPACITY);
      }
    }
    std::move(other.begin(), other.end(), begin());
  }

  template <size_t OTHER_CAPACITY>
  constexpr auto
  operator=(const StaticVector<Element, OTHER_CAPACITY>& other) noexcept -> StaticVector& {
    if constexpr (CAPACITY < OTHER_CAPACITY) {
      if (other.size() > CAPACITY) {
        Igor::Panic("Size {} is greater than capacity {}.", other.size(), CAPACITY);
      }
    }
    m_size = other.size();
    std::copy(other.begin(), other.end(), begin());

    return *this;
  }

  template <size_t OTHER_CAPACITY>
  constexpr auto
  operator=(StaticVector<Element, OTHER_CAPACITY>&& other) noexcept -> StaticVector& {
    if constexpr (CAPACITY < OTHER_CAPACITY) {
      if (other.size() > CAPACITY) {
        Igor::Panic("Size {} is greater than capacity {}.", other.size(), CAPACITY);
      }
    }
    m_size = other.size();
    std::move(other.begin(), other.end(), begin());

    return *this;
  }

  constexpr ~StaticVector() noexcept = default;

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto operator[](size_t idx) noexcept -> reference {
    IGOR_ASSERT(idx < m_size, "Index {} is out of bounds for vector of size {}.", idx, m_size);
    return m_data[idx];
  }

  [[nodiscard]] constexpr auto operator[](size_t idx) const noexcept -> const_reference {
    IGOR_ASSERT(idx < m_size, "Index {} is out of bounds for vector of size {}.", idx, m_size);
    return m_data[idx];
  }

  [[nodiscard]] constexpr auto at(size_t idx) -> reference {
    if (idx >= m_size) {
      throw std::out_of_range(
          std::format("Index {} is out of bounds for vector of size {}.", idx, m_size));
    }
    return m_data[idx];
  }

  [[nodiscard]] constexpr auto at(size_t idx) const -> const_reference {
    if (idx >= m_size) {
      throw std::out_of_range(
          std::format("Index {} is out of bounds for vector of size {}.", idx, m_size));
    }
    return m_data[idx];
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto begin() noexcept -> iterator { return m_data.begin(); }
  [[nodiscard]] constexpr auto begin() const noexcept -> const_iterator { return m_data.cbegin(); }
  [[nodiscard]] constexpr auto cbegin() const noexcept -> const_iterator { return m_data.cbegin(); }

  [[nodiscard]] constexpr auto end() noexcept -> iterator {
    return std::next(m_data.begin(), static_cast<difference_type>(m_size));
  }
  [[nodiscard]] constexpr auto end() const noexcept -> const_iterator {
    return std::next(m_data.cbegin(), static_cast<difference_type>(m_size));
  }
  [[nodiscard]] constexpr auto cend() const noexcept -> const_iterator {
    return std::next(m_data.cbegin(), static_cast<difference_type>(m_size));
  }

  [[nodiscard]] constexpr auto rbegin() noexcept -> reverse_iterator {
    return std::next(m_data.rbegin(), static_cast<difference_type>(CAPACITY - m_size));
  }
  [[nodiscard]] constexpr auto rbegin() const noexcept -> const_reverse_iterator {
    return std::next(m_data.crbegin(), static_cast<difference_type>(CAPACITY - m_size));
  }
  [[nodiscard]] constexpr auto crbegin() const noexcept -> const_reverse_iterator {
    return std::next(m_data.crbegin(), static_cast<difference_type>(CAPACITY - m_size));
  }

  [[nodiscard]] constexpr auto rend() noexcept -> reverse_iterator {
    return std::next(rbegin(), static_cast<difference_type>(m_size));
  }
  [[nodiscard]] constexpr auto rend() const noexcept -> const_reverse_iterator {
    return std::next(crbegin(), static_cast<difference_type>(m_size));
  }
  [[nodiscard]] constexpr auto crend() const noexcept -> const_reverse_iterator {
    return std::next(crbegin(), static_cast<difference_type>(m_size));
  }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto front() noexcept -> reference {
    IGOR_ASSERT(m_size > 0, "StaticVector is empty.");
    return m_data[0];
  }
  [[nodiscard]] constexpr auto front() const noexcept -> const_reference {
    IGOR_ASSERT(m_size > 0, "StaticVector is empty.");
    return m_data[0];
  }
  [[nodiscard]] constexpr auto back() noexcept -> reference {
    IGOR_ASSERT(m_size > 0, "StaticVector is empty.");
    return m_data[m_size - 1];
  }
  [[nodiscard]] constexpr auto back() const noexcept -> const_reference {
    IGOR_ASSERT(m_size > 0, "StaticVector is empty.");
    return m_data[m_size - 1];
  }

  [[nodiscard]] constexpr auto data() noexcept -> pointer { return m_data.data(); }
  [[nodiscard]] constexpr auto data() const noexcept -> const_pointer { return m_data.data(); }

  // -----------------------------------------------------------------------------------------------
  [[nodiscard]] constexpr auto empty() const noexcept -> bool { return m_size == 0; }

  [[nodiscard]] constexpr auto size() const noexcept -> size_t { return m_size; }

  [[nodiscard]] constexpr auto max_size() const noexcept -> size_t { return CAPACITY; }

  [[nodiscard]] constexpr auto capacity() const noexcept -> size_t { return CAPACITY; }

  constexpr void reserve(size_t new_capacity) const noexcept {
    IGOR_ASSERT(new_capacity <= CAPACITY,
                "Requested capacity {} cannot fit into StaticVector with capacity {}",
                new_capacity,
                CAPACITY);
  }

  constexpr void shrink_to_fit() const noexcept {}

  // -----------------------------------------------------------------------------------------------
  constexpr void clear() noexcept { m_size = 0; }

  constexpr void push_back(const Element& e) noexcept {
    IGOR_ASSERT(m_size + 1 <= CAPACITY, "StaticVector is already at max. capacity {}", CAPACITY);
    m_data[m_size++] = e;
  }

  constexpr void push_back(Element&& e) noexcept {
    IGOR_ASSERT(m_size + 1 <= CAPACITY, "StaticVector is already at max. capacity {}", CAPACITY);
    m_data[m_size++] = e;
  }

  template <typename... Params>
  constexpr void emplace_back(Params&&... params) noexcept {
    IGOR_ASSERT(m_size + 1 <= CAPACITY, "StaticVector is already at max. capacity {}", CAPACITY);
    m_data[m_size++] = value_type{std::forward<Params>(params)...};
  }

  // -----------------------------------------------------------------------------------------------
  constexpr auto erase(iterator pos) noexcept -> iterator { return erase(pos, std::next(pos)); }
  constexpr auto erase(iterator first, iterator last) noexcept -> iterator {
    std::move(last, end(), first);

    const auto num_elems_removed = std::distance(first, last);
    IGOR_ASSERT(num_elems_removed >= 0,
                "first > last, invalid iterator pair; num_elems_removed={}",
                num_elems_removed);
    IGOR_ASSERT(
        m_size >= static_cast<size_t>(num_elems_removed),
        "Try to remove more elements from StaticVector than its capacity; num_elems_removed={}",
        num_elems_removed);
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
    IGOR_ASSERT(num_elems_removed >= 0,
                "first > last, invalid iterator pair; num_elems_removed={}",
                num_elems_removed);
    IGOR_ASSERT(
        m_size >= static_cast<size_t>(num_elems_removed),
        "Try to remove more elements from StaticVector than its capacity; num_elems_removed={}",
        num_elems_removed);
    m_size -= static_cast<size_t>(num_elems_removed);

    return first;
  }
};

}  // namespace Igor

#endif  // IGOR_STATIC_VECTOR_HPP_
