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

#ifndef IGOR_PROGRESS_BAR_HPP_
#define IGOR_PROGRESS_BAR_HPP_

#include <cassert>
#include <iomanip>
#include <iostream>

namespace Igor {

class ProgressBar {
  std::size_t m_max_progress;
  std::size_t m_length;
  std::size_t m_progress = 0UZ;

  enum : char {
    DONE_CHAR     = '#',
    NOT_DONE_CHAR = '.',
  };

 public:
  constexpr ProgressBar(std::size_t max_progress, std::size_t length) noexcept
      : m_max_progress(max_progress),
        m_length(length - 5UZ) {
    assert(length < std::numeric_limits<int>::max());
  }

  constexpr void update() noexcept {
    m_progress = std::min(m_progress + 1, m_max_progress);
    show();
  }

  void show() const noexcept {
    const auto done_length = (m_length * m_progress) / m_max_progress;
    const auto done_prct   = (100UZ * m_progress) / m_max_progress;
    std::cout << "\r[";
    std::cout << std::string(done_length, DONE_CHAR);
    std::cout << std::string(m_length - done_length, NOT_DONE_CHAR);
    std::cout << "] ";
    std::cout << std::setw(3) << done_prct << "%" << std::flush;
  }
};

}  // namespace Igor

#endif  // IGOR_PROGRESS_BAR_HPP_
