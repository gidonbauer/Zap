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

#ifndef IGOR_MEMORY_TO_STRING_HPP_
#define IGOR_MEMORY_TO_STRING_HPP_

#include <cstdint>
#include <string>

namespace Igor {

[[nodiscard]] auto memory_to_string(uint64_t mem_in_bytes) noexcept -> std::string {
  using namespace std::string_literals;
  constexpr uint64_t step_factor = 1024;

  if (mem_in_bytes < step_factor) {
    return std::to_string(mem_in_bytes) + " B"s;
  }
  if (mem_in_bytes < step_factor * step_factor) {
    return std::to_string(mem_in_bytes / step_factor) + " kB"s;
  }
  if (mem_in_bytes < step_factor * step_factor * step_factor) {
    return std::to_string(mem_in_bytes / (step_factor * step_factor)) + " MB"s;
  }
  return std::to_string(mem_in_bytes / (step_factor * step_factor * step_factor)) + " GB"s;
}

}  // namespace Igor

#endif  // IGOR_MEMORY_TO_STRING_HPP_
