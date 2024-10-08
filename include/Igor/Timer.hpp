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

#ifndef IGOR_TIMER_HPP_
#define IGOR_TIMER_HPP_

#include <chrono>

#include "./Logging.hpp"
#include "./Macros.hpp"

namespace Igor {

class ScopeTimer {
  std::string m_scope_name;
  std::chrono::high_resolution_clock::time_point m_t_begin;

 public:
  [[nodiscard]] ScopeTimer(std::string scope_name = "Scope") noexcept
      : m_scope_name(std::move(scope_name)),
        m_t_begin(std::chrono::high_resolution_clock::now()) {}

  ScopeTimer(const ScopeTimer& other) noexcept                    = delete;
  ScopeTimer(ScopeTimer&& other) noexcept                         = delete;
  auto operator=(const ScopeTimer& other) noexcept -> ScopeTimer& = delete;
  auto operator=(ScopeTimer&& other) noexcept -> ScopeTimer&      = delete;
  ~ScopeTimer() noexcept {
    const auto t_duration =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - m_t_begin);
    detail::Time("{} took {}.", m_scope_name, t_duration);
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define IGOR_TIME_SCOPE(...)                                                                       \
  if constexpr (const auto IGOR_COMBINE(IGOR__SCOPE__TIMER__NAME__, __LINE__) =                    \
                    Igor::ScopeTimer{__VA_ARGS__};                                                 \
                true)

}  // namespace Igor

#endif  // IGOR_TIMER_HPP_
