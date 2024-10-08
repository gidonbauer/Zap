#ifndef IGOR_MATH_HPP_
#define IGOR_MATH_HPP_

#include <cmath>
#include <concepts>
#include <limits>

namespace Igor {

template <std::floating_point T>
constexpr auto constexpr_sqrt(T x) noexcept -> T {
  if consteval {
    constexpr auto sqrt_newton_raphson = [](const auto& local_sqrt_newton_raphson,
                                            const T& local_x,
                                            const T& curr,
                                            const T& prev) -> T {
      return curr == prev ? curr
                          : local_sqrt_newton_raphson(local_sqrt_newton_raphson,
                                                      local_x,
                                                      static_cast<T>(0.5) * (curr + local_x / curr),
                                                      curr);
    };

    return x >= static_cast<T>(0) && x < std::numeric_limits<T>::infinity()
               ? sqrt_newton_raphson(sqrt_newton_raphson, x, x, static_cast<T>(0))
               : std::numeric_limits<T>::quiet_NaN();
  } else {
    return std::sqrt(x);
  }
}

}  // namespace Igor

#endif  // IGOR_MATH_HPP_
