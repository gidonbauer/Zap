#ifndef AD_STD_COMPATIBILITY_HPP_
#define AD_STD_COMPATIBILITY_HPP_

namespace std {

using ad::internal::ceil;
using ad::internal::floor;

using ad::internal::max;
using ad::internal::min;

using ad::internal::abs;
using ad::internal::acos;
using ad::internal::acosh;
using ad::internal::asin;
using ad::internal::asinh;
using ad::internal::atan;
using ad::internal::atan2;
using ad::internal::atanh;
using ad::internal::cos;
using ad::internal::cosh;
using ad::internal::erf;
using ad::internal::erfc;
using ad::internal::exp;
using ad::internal::expm1;
using ad::internal::fabs;
// using ad::internal::frexp;
using ad::internal::hypot;
using ad::internal::isfinite;
using ad::internal::isinf;
using ad::internal::isnan;
// using ad::internal::isnormal;
// using ad::internal::ldexp;
using ad::internal::log;
using ad::internal::log10;
using ad::internal::log1p;
// using ad::internal::lround;
using ad::internal::pow;
// using ad::internal::round;
using ad::internal::sin;
using ad::internal::sinh;
using ad::internal::sqrt;
using ad::internal::tan;
using ad::internal::tanh;

template <class AD_TAPE_REAL, class DATA_HANDLER>
struct numeric_limits<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER>>
    : public numeric_limits<AD_TAPE_REAL> {
  using AD_T = ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER>;
  static constexpr auto min() -> AD_T { return numeric_limits<typename AD_T::VALUE_TYPE>::min(); }
  static constexpr auto max() -> AD_T { return numeric_limits<typename AD_T::VALUE_TYPE>::max(); }
  static constexpr auto epsilon() -> AD_T {
    return numeric_limits<typename AD_T::VALUE_TYPE>::epsilon();
  }
  static constexpr auto round_error() -> AD_T {
    return numeric_limits<typename AD_T::VALUE_TYPE>::round_error();
  }
  static constexpr auto infinity() -> AD_T {
    return numeric_limits<typename AD_T::VALUE_TYPE>::infinity();
  }
  static constexpr auto quiet_NaN() -> AD_T {
    return numeric_limits<typename AD_T::VALUE_TYPE>::quiet_NaN();
  }
  static constexpr auto signaling_NaN() -> AD_T {
    return numeric_limits<typename AD_T::VALUE_TYPE>::signaling_NaN();
  }
  static constexpr auto denorm_min() -> AD_T {
    return numeric_limits<typename AD_T::VALUE_TYPE>::denorm_min();
  }
  static constexpr auto lowest() -> AD_T {
    return numeric_limits<typename AD_T::VALUE_TYPE>::lowest();
  }
};

}  // namespace std

#endif  // AD_STD_COMPATIBILITY_HPP_
