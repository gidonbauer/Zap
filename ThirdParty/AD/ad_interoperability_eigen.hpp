#ifndef AD_INTEROPERABILITY_EIGEN_HPP_
#define AD_INTEROPERABILITY_EIGEN_HPP_

#include <Eigen/Core>

template <class AD_TAPE_REAL, class DATA_HANDLER>
struct Eigen::NumTraits<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER>>
    : Eigen::GenericNumTraits<
          typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER>::VALUE_TYPE> {
  using Real       = ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER>;
  using NonInteger = Real;
  using Literal    = Real;
  using Nested     = Real;

  using value_type = typename ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER>::VALUE_TYPE;

  EIGEN_DEVICE_FUNC static inline auto epsilon() -> Real {
    return NumTraits<value_type>::epsilon();
  }
  EIGEN_DEVICE_FUNC static inline auto dummy_precision() -> Real {
    return NumTraits<value_type>::dummy_precision();
  }

  enum {
    IsComplex             = NumTraits<value_type>::IsComplex,
    IsInteger             = NumTraits<value_type>::IsInteger,
    ReadCost              = NumTraits<value_type>::ReadCost,
    AddCost               = NumTraits<value_type>::AddCost,
    MulCost               = NumTraits<value_type>::MulCost,
    IsSigned              = NumTraits<value_type>::IsSigned,
    RequireInitialization = 1
  };

  EIGEN_DEVICE_FUNC static inline auto lowest() -> Real {
    if (IsInteger) {
      return (numext::numeric_limits<value_type>::min)();
    } else {
      return -(numext::numeric_limits<value_type>::max)();
    }
  }
};

template <typename BinaryOp, typename REAL, typename DATA>
struct Eigen::ScalarBinaryOpTraits<ad::internal::active_type<REAL, DATA>, REAL, BinaryOp> {
  using ReturnType = ad::internal::active_type<REAL, DATA>;
};

template <typename BinaryOp, typename REAL, typename DATA>
struct Eigen::ScalarBinaryOpTraits<REAL, ad::internal::active_type<REAL, DATA>, BinaryOp> {
  using ReturnType = ad::internal::active_type<REAL, DATA>;
};

#endif  // AD_INTEROPERABILITY_EIGEN_HPP_
