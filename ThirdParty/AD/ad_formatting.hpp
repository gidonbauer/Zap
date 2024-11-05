#ifndef AD_FORMATTING_HPP_
#define AD_FORMATTING_HPP_

template <typename T>
struct fmt::formatter<ad::internal::active_type<T, ad::internal::ts_data<T>>>
    : public fmt::formatter<T> {
  template <typename FormatContext>
  constexpr auto format(const ad::internal::active_type<T, ad::internal::ts_data<T>>& x,
                        FormatContext& ctx) const noexcept {
    ctx.advance_to(fmt::formatter<T>::format(ad::value(x), ctx));
    ctx.advance_to(fmt::format_to(ctx.out(), " ("));
    ctx.advance_to(fmt::formatter<T>::format(ad::derivative(x), ctx));
    ctx.advance_to(fmt::format_to(ctx.out(), ")"));
    return ctx.out();
  }
};

#endif  // AD_FORMATTING_HPP_
