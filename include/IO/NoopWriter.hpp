#ifndef ZAP_IO_NOOP_WRITER_HPP_
#define ZAP_IO_NOOP_WRITER_HPP_

namespace Zap::IO {

class NoopWriter {
 public:
  [[nodiscard]] constexpr auto write_data(const auto& /*ignored*/) noexcept -> bool { return true; }
};

}  // namespace Zap::IO

#endif  // ZAP_IO_NOOP_WRITER_HPP_
