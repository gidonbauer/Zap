#include "IO/IncMatrixWriter.hpp"
#include "Scalar/BoundaryConditions.hpp"
#include "Scalar/Solver.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR)

auto main(int argc, char** argv) -> int {
  if (argc < 3) {
    Igor::Warn("Usage: {} <nx> <ny>", *argv);
    return 1;
  }

  using Float          = double;
  constexpr auto x_min = static_cast<Float>(0.0);
  constexpr auto x_max = static_cast<Float>(10.0);
  constexpr auto y_min = static_cast<Float>(0.0);
  constexpr auto y_max = static_cast<Float>(10.0);

  auto parse_size_t = [](const char* cstr) -> size_t {
    char* end        = nullptr;
    const size_t val = std::strtoul(cstr, &end, 10);
    if (end != cstr + std::strlen(cstr)) {  // NOLINT
      Igor::Panic("String `{}` contains non-digits.", cstr);
    }
    if (val == 0UL) {
      Igor::Panic("Could not parse string `{}` to size_t.", cstr);
    }
    if (val == std::numeric_limits<unsigned long>::max()) {
      Igor::Panic("Could not parse string `{}` to size_t: {}", cstr, std::strerror(errno));
    }
    return val;
  };

  const auto nx = parse_size_t(argv[1]);  // NOLINT
  const auto ny = parse_size_t(argv[2]);  // NOLINT
  assert(nx < std::numeric_limits<int>::max());
  assert(ny < std::numeric_limits<int>::max());

  Igor::Info("nx = {}", nx);
  Igor::Info("ny = {}", ny);

  const auto dx = (x_max - x_min) / static_cast<Float>(nx);
  const auto dy = (y_max - y_min) / static_cast<Float>(ny);

  Zap::Vector<Float> x(nx);
  for (int i = 0; i < static_cast<int>(nx); ++i) {
    x(i) = x_min + dx * (static_cast<Float>(i) + static_cast<Float>(0.5));
  }
  Zap::Vector<Float> y(static_cast<int>(ny));
  for (int i = 0; i < static_cast<int>(ny); ++i) {
    y(i) = y_min + dy * (static_cast<Float>(i) + static_cast<Float>(0.5));
  }

  Zap::Matrix<Float> u0(static_cast<int>(ny), static_cast<int>(nx));
  for (int yi = 0; yi < static_cast<int>(ny); ++yi) {
    for (int xi = 0; xi < static_cast<int>(nx); ++xi) {
      // u0(yi, xi) = static_cast<Float>(x(xi) >= (x_min + x_max) / 2);
      // u0(yi, xi) = static_cast<Float>(std::pow(x(xi) - (x_min + x_max) / 2, 2) +
      //                                     std::pow(y(yi) - (y_min + y_max) / 2, 2) <=
      //                                 1.0 * 1.0);
      // u0(yi, xi) =
      //     (x(xi) + y(yi)) * static_cast<Float>((x(xi) * x(xi) + y(yi) * y(yi)) <=
      //                                          std::pow((x_min + x_max + y_min + y_max) / 4, 2));
      // u0(yi, xi) = x(xi) * static_cast<Float>(x(xi) <= (x_min + x_max) / 2);
      // u0(yi, xi) = y(yi) * static_cast<Float>(y(yi) <= (y_min + y_max) / 2);
      // u0(yi, xi) = 2 * static_cast<Float>(x(xi) <= (x_min + x_max) / 2);
      u0(yi, xi) = 2 * static_cast<Float>(x(xi) >= (x_min + x_max) / 2);
    }
  }

  constexpr auto tend = static_cast<Float>(4.0);

  auto boundary = [](Zap::Matrix<Float>& u_next,
                     const Zap::Matrix<Float>& u_curr,
                     Float dt,
                     Float dx,
                     Float dy,
                     auto numerical_flux) {
    // return Zap::Scalar::zero_flux_boundary(u_next, u_curr, dt, dx, dy, numerical_flux);
    // return Zap::Scalar::periodic_boundary(u_next, u_curr, dt, dx, dy, numerical_flux);
    return Zap::Scalar::equal_value_boundary(u_next, u_curr, dt, dx, dy, numerical_flux);
  };

  constexpr auto u_filename = OUTPUT_DIR "u.dat";
  Zap::IO::IncMatrixWriter u_writer(u_filename, u0);

  constexpr auto t_filename = OUTPUT_DIR "t.dat";
  Zap::IO::IncMatrixWriter<Float, 1, 1, 0> t_writer(t_filename, 1, 1, 0);

  const auto u_res = Zap::Scalar::solve_2d_burgers(x, y, u0, tend, boundary, u_writer, t_writer);

  Igor::Info("Wrote solution to `{}`", u_filename);
  Igor::Info("Wrote time steps to `{}`", t_filename);
}
