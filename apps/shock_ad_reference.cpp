#include <filesystem>

#include "IO/ToNpy.hpp"
#include "ShockAD_1D/Solver.hpp"

using namespace Zap::ShockAD_1D;
using namespace Zap::IO;

#include "Igor/Logging.hpp"
#include "Igor/Macros.hpp"
#include "Igor/Timer.hpp"

#include "Common.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "shock_ad_reference/"

using PassiveFloat           = double;
using ActiveFloat            = ad::gt1s<PassiveFloat>::type;
constexpr PassiveFloat X_MIN = 0.0;
constexpr PassiveFloat X_MAX = 2.0;

// - Available command line options ----------------------------------------------------------------
struct Args {
  size_t nx          = 250;
  PassiveFloat tend  = 1.0;
  PassiveFloat eps   = 0.0;
  PassiveFloat C     = 25.0;  // 50.0;  // 1.0;
  PassiveFloat alpha = 1.0;   // 1.0;   // 0.5;
  bool run_benchmark = false;
};

// - Print usage -----------------------------------------------------------------------------------
template <Igor::detail::Level level>
void usage(std::string_view prog, std::ostream& out) noexcept {
  Args args{};

  out << Igor::detail::level_repr(level) << "Usage: " << prog
      << " [--nx nx] [--tend tend] [--eps eps] [--C C] [--alpha alpha] [--run-benchmark]\n";
  out << "\t--nx                 Number of cells in x-direction, default is " << args.nx << '\n';
  out << "\t--tend               Final time for simulation, default is " << args.tend << '\n';
  out << "\t--eps                Initial perturbation, default is " << args.eps << '\n';
  out << "\t--C                  First un-smoothing parameter, default is " << args.C << '\n';
  out << "\t--alpha              Second un-smoothing parameter, default is " << args.alpha << '\n';
  out << "\t--run-benchmark      Run the benchmark, default is " << std::boolalpha
      << args.run_benchmark << '\n';
}

// - Parse command line arguments ------------------------------------------------------------------
[[nodiscard]] auto parse_args(int argc, char** argv) noexcept -> std::optional<Args> {
  Args args{};

  const auto prog = next_arg(argc, argv);
  IGOR_ASSERT(prog.has_value(), "Could not get program name from command line arguments.");

  for (auto arg = next_arg(argc, argv); arg.has_value(); arg = next_arg(argc, argv)) {
    using namespace std::string_view_literals;
    if (arg == "--nx"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --nx.\n";
        return std::nullopt;
      }
      args.nx = parse_size_t(*arg);
    } else if (arg == "--tend"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --tend.\n";
        return std::nullopt;
      }
      args.tend = parse_double(*arg);
    } else if (arg == "--eps"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --eps.\n";
        return std::nullopt;
      }
      args.eps = parse_double(*arg);
    } else if (arg == "--C"sv || arg == "-C"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --C.\n";
        return std::nullopt;
      }
      args.C = parse_double(*arg);
    } else if (arg == "--alpha"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --alpha.\n";
        return std::nullopt;
      }
      args.alpha = parse_double(*arg);
    } else if (arg == "--run-benchmark"sv || arg == "--run-benchmarks"sv) {
      args.run_benchmark = true;
    } else if (arg == "-h"sv || arg == "--help"sv) {
      usage<Igor::detail::Level::INFO>(*prog, std::cout);
      return std::nullopt;
    } else {
      usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
      std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC) << "Unknown argument `"
                << *arg << "`.\n";
      return std::nullopt;
    }
  }
  return args;
}

// -------------------------------------------------------------------------------------------------
// Chi function: 1 if x in [lo, up], 0 otherwise
template <typename T>
[[nodiscard]] constexpr auto chi(const T& x, const T& lo, const T& up) noexcept -> T {
  return static_cast<T>(x >= lo && x <= up);
}

// -------------------------------------------------------------------------------------------------
template <typename AT>
[[nodiscard]] constexpr auto x_shock(PassiveFloat t, AT eps) noexcept -> AT {
  return std::sqrt(1 + (1 + eps) * t);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto xi_shock(PassiveFloat t) noexcept -> PassiveFloat {
  return t / (2 * std::sqrt(1 + t));
}

// -------------------------------------------------------------------------------------------------
template <typename AT>
[[nodiscard]] constexpr auto u_eps(PassiveFloat x, PassiveFloat t, AT eps) noexcept -> AT {
  return (1 + eps) * x / (1 + (1 + eps) * t) *
         chi(ad::value(x), PassiveFloat{0}, ad::value(x_shock(t, eps)));
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto v(PassiveFloat x, PassiveFloat t) noexcept -> PassiveFloat {
  return x / ((1 + t) * (1 + t)) * chi(x, PassiveFloat{0}, x_shock(t, PassiveFloat{0}));
}

// -------------------------------------------------------------------------------------------------
// Perturbed initial value for Burgers equation
//   -> u_eps = u + eps * v + shock offsets
//   -> u = x * chi(x, 0, 1)
//   -> v = x * chi(x, 0, 1)
//   -> no shock offsets
template <typename AT, typename PT>
[[nodiscard]] constexpr auto
u0_eps_func(const std::vector<PT>& x,  // x grid
            const AT& eps  // perturbation factor (we want to differentiate w.r.t. eps later)
            ) -> std::vector<AT> {
  std::vector<AT> u0(x.size());
  std::transform(std::cbegin(x), std::cend(x), std::begin(u0), [&](const PT& xi) -> AT {
    return u_eps(xi, PassiveFloat{0}, eps);
  });
  return u0;
}

// -------------------------------------------------------------------------------------------------
auto main(int argc, char** argv) -> int {
  const auto args = parse_args(argc, argv);
  if (!args.has_value()) { return 1; }

  {
    std::error_code ec;
    std::filesystem::create_directories(OUTPUT_DIR, ec);
    if (ec) {
      Igor::Warn("Could not create directory `{}`: {}", OUTPUT_DIR, ec.message());
      return 1;
    }
  }

  // - Parameters ---------------------
  Igor::Info("tend  = {}", args->tend);
  Igor::Info("eps   = {}", args->eps);
  Igor::Info("C     = {}", args->C);
  Igor::Info("alpha = {}", args->alpha);

  // Perturbation factor, we want to differentiate w.r.t. eps later
  const ActiveFloat eps = args->eps;
  ad::derivative(eps)   = 1.0;

  // Initial shock loaction at time t=0
  constexpr PassiveFloat x1_0 = 1.0;
  // - Parameters ---------------------

  if (!args->run_benchmark) {
    Igor::Info("nx    = {}", args->nx);

    const PassiveFloat dx = (X_MAX - X_MIN) / static_cast<PassiveFloat>(args->nx - 1UZ);
    std::vector<PassiveFloat> x(args->nx);
    std::generate(std::begin(x), std::end(x), [i = 0UZ, dx]() mutable {
      return X_MIN + static_cast<PassiveFloat>(i++) * dx;
    });

    // Perturbed initial value for Burgers equation at t=0
    const std::vector<ActiveFloat> u0 = u0_eps_func(x, eps);

    Igor::ScopeTimer timer("Lax Friedrichs solver");
    const auto [t, u, x1] = solve_1d_burgers(x, u0, x1_0, args->tend, args->C, args->alpha);

    Igor::Info("Solution at time {}:", t.back());
    Igor::Info("x1 = {}", ad::value(x1.back()));
    Igor::Info("xi1 = {}", ad::derivative(x1.back()));

    constexpr auto x_filename   = OUTPUT_DIR "x.npy";
    constexpr auto t_filename   = OUTPUT_DIR "t.npy";
    constexpr auto x1_filename  = OUTPUT_DIR "x1.npy";
    constexpr auto xi1_filename = OUTPUT_DIR "xi1.npy";
    constexpr auto u_filename   = OUTPUT_DIR "u.npy";
    constexpr auto v_filename   = OUTPUT_DIR "v.npy";

    if (!vector_to_npy(x_filename, x)) { return 1; }
    Igor::Info("Saved grid to {}", x_filename);

    if (!vector_to_npy(t_filename, t)) { return 1; }
    Igor::Info("Saved time steps to {}", t_filename);

    if (!vector_to_npy(x1_filename, xi1_filename, x1)) { return 1; }
    Igor::Info("Saved shock position to {}", x1_filename);
    Igor::Info("Saved change in shock position to {}", xi1_filename);

    if (!matrix_to_npy(u_filename, v_filename, u)) { return 1; }
    Igor::Info("Saved solution to {}", u_filename);
    Igor::Info("Saved derivative to {}", v_filename);
  } else {
    Igor::ScopeTimer timer("Running benchmark");
    constexpr auto output_file = "./shock_ad_results.txt";
    std::ofstream out(output_file);
    if (!out) {
      Igor::Warn("Could not open output file `{}`: {}", output_file, std::strerror(errno));
    }
    Igor::Info("Write results to `{}`.", output_file);

    constexpr size_t N_SIMPSON = 10'000UZ;
    for (size_t n = 10; n <= 400UZ; n += 10) {
      const PassiveFloat dx = (X_MAX - X_MIN) / static_cast<PassiveFloat>(n - 1UZ);
      std::vector<PassiveFloat> x(n);
      std::generate(std::begin(x), std::end(x), [i = 0UZ, dx]() mutable {
        return X_MIN + static_cast<PassiveFloat>(i++) * dx;
      });

      // Perturbed initial value for Burgers equation at t=0
      const std::vector<ActiveFloat> u0 = u0_eps_func(x, eps);

      const auto [t, u, x1] = solve_1d_burgers(x, u0, x1_0, args->tend, args->C, args->alpha);

      const auto x1_expect  = x_shock(args->tend, 0.0);
      const auto xi1_expect = xi_shock(args->tend);

      out << "nx                = " << n << '\n';

      out << "x_shock           = " << std::setprecision(8) << ad::value(x1.back()) << '\n';
      out << "expected_x_shock  = " << std::setprecision(8) << x1_expect << '\n';
      out << "abs_err_x_shock   = " << std::setprecision(8)
          << std::abs(ad::value(x1.back()) - x1_expect) << '\n';
      out << "rel_err_x_shock   = " << std::setprecision(8)
          << std::abs((x1.back() - x1_expect) / x1_expect) << '\n';
      out << "xi_shock          = " << std::setprecision(8) << ad::derivative(x1.back()) << '\n';
      out << "expected_xi_shock = " << std::setprecision(8) << xi1_expect << '\n';
      out << "abs_err_xi_shock  = " << std::setprecision(8)
          << std::abs(ad::derivative(x1.back()) - xi1_expect) << '\n';
      out << "rel_err_xi_shock  = " << std::setprecision(8)
          << std::abs((ad::derivative(x1.back()) - xi1_expect) / xi1_expect) << '\n';

      auto L1_func_u = [&](PassiveFloat x_eval) -> PassiveFloat {
        return std::abs(ad::value(piecewise_linear(x, u.back(), ActiveFloat{x_eval})) -
                        u_eps(x_eval, args->tend, PassiveFloat{0}));
      };
      const auto L1_err_u = simpsons_rule_1d(L1_func_u, X_MIN, X_MAX, N_SIMPSON);

      auto L1_func_v = [&](PassiveFloat x_eval) -> PassiveFloat {
        return std::abs(ad::derivative(piecewise_linear(x, u.back(), ActiveFloat{x_eval})) -
                        v(x_eval, args->tend));
      };

      // TODO: What is the correct(?) offset, fixed, based on dx, zero? -> Shock tracking solver has offset 0
      const auto offset   = 0.35;  // args->C * std::pow(dx, args->alpha);
      const auto L1_err_v = simpsons_rule_1d(L1_func_v, X_MIN, x1_expect - offset, N_SIMPSON / 2) +
                            simpsons_rule_1d(L1_func_v, x1_expect + offset, X_MAX, N_SIMPSON / 2);

      out << "L1_err_u          = " << std::setprecision(8) << L1_err_u << '\n';
      out << "L1_err_v          = " << std::setprecision(8) << L1_err_v << '\n';
      out << "------------------------------------------------------------\n";
    }
  }
}
