#include <filesystem>
#include <numbers>

#include <AD/ad.hpp>

#define ZAP_NO_CHUNK_INFO
#include "CellBased/PredictShock.hpp"
#include "CellBased/Solver.hpp"

#include "IO/NoopWriter.hpp"
#include "IO/ToNpy.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Macros.hpp"
#include "Igor/Timer.hpp"

#include "Common.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "predict_shock/"

using namespace Zap::CellBased;
using namespace Zap::IO;
using namespace std::string_view_literals;

// - Setup -----------------------------------------------------------------------------------------
using PassiveFloat           = double;
using ActiveFloat            = ad::gt1s<PassiveFloat>::type;
constexpr PassiveFloat X_MIN = 0.0;
constexpr PassiveFloat X_MAX = 5.0;
constexpr PassiveFloat Y_MIN = 0.0;
constexpr PassiveFloat Y_MAX = 5.0;

// - Available command line options ----------------------------------------------------------------
struct Args {
  size_t nx                      = 25;
  size_t ny                      = 25;
  PassiveFloat tend              = 1.0;
  PassiveFloat CFL_safety_factor = 0.25;
  PassiveFloat delta_eps         = 1e-8;
  bool first_order               = false;
  bool simple_dif_func           = false;
};

// - Print usage -----------------------------------------------------------------------------------
template <Igor::detail::Level level>
void usage(std::string_view prog, std::ostream& out) noexcept {
  Args args{};

  out << Igor::detail::level_repr(level) << "Usage: " << prog
      << " [--nx nx] [--ny ny] [--tend tend] [--CFL CFL-safety-factor] [--delta-eps delta_eps] "
         "[--first-order] [--simple-dif-func]\n";
  out << "\t--nx              Number of cells in x-direction, default is " << args.nx << '\n';
  out << "\t--ny              Number of cells in y-direction, default is " << args.ny << '\n';
  out << "\t--tend            Final time for simulation, default is " << args.tend << '\n';
  out << "\t--CFL             Safety factor for CFL condition, must be in (0, 1), default is "
      << args.CFL_safety_factor << '\n';
  out << "\t--delta-eps       Delta epsilon for finite differences, default is " << args.delta_eps
      << '\n';
  out << "\t--first-order     Use first order finite differences, default is " << std::boolalpha
      << args.first_order << '\n';
  out << "\t--simple-dif-func Use a simplified version of the dif-function to compare the curves, "
         "default is "
      << std::boolalpha << args.simple_dif_func << '\n';
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
    } else if (arg == "--ny"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --ny.\n";
        return std::nullopt;
      }
      args.ny = parse_size_t(*arg);
    } else if (arg == "--tend"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --tend.\n";
        return std::nullopt;
      }
      args.tend = parse_double(*arg);
    } else if (arg == "--CFL"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --CFL.\n";
        return std::nullopt;
      }
      args.CFL_safety_factor = parse_double(*arg);
    } else if (arg == "--delta-eps"sv || arg == "--eps"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --delta-eps.\n";
        return std::nullopt;
      }
      args.delta_eps = parse_double(*arg);
    } else if (arg == "--first-order"sv) {
      args.first_order = true;
    } else if (arg == "--simple-dif-func"sv || arg == "--simple"sv) {
      args.simple_dif_func = true;
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
[[nodiscard]] auto run(size_t nx,
                       size_t ny,
                       PassiveFloat tend,
                       PassiveFloat eps_value,
                       PassiveFloat CFL_safety_factor) noexcept
    -> std::optional<UniformGrid<ActiveFloat, PassiveFloat>> {
  // Igor::ScopeTimer timer(std::format("Solver with nx={}, ny={}, eps={}", nx, ny, eps_value));

  ActiveFloat eps     = eps_value;
  ad::derivative(eps) = 1;

  constexpr PassiveFloat r = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;
  auto u0                  = [&](ActiveFloat x, ActiveFloat y) -> ActiveFloat {
    return (1 + eps) * (std::pow(x - X_MIN, 2) + std::pow(y - Y_MIN, 2)) *
           static_cast<ActiveFloat>((std::pow(x - X_MIN, 2) + std::pow(y - Y_MIN, 2)) <=
                                    std::pow(r, 2));
  };

  UniformGrid<ActiveFloat, PassiveFloat> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.periodic_boundary();

  auto init_shock = [=]<typename T>(T t) -> SimCoord<T> {
    IGOR_ASSERT(t >= 0 && t <= 1, "t must be in [0, 1] but is {}", t);
    return {
        r * std::cos(std::numbers::pi_v<PassiveFloat> / 2 * t),
        r * std::sin(std::numbers::pi_v<PassiveFloat> / 2 * t),
    };
  };

  if (!grid.cut_curve(init_shock)) { return std::nullopt; }
  grid.fill_quad(u0);

  NoopWriter grid_writer;
  NoopWriter t_writer;
  return solve_2d_burgers<ExtendType::NEAREST>(
      grid, tend, grid_writer, t_writer, CFL_safety_factor);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto run_passive(size_t nx,
                               size_t ny,
                               PassiveFloat tend,
                               PassiveFloat eps_value,
                               PassiveFloat CFL_safety_factor) noexcept
    -> std::optional<UniformGrid<PassiveFloat, PassiveFloat>> {
  // Igor::ScopeTimer timer(std::format("Solver with nx={}, ny={}, eps={}", nx, ny, eps_value));

  constexpr PassiveFloat r = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;
  auto u0                  = [=](PassiveFloat x, PassiveFloat y) -> PassiveFloat {
    return (1 + eps_value) * (std::pow(x - X_MIN, 2) + std::pow(y - Y_MIN, 2)) *
           static_cast<PassiveFloat>((std::pow(x - X_MIN, 2) + std::pow(y - Y_MIN, 2)) <=
                                     std::pow(r, 2));
  };

  UniformGrid<PassiveFloat, PassiveFloat> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.periodic_boundary();

  auto init_shock = [=]<typename T>(T t) -> SimCoord<T> {
    IGOR_ASSERT(t >= 0 && t <= 1, "t must be in [0, 1] but is {}", t);
    return {
        r * std::cos(std::numbers::pi_v<PassiveFloat> / 2 * t),
        r * std::sin(std::numbers::pi_v<PassiveFloat> / 2 * t),
    };
  };

  if (!grid.cut_curve(init_shock)) { return std::nullopt; }
  grid.fill_quad(u0);

  NoopWriter grid_writer;
  NoopWriter t_writer;
  return solve_2d_burgers<ExtendType::NEAREST>(
      grid, tend, grid_writer, t_writer, CFL_safety_factor);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto fd_derivative(size_t nx,
                                 size_t ny,
                                 PassiveFloat tend,
                                 PassiveFloat CFL_safety_factor,
                                 bool first_order,
                                 size_t expected_num_cuts,
                                 PassiveFloat delta_eps) noexcept {
  const auto res_plus_eps = run_passive(nx, ny, tend, delta_eps, CFL_safety_factor);
  if (!res_plus_eps.has_value()) { Igor::Panic("Solver for eps={} failed.", delta_eps); }

  const auto res_minus_eps =
      run_passive(nx, ny, tend, first_order ? PassiveFloat{0} : -delta_eps, CFL_safety_factor);
  if (!res_minus_eps.has_value()) {
    Igor::Panic("Solver for eps={} failed.", first_order ? PassiveFloat{0} : -delta_eps);
  }

  const auto shock_plus_eps  = res_plus_eps->get_shock_curve();
  const auto shock_minus_eps = res_minus_eps->get_shock_curve();

  if (shock_plus_eps.size() != shock_minus_eps.size() ||
      shock_plus_eps.size() != expected_num_cuts) {
    Igor::Warn("Incompatible shock curve sizes: eps=0 => {}, eps={} => {}, eps={} => {}",
               expected_num_cuts,
               delta_eps,
               shock_plus_eps.size(),
               first_order ? PassiveFloat{0} : -delta_eps,
               shock_minus_eps.size());
  }

  std::vector<SimCoord<PassiveFloat>> der_shock_fd(expected_num_cuts);
  for (size_t i = 0; i < expected_num_cuts; ++i) {
    der_shock_fd[i].x =
        (shock_plus_eps[i].x - shock_minus_eps[i].x) / (first_order ? delta_eps : 2 * delta_eps);
    der_shock_fd[i].y =
        (shock_plus_eps[i].y - shock_minus_eps[i].y) / (first_order ? delta_eps : 2 * delta_eps);
  }

  return std::make_tuple(shock_plus_eps, shock_minus_eps, der_shock_fd);
}
// -------------------------------------------------------------------------------------------------
auto main(int argc, char** argv) -> int {
  const auto args = parse_args(argc, argv);
  if (!args.has_value()) { return 1; }
  assert(args->nx < std::numeric_limits<int>::max());
  assert(args->ny < std::numeric_limits<int>::max());

  {
    std::error_code ec;
    std::filesystem::create_directories(OUTPUT_DIR, ec);
    if (ec) {
      Igor::Warn("Could not create directory `{}`: {}", OUTPUT_DIR, ec.message());
      return 1;
    }
  }

  if (args->tend <= 0.0) {
    Igor::Warn("tend must be larger than 0, but is {}.", args->tend);
    return 1;
  }

  Igor::Info("nx              = {}", args->nx);
  Igor::Info("ny              = {}", args->ny);
  Igor::Info("tend            = {}", args->tend);
  Igor::Info("CFL             = {}", args->CFL_safety_factor);
  Igor::Info("delta-eps       = {}", args->delta_eps);
  Igor::Info("first-order     = {}", args->first_order);
  Igor::Info("simple-dif-func = {}", args->simple_dif_func);
  std::cout << "------------------------------------------------------------\n";

  const auto res_eps_zero =
      run(args->nx, args->ny, args->tend, PassiveFloat{0}, args->CFL_safety_factor);
  if (!res_eps_zero.has_value()) {
    Igor::Warn("Solver for eps=0 failed.");
    return 1;
  }
  const auto [true_shock_zero, ad_derivative_zero] = [&] {
    const auto active_shock = res_eps_zero->get_shock_curve();
    std::vector<SimCoord<PassiveFloat>> shock(active_shock.size());
    std::vector<SimCoord<PassiveFloat>> derivative(active_shock.size());

    for (size_t i = 0; i < active_shock.size(); ++i) {
      const auto& p = active_shock[i];
      shock[i]      = {ad::value(p.x), ad::value(p.y)};
      derivative[i] = {ad::derivative(p.x), ad::derivative(p.y)};
    }

    return std::make_pair(shock, derivative);
  }();

  const auto fd_res              = fd_derivative(args->nx,
                                    args->ny,
                                    args->tend,
                                    args->CFL_safety_factor,
                                    args->first_order,
                                    true_shock_zero.size(),
                                    args->delta_eps);
  const auto& fd_derivative_zero = std::get<2>(fd_res);

  // Save shock curve for eps=0
  {
    const auto filename = fmt::format("{}shock_{:.6f}.npy", OUTPUT_DIR, 0.0);
    if (!points_to_npy(filename, true_shock_zero)) {
      Igor::Warn("Could not save shock to `{}`.", filename);
      return 1;
    }
  }

  // Save AD derivative
  {
    const auto filename = fmt::format("{}ad_derivative.npy", OUTPUT_DIR);
    if (!points_to_npy(filename, ad_derivative_zero)) {
      Igor::Warn("Could not save shock to `{}`.", filename);
      return 1;
    }
  }

  // Save FD derivative
  {
    const auto filename = fmt::format("{}fd_derivative.npy", OUTPUT_DIR);
    if (!points_to_npy(filename, fd_derivative_zero)) {
      Igor::Warn("Could not save shock to `{}`.", filename);
      return 1;
    }
  }

  const auto ad_results_filename = fmt::format("{}shock_pred_errors_ad.txt", OUTPUT_DIR);
  std::ofstream ad_res_out(ad_results_filename);
  if (!ad_res_out) {
    Igor::Warn(
        "Could not open file `{}` for writing: {}", ad_results_filename, std::strerror(errno));
    return 1;
  }

  const auto fd_results_filename = fmt::format("{}shock_pred_errors_fd.txt", OUTPUT_DIR);
  std::ofstream fd_res_out(fd_results_filename);
  if (!fd_res_out) {
    Igor::Warn(
        "Could not open file `{}` for writing: {}", fd_results_filename, std::strerror(errno));
    return 1;
  }

  const auto sanity_check_filename = fmt::format("{}shock_zero_errors.txt", OUTPUT_DIR);
  std::ofstream sanity_out(sanity_check_filename);
  if (!sanity_out) {
    Igor::Warn(
        "Could not open file `{}` for writing: {}", sanity_check_filename, std::strerror(errno));
    return 1;
  }

  // clang-format off
  constexpr std::array epss = {
      1e-5, 2e-5, 5e-5,
      1e-4, 2e-4, 5e-4,
      1e-3, 2e-3, 5e-3,
      1e-2, 2e-2, 5e-2,
      1e-1, 2e-1, 5e-1,
  };
  // clang-format on
  for (auto eps : epss) {
    const auto res_eps = run_passive(args->nx, args->ny, args->tend, eps, args->CFL_safety_factor);
    if (!res_eps_zero.has_value()) {
      Igor::Warn("Solver for eps={} failed.", eps);
      return 1;
    }

    const auto pred_shock_eps_ad = predict_shock_front(true_shock_zero, ad_derivative_zero, eps);
    {
      const auto filename = fmt::format("{}pred_shock_ad_{:.6f}.npy", OUTPUT_DIR, eps);
      if (!points_to_npy(filename, pred_shock_eps_ad)) {
        Igor::Warn("Could not save shock to `{}`.", filename);
        return 1;
      }
    }

    const auto pred_shock_eps_fd = predict_shock_front(true_shock_zero, fd_derivative_zero, eps);
    {
      const auto filename = fmt::format("{}pred_shock_fd_{:.6f}.npy", OUTPUT_DIR, eps);
      if (!points_to_npy(filename, pred_shock_eps_fd)) {
        Igor::Warn("Could not save shock to `{}`.", filename);
        return 1;
      }
    }

    const auto true_shock_eps = res_eps->get_shock_curve();
    {
      const auto filename = fmt::format("{}shock_{:.6f}.npy", OUTPUT_DIR, eps);
      if (!points_to_npy(filename, true_shock_eps)) {
        Igor::Warn("Could not save shock to `{}`.", filename);
        return 1;
      }
    }

    constexpr auto NUM_SLICES = 10UZ;
    const auto dif_ad =
        compare_curves<64UZ>(true_shock_eps, pred_shock_eps_ad, NUM_SLICES, args->simple_dif_func);

    const auto dif_fd =
        compare_curves<64UZ>(true_shock_eps, pred_shock_eps_fd, NUM_SLICES, args->simple_dif_func);

    const auto dif_no_pred =
        compare_curves<64UZ>(true_shock_eps, true_shock_zero, NUM_SLICES, args->simple_dif_func);

    ad_res_out << eps << " = " << dif_ad << '\n';
    fd_res_out << eps << " = " << dif_fd << '\n';
    sanity_out << eps << " = " << dif_no_pred << '\n';

    Igor::Info("{:.5f} => {:.12f} (AD)", eps, dif_ad);
    Igor::Info("           {:.12f} (FD)", dif_fd);
    Igor::Info("           {:.12f} (No predict)", dif_no_pred);
  }
  std::cout << "------------------------------------------------------------\n";
  Igor::Info("Saved shock fronts in `{}`.", OUTPUT_DIR);
  Igor::Info("Saved results in `{}` and `{}`.", ad_results_filename, fd_results_filename);
  Igor::Info("Saved sanity check in `{}`.", sanity_check_filename);
}
