#include <cmath>
#include <filesystem>

// #define ZAP_STATIC_CUT
#include "CellBased/EigenDecomp.hpp"
#include "CellBased/Grid.hpp"
#include "CellBased/Solver.hpp"

#include "MatBased/BoundaryConditions.hpp"
#include "MatBased/Solver.hpp"

#include "IO/IncCellWriter.hpp"
#include "IO/IncMatrixWriter.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Macros.hpp"

#include "Common.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "example_quarter_circle/"

// - Setup -----------------------------------------------------------------------------------------
using ActiveFloat            = double;
using PassiveFloat           = double;
constexpr size_t DIM         = 1;
constexpr PassiveFloat X_MIN = 0.0;
constexpr PassiveFloat X_MAX = 5.0;
constexpr PassiveFloat Y_MIN = 0.0;
constexpr PassiveFloat Y_MAX = 5.0;

// - Available command line options ----------------------------------------------------------------
struct Args {
  bool run_benchmark             = false;
  size_t nx                      = 25;
  size_t ny                      = 25;
  PassiveFloat tend              = 1.0;
  PassiveFloat CFL_safety_factor = 0.5;
};

// - Print usage -----------------------------------------------------------------------------------
template <Igor::detail::Level level>
void usage(std::string_view prog, std::ostream& out) noexcept {
  Args args{};

  out << Igor::detail::level_repr(level) << "Usage: " << prog
      << " [--run-benchmark] [--nx nx] [--ny ny] [--tend tend] [--CFL CFL-safety-factor]\n";
  out << "\t--run-benchmark   Run the solver for different grid sizes and compare result against "
         "high-resolution Godunov solver, default is "
      << std::boolalpha << args.run_benchmark << '\n';
  out << "\t--nx              Number of cells in x-direction, default is " << args.nx << '\n';
  out << "\t--ny              Number of cells in y-direction, default is " << args.ny << '\n';
  out << "\t--tend            Final time for simulation, default is " << args.tend << '\n';
  out << "\t--CFL             Safety factor for CFL condition, must be in (0, 1], default is "
      << args.CFL_safety_factor << '\n';
}

// - Parse command line arguments ------------------------------------------------------------------
[[nodiscard]] auto parse_args(int& argc, char**& argv) noexcept -> std::optional<Args> {
  Args args{};

  const auto prog = next_arg(argc, argv);
  IGOR_ASSERT(prog.has_value(), "Could not get program name from command line arguments.");

  for (auto arg = next_arg(argc, argv); arg.has_value(); arg = next_arg(argc, argv)) {
    using namespace std::string_view_literals;
    if (arg == "--run-benchmark"sv) {
      args.run_benchmark = true;
    } else if (arg == "--nx"sv) {
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
    } else if (arg == "--CFL"sv || arg == "cfl"sv) {
      arg = next_arg(argc, argv);
      if (!arg.has_value()) {
        usage<Igor::detail::Level::PANIC>(*prog, std::cerr);
        std::cerr << Igor::detail::level_repr(Igor::detail::Level::PANIC)
                  << "Did not provide number for option --CFL.\n";
        return std::nullopt;
      }
      args.CFL_safety_factor = parse_double(*arg);
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
[[nodiscard]] constexpr auto sqr(auto x) noexcept { return x * x; }

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto u0(PassiveFloat x, PassiveFloat y) noexcept -> ActiveFloat {
  static_assert(DIM == 1);
  return (sqr(x - X_MIN) + sqr(y - Y_MIN)) *
         static_cast<ActiveFloat>((sqr(x - X_MIN) + sqr(y - Y_MIN)) <=
                                  sqr((X_MIN + X_MAX + Y_MIN + Y_MAX) / 4));
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto init_shock(PassiveFloat t) noexcept -> Zap::CellBased::SimCoord<ActiveFloat> {
  const auto r = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;
  return {
      .x = r * std::cos(std::numbers::pi_v<PassiveFloat> / 2 * t),
      .y = r * std::sin(std::numbers::pi_v<PassiveFloat> / 2 * t),
  };
};

// -------------------------------------------------------------------------------------------------
auto run_cell_based(size_t nx,
                    size_t ny,
                    PassiveFloat tend,
                    PassiveFloat CFL_safety_factor) noexcept
    -> std::optional<Zap::CellBased::UniformGrid<ActiveFloat, PassiveFloat, DIM>> {
  Zap::CellBased::UniformGrid<ActiveFloat, PassiveFloat, DIM> grid(
      X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  // grid.same_value_boundary();
  grid.periodic_boundary();
  if (!grid.cut_curve(init_shock)) { return std::nullopt; }
  grid.fill_four_point(u0);

  auto solver = Zap::CellBased::make_solver<Zap::CellBased::ExtendType::NEAREST>(
      Zap::CellBased::SingleEq::A{}, Zap::CellBased::SingleEq::B{});

  const std::string u_filename =
      OUTPUT_DIR "u_cell_based_" + std::to_string(nx) + "x" + std::to_string(ny) + ".grid";
  Zap::IO::IncCellWriter u_writer(u_filename, grid);
  const std::string t_filename =
      OUTPUT_DIR "t_cell_based_" + std::to_string(nx) + "x" + std::to_string(ny) + ".mat";
  Zap::IO::IncMatrixWriter<PassiveFloat, 1, 1, 0> t_writer(t_filename, 1, 1, 0);

  return solver.solve(grid, tend, u_writer, t_writer, CFL_safety_factor);
}

// -------------------------------------------------------------------------------------------------
auto run_mat_based(size_t nx, size_t ny, PassiveFloat tend, PassiveFloat CFL_safety_factor) noexcept
    -> std::optional<std::tuple<Zap::MatBased::Vector<PassiveFloat>,
                                Zap::MatBased::Vector<PassiveFloat>,
                                Zap::MatBased::Matrix<PassiveFloat>>> {
  const auto dx = (X_MAX - X_MIN) / static_cast<PassiveFloat>(nx);
  const auto dy = (Y_MAX - Y_MIN) / static_cast<PassiveFloat>(ny);

  Zap::MatBased::Vector<PassiveFloat> x(nx);
  std::generate(std::begin(x), std::end(x), [i = 0, dx]() mutable {
    return X_MIN + dx * (static_cast<PassiveFloat>(i++) + static_cast<PassiveFloat>(0.5));
  });
  Zap::MatBased::Vector<PassiveFloat> y(static_cast<int>(ny));
  std::generate(std::begin(y), std::end(y), [i = 0, dy]() mutable {
    return Y_MIN + dy * (static_cast<PassiveFloat>(i++) + static_cast<PassiveFloat>(0.5));
  });

  Zap::MatBased::Matrix<PassiveFloat> u0_grid(static_cast<int>(ny), static_cast<int>(nx));
  for (int yi = 0; yi < static_cast<int>(ny); ++yi) {
    for (int xi = 0; xi < static_cast<int>(nx); ++xi) {
      // - Fill center -----------------
      // u0_grid(yi, xi) = u0(x(xi), y(yi));

      // - Fill four point -------------
      u0_grid(yi, xi) = (u0(x(xi) - dx / 4, y(yi) - dy / 4) + u0(x(xi) + dx / 4, y(yi) - dy / 4) +
                         u0(x(xi) - dx / 4, y(yi) + dy / 4) + u0(x(xi) + dx / 4, y(yi) + dy / 4)) /
                        4;
    }
  }

  auto boundary = [](Zap::MatBased::Matrix<PassiveFloat>& u_next,
                     const Zap::MatBased::Matrix<PassiveFloat>& u_curr,
                     PassiveFloat dt,
                     PassiveFloat local_dx,
                     PassiveFloat local_dy,
                     auto numerical_flux_x,
                     auto numerical_flux_y) {
    return Zap::MatBased::periodic_boundary(
        u_next, u_curr, dt, local_dx, local_dy, numerical_flux_x, numerical_flux_y);
    // return Zap::MatBased::equal_value_boundary(
    //     u_next, u_curr, dt, local_dx, local_dy, numerical_flux_x, numerical_flux_y);
  };

  const std::string u_filename =
      OUTPUT_DIR "u_mat_based_" + std::to_string(nx) + "x" + std::to_string(ny) + ".mat";
  Zap::IO::IncMatrixWriter u_writer(u_filename, u0_grid);

  const std::string t_filename =
      OUTPUT_DIR "t_mat_based_" + std::to_string(nx) + "x" + std::to_string(ny) + ".mat";
  Zap::IO::IncMatrixWriter<PassiveFloat, 1, 1, 0> t_writer(t_filename, 1, 1, 0);

  const auto res = Zap::MatBased::solve_2d_burgers(
      x, y, u0_grid, tend, boundary, u_writer, t_writer, CFL_safety_factor);

  if (!res.has_value()) { return std::nullopt; }
  return std::make_tuple(x, y, *res);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto min_max_cell_value(
    const std::vector<Zap::CellBased::Cell<ActiveFloat, PassiveFloat, DIM>>& cells) noexcept
    -> std::pair<Eigen::Vector<ActiveFloat, DIM>, Eigen::Vector<ActiveFloat, DIM>> {
  Eigen::Vector<ActiveFloat, DIM> min = [cell = cells.front()] {
    assert(cell.is_cartesian() || cell.is_cut());
    if (cell.is_cartesian()) {
      return cell.get_cartesian().value;
    } else {
      Eigen::Vector<ActiveFloat, DIM> res;
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        res(i) = std::min(cell.get_cut().left_value(i), cell.get_cut().right_value(i));
      }
      return res;
    }
  }();

  Eigen::Vector<ActiveFloat, DIM> max = [cell = cells.front()] {
    assert(cell.is_cartesian() || cell.is_cut());
    if (cell.is_cartesian()) {
      return cell.get_cartesian().value;
    } else {
      Eigen::Vector<ActiveFloat, DIM> res;
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        res(i) = std::min(cell.get_cut().left_value(i), cell.get_cut().right_value(i));
      }
      return res;
    }
  }();

  for (const auto& cell : cells) {
    if (cell.is_cartesian()) {
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        min(i) = std::min(min(i), cell.get_cartesian().value(i));
        max(i) = std::max(max(i), cell.get_cartesian().value(i));
      }
    } else {
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        min(i) = std::min({min(i), cell.get_cut().left_value(i), cell.get_cut().right_value(i)});
        max(i) = std::max({max(i), cell.get_cut().left_value(i), cell.get_cut().right_value(i)});
      }
    }
  }

  return std::make_pair(min, max);
}

// -------------------------------------------------------------------------------------------------
constexpr auto eval_mat_at(PassiveFloat x,
                           PassiveFloat y,
                           const Zap::MatBased::Vector<PassiveFloat>& xs,
                           const Zap::MatBased::Vector<PassiveFloat>& ys,
                           const Zap::MatBased::Matrix<PassiveFloat>& u) noexcept -> PassiveFloat {
  assert(ys.size() == u.rows());
  assert(xs.size() == u.cols());
  assert(xs.size() >= 2);
  assert(ys.size() >= 2);

  if (std::isnan(x) || std::isnan(y)) { return std::numeric_limits<PassiveFloat>::quiet_NaN(); }

  const auto dx = xs(1) - xs(0);
  const auto x_idx =
      static_cast<Eigen::Index>(std::clamp(std::round((x - xs(0)) / dx),
                                           static_cast<PassiveFloat>(0),
                                           static_cast<PassiveFloat>(xs.rows() - 1)));

  const auto dy = ys(1) - ys(0);
  const auto y_idx =
      static_cast<Eigen::Index>(std::clamp(std::round((y - ys(0)) / dy),
                                           static_cast<PassiveFloat>(0),
                                           static_cast<PassiveFloat>(ys.rows() - 1)));

  return u(y_idx, x_idx);
}

// -------------------------------------------------------------------------------------------------
template <typename FUNC>
[[nodiscard]] auto simpsons_rule_2d(FUNC f,
                                    PassiveFloat x_min,
                                    PassiveFloat x_max,
                                    PassiveFloat y_min,
                                    PassiveFloat y_max,
                                    size_t nx,
                                    size_t ny) noexcept -> PassiveFloat {
  auto get_w = [](size_t idx, size_t end_idx) -> PassiveFloat {
    if (idx == 0 || idx == end_idx) {
      return 1;
    } else {
      if (idx % 2 == 1) {
        return 4;
      } else {
        return 2;
      }
    }
  };

  if (nx % 2 == 1) { nx += 1; }
  if (ny % 2 == 1) { ny += 1; }

  const PassiveFloat dx = (x_max - x_min) / static_cast<PassiveFloat>(nx);
  const PassiveFloat dy = (y_max - y_min) / static_cast<PassiveFloat>(ny);

  PassiveFloat res = 0;
#pragma omp parallel for reduction(+ : res)
  for (size_t i = 0; i < (nx + 1) * (ny + 1); ++i) {
    const size_t ix       = i / (nx + 1);
    const size_t iy       = i % (nx + 1);
    const PassiveFloat x  = x_min + static_cast<PassiveFloat>(ix) * dx;
    const PassiveFloat y  = y_min + static_cast<PassiveFloat>(iy) * dy;
    const PassiveFloat wx = get_w(ix, nx);
    const PassiveFloat wy = get_w(iy, ny);
    res += wx * wy * f(x, y);
  }
  return res * dx * dy / 9;
}

// -------------------------------------------------------------------------------------------------
auto compare(size_t nx,
             size_t ny,
             PassiveFloat tend,
             const std::tuple<Zap::MatBased::Vector<PassiveFloat>,
                              Zap::MatBased::Vector<PassiveFloat>,
                              Zap::MatBased::Matrix<PassiveFloat>>& hr_res,
             PassiveFloat L1_hr,
             PassiveFloat CFL_safety_factor,
             std::ostream& out) noexcept -> bool {
  const auto cell_res = run_cell_based(nx, ny, tend, CFL_safety_factor);
  const auto mat_res  = run_mat_based(nx, ny, tend, CFL_safety_factor);

  if (!cell_res.has_value()) {
    Igor::Warn("Cell based solver failed for grid of size {}x{}", nx, ny);
    return false;
  }
  if (!mat_res.has_value()) {
    Igor::Warn("Matrix based solver failed for grid of size {}x{}", nx, ny);
    return false;
  }

  const auto [min_cell, max_cell] = min_max_cell_value(cell_res->cells());

  const auto [mat_x, mat_y, mat_u] = *mat_res;
  PassiveFloat min_mat             = mat_u(0, 0);
  PassiveFloat max_mat             = mat_u(0, 0);
  for (Eigen::Index row = 0; row < mat_u.rows(); ++row) {
    for (Eigen::Index col = 0; col < mat_u.cols(); ++col) {
      min_mat = std::min(min_mat, mat_u(row, col));
      max_mat = std::max(max_mat, mat_u(row, col));
    }
  }

  out << nx << 'x' << ny << ":\n";
  out << "min_cell = " << min_cell(0) << '\n';
  out << "max_cell = " << max_cell(0) << '\n';
  out << "min_mat  = " << min_mat << '\n';
  out << "max_mat  = " << max_mat << '\n';

  // - L1 of matrix-based --------------------------------------------------------------------------
  {
    auto f = [&](PassiveFloat x, PassiveFloat y) {
      return std::abs(
          eval_mat_at(x, y, std::get<0>(*mat_res), std::get<1>(*mat_res), std::get<2>(*mat_res)));
    };
    const auto L1 = simpsons_rule_2d(f, X_MIN, X_MAX, Y_MIN, Y_MAX, 2'000, 2'000);
    out << "L1_mat = " << L1 << '\n';
  }

  // - L1 of cell-based ----------------------------------------------------------------------------
  {
    auto f = [&](PassiveFloat x, PassiveFloat y) {
      return std::abs(cell_res->eval(Zap::CellBased::SimCoord<PassiveFloat>{.x = x, .y = y})(0));
    };
    const auto L1 = simpsons_rule_2d(f, X_MIN, X_MAX, Y_MIN, Y_MAX, 2'000, 2'000);
    out << "L1_cell = " << L1 << '\n';
  }

  // - L1 of high-resolution -----------------------------------------------------------------------
  { out << "L1_hr = " << L1_hr << '\n'; }

  // - L1 error mat-based to cell-based ------------------------------------------------------------
  {
    auto f = [&](PassiveFloat x, PassiveFloat y) {
      return std::abs(
          cell_res->eval(Zap::CellBased::SimCoord<PassiveFloat>{.x = x, .y = y})(0) -
          eval_mat_at(x, y, std::get<0>(*mat_res), std::get<1>(*mat_res), std::get<2>(*mat_res)));
    };
    const auto L1_error = simpsons_rule_2d(f, X_MIN, X_MAX, Y_MIN, Y_MAX, 2'000, 2'000);
    out << "L1_error_mat_cell = " << L1_error << '\n';
  }

  // - L1 error high-resolution to cell-based ------------------------------------------------------
  {
    auto f = [&](PassiveFloat x, PassiveFloat y) {
      return std::abs(
          cell_res->eval(Zap::CellBased::SimCoord<PassiveFloat>{.x = x, .y = y})(0) -
          eval_mat_at(x, y, std::get<0>(hr_res), std::get<1>(hr_res), std::get<2>(hr_res)));
    };
    const auto L1_error = simpsons_rule_2d(f, X_MIN, X_MAX, Y_MIN, Y_MAX, 2'000, 2'000);
    out << "L1_error_hr_cell = " << L1_error << '\n';
  }

  // - L1 error high-resolution to mat-based -------------------------------------------------------
  {
    auto f = [&](PassiveFloat x, PassiveFloat y) {
      return std::abs(
          eval_mat_at(x, y, std::get<0>(*mat_res), std::get<1>(*mat_res), std::get<2>(*mat_res)) -
          eval_mat_at(x, y, std::get<0>(hr_res), std::get<1>(hr_res), std::get<2>(hr_res)));
    };
    const auto L1_error = simpsons_rule_2d(f, X_MIN, X_MAX, Y_MIN, Y_MAX, 2'000, 2'000);
    out << "L1_error_hr_mat = " << L1_error << '\n';
  }

  out << "--------------------------------------------------------------------------------\n"
      << std::flush;

  return true;
}

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

  if (args->tend <= 0.0) {
    Igor::Warn("tend must be larger than 0, but is {}.", args->tend);
    return 1;
  }

  Igor::Info("Save results in `" OUTPUT_DIR "`.");
  Igor::Info("run-benchmark     = {}", args->run_benchmark);
  Igor::Info("tend              = {}", args->tend);
  Igor::Info("CFL-safety-factor = {}", args->CFL_safety_factor);

  // - Run high-resolution Godunov solver ----------------------------------------------------------
  const auto hr_nx  = 800;
  const auto hr_ny  = 800;
  const auto hr_res = run_mat_based(hr_nx, hr_ny, args->tend, args->CFL_safety_factor);
  if (!hr_res.has_value()) {
    Igor::Warn("High-resultion ({}x{}) solver failed.", hr_nx, hr_ny);
    return 1;
  }

  auto f_hr = [&](PassiveFloat x, PassiveFloat y) {
    return std::abs(
        eval_mat_at(x, y, std::get<0>(*hr_res), std::get<1>(*hr_res), std::get<2>(*hr_res)));
  };
  const auto L1_hr = simpsons_rule_2d(f_hr, X_MIN, X_MAX, Y_MIN, Y_MAX, 2'000, 2'000);
  // - Run high-resolution Godunov solver ----------------------------------------------------------

  if (args->run_benchmark) {
#ifdef ZAP_STATIC_CUT
    constexpr auto output_file = "./static_cut_results.txt";
#else
    constexpr auto output_file = "./moving_cut_results.txt";
#endif  // ZAP_STATIC_CUT
    std::ofstream out(output_file);
    if (!out) {
      Igor::Warn("Could not open output file `{}`: {}", output_file, std::strerror(errno));
    }

    bool all_success = true;
    // TODO: Find out why 71 does not work
    // TODO: Investigate why 161 has constant max. value 6.24992
    //        -> maybe because a subcell becomes very small
#ifdef ZAP_STATIC_CUT
    constexpr std::array ns = {
        3UZ,  5UZ,  7UZ,   9UZ,   11UZ,  15UZ,  21UZ,  31UZ,  41UZ,  51UZ,  61UZ, /* 71UZ, */
        81UZ, 91UZ, 101UZ, 111UZ, 121UZ, 131UZ, 141UZ, 151UZ, 161UZ, 171UZ, 181UZ, 191UZ,
    };
#else
    constexpr std::array ns = {
        3UZ,  5UZ,  7UZ,   9UZ,   11UZ,  15UZ,  21UZ,  31UZ,  41UZ,  51UZ,  61UZ,  71UZ,
        81UZ, 91UZ, 101UZ, 111UZ, 121UZ, 131UZ, 141UZ, 151UZ, 161UZ, 171UZ, 181UZ, 191UZ,
    };
#endif  // ZAP_STATIC_CUT

    for (size_t n : ns) {
      bool success;
      if (n < 10) {
        success = compare(n, n, args->tend, *hr_res, L1_hr, 0.5, out);
      } else {
        success = compare(n, n, args->tend, *hr_res, L1_hr, args->CFL_safety_factor, out);
      }
      if (success) {
        Igor::Info("Solver for {}x{}-grid finished successfully.", n, n);
      } else {
        Igor::Warn("Solver for {}x{}-grid failed.", n, n);
      }
      all_success = all_success && success;
    }

    Igor::Info("Saved errors to `{}`.", output_file);

    return all_success ? 0 : 1;
  } else {
    Igor::Info("nx                = {}", args->nx);
    Igor::Info("ny                = {}", args->ny);

    const auto success =
        compare(args->nx, args->ny, args->tend, *hr_res, L1_hr, args->CFL_safety_factor, std::cout);
    return success ? 0 : 1;
  }
}