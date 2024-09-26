#include <atomic>
#include <cmath>
#include <filesystem>

#define ZAP_STATIC_CUT
#include "CellBased/EigenDecomp.hpp"
#include "CellBased/Grid.hpp"
#include "CellBased/Solver.hpp"

#include "MatBased/BoundaryConditions.hpp"
#include "MatBased/Solver.hpp"

#include "IO/IncCellWriter.hpp"
#include "IO/IncMatrixWriter.hpp"

#include "Igor.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "example_static_shock/"

// - Setup -----------------------------------------------------------------------------------------
using Float           = double;
constexpr size_t DIM  = 1;
constexpr Float X_MIN = 0.0;
constexpr Float X_MAX = 5.0;
constexpr Float Y_MIN = 0.0;
constexpr Float Y_MAX = 5.0;

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto parse_double(const char* cstr) -> double {
  char* end        = nullptr;
  const double val = std::strtod(cstr, &end);
  if (val == HUGE_VAL) {
    Igor::Panic("Could not parse string `{}` to double: Out of range.", cstr);
  }
  if (cstr == end) { Igor::Panic("Could not parse string `{}` to double.", cstr); }
  return val;
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto sqr(Float x) noexcept -> Float { return x * x; }

// -------------------------------------------------------------------------------------------------
[[nodiscard]] constexpr auto u0(Float x, Float y) noexcept -> Float {
  static_assert(DIM == 1);
  return (sqr(x - X_MIN) + sqr(y - Y_MIN)) *
         static_cast<Float>((sqr(x - X_MIN) + sqr(y - Y_MIN)) <=
                            sqr((X_MIN + X_MAX + Y_MIN + Y_MAX) / 4));
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto init_shock(Float t) noexcept -> Zap::CellBased::Point<Float> {
  const auto r = (X_MIN + X_MAX + Y_MIN + Y_MAX) / 4;
  return Zap::CellBased::Point<Float>{
      r * std::cos(std::numbers::pi_v<Float> / 2 * t),
      r * std::sin(std::numbers::pi_v<Float> / 2 * t),
  };
};

// -------------------------------------------------------------------------------------------------
auto run_cell_based(size_t nx,
                    size_t ny,
                    Float tend) noexcept -> std::optional<Zap::CellBased::UniformGrid<Float, DIM>> {
  Zap::CellBased::UniformGrid<Float, DIM> grid(X_MIN, X_MAX, nx, Y_MIN, Y_MAX, ny);
  grid.same_value_boundary();
  if (!grid.cut_curve(init_shock)) { return std::nullopt; }
  grid.fill_four_point(u0);

  Zap::CellBased::Solver solver(Zap::CellBased::SingleEq::A{}, Zap::CellBased::SingleEq::B{});

  const std::string u_filename =
      OUTPUT_DIR "u_cell_based_" + std::to_string(nx) + "_" + std::to_string(ny) + ".grid";
  Zap::IO::IncCellWriter u_writer(u_filename, grid);
  const std::string t_filename =
      OUTPUT_DIR "t_cell_based_" + std::to_string(nx) + "_" + std::to_string(ny) + ".mat";
  Zap::IO::IncMatrixWriter<Float, 1, 1, 0> t_writer(t_filename, 1, 1, 0);

  return solver.solve(grid, tend, u_writer, t_writer);
}

// -------------------------------------------------------------------------------------------------
auto run_mat_based(size_t nx,
                   size_t ny,
                   Float tend) noexcept -> std::optional<std::tuple<Zap::MatBased::Vector<Float>,
                                                                    Zap::MatBased::Vector<Float>,
                                                                    Zap::MatBased::Matrix<Float>>> {
  const auto dx = (X_MAX - X_MIN) / static_cast<Float>(nx);
  const auto dy = (Y_MAX - Y_MIN) / static_cast<Float>(ny);

  Zap::MatBased::Vector<Float> x(nx);
  std::generate(std::begin(x), std::end(x), [i = 0, dx]() mutable {
    return X_MIN + dx * (static_cast<Float>(i++) + static_cast<Float>(0.5));
  });
  Zap::MatBased::Vector<Float> y(static_cast<int>(ny));
  std::generate(std::begin(y), std::end(y), [i = 0, dy]() mutable {
    return Y_MIN + dy * (static_cast<Float>(i++) + static_cast<Float>(0.5));
  });

  Zap::MatBased::Matrix<Float> u0_grid(static_cast<int>(ny), static_cast<int>(nx));
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

  auto boundary = [](Zap::MatBased::Matrix<Float>& u_next,
                     const Zap::MatBased::Matrix<Float>& u_curr,
                     Float dt,
                     Float local_dx,
                     Float local_dy,
                     auto numerical_flux_x,
                     auto numerical_flux_y) {
    return Zap::MatBased::equal_value_boundary(
        u_next, u_curr, dt, local_dx, local_dy, numerical_flux_x, numerical_flux_y);
  };

  const std::string u_filename =
      OUTPUT_DIR "u_mat_based_" + std::to_string(nx) + "_" + std::to_string(ny) + ".mat";
  Zap::IO::IncMatrixWriter u_writer(u_filename, u0_grid);

  const std::string t_filename =
      OUTPUT_DIR "t_mat_based_" + std::to_string(nx) + "_" + std::to_string(ny) + ".mat";
  Zap::IO::IncMatrixWriter<Float, 1, 1, 0> t_writer(t_filename, 1, 1, 0);

  const auto res =
      Zap::MatBased::solve_2d_burgers(x, y, u0_grid, tend, boundary, u_writer, t_writer);

  if (!res.has_value()) { return std::nullopt; }
  return std::make_tuple(x, y, *res);
}

// -------------------------------------------------------------------------------------------------
[[nodiscard]] auto
min_max_cell_value(const std::vector<Zap::CellBased::Cell<Float, DIM>>& cells) noexcept
    -> std::pair<Eigen::Vector<Float, DIM>, Eigen::Vector<Float, DIM>> {
  Eigen::Vector<Float, DIM> min = [cell = cells.front()] {
    assert(cell.is_cartesian() || cell.is_cut());
    if (cell.is_cartesian()) {
      return cell.get_cartesian().value;
    } else {
      Eigen::Vector<Float, DIM> res;
      for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(DIM); ++i) {
        res(i) = std::min(cell.get_cut().left_value(i), cell.get_cut().right_value(i));
      }
      return res;
    }
  }();

  Eigen::Vector<Float, DIM> max = [cell = cells.front()] {
    assert(cell.is_cartesian() || cell.is_cut());
    if (cell.is_cartesian()) {
      return cell.get_cartesian().value;
    } else {
      Eigen::Vector<Float, DIM> res;
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
        min(i) =
            std::min(min(i), std::min(cell.get_cut().left_value(i), cell.get_cut().right_value(i)));
        max(i) =
            std::max(max(i), std::max(cell.get_cut().left_value(i), cell.get_cut().right_value(i)));
      }
    }
  }

  return std::make_pair(min, max);
}

// -------------------------------------------------------------------------------------------------
constexpr auto eval_mat_at(Float x,
                           Float y,
                           const Zap::MatBased::Vector<Float>& xs,
                           const Zap::MatBased::Vector<Float>& ys,
                           const Zap::MatBased::Matrix<Float>& u) noexcept -> Float {
  assert(ys.size() == u.rows());
  assert(xs.size() == u.cols());
  assert(xs.size() >= 2);
  assert(ys.size() >= 2);

  if (std::isnan(x) || std::isnan(y)) { return std::numeric_limits<Float>::quiet_NaN(); }

  const auto dx    = xs(1) - xs(0);
  const auto x_idx = static_cast<Eigen::Index>(std::clamp(
      std::round((x - xs(0)) / dx), static_cast<Float>(0), static_cast<Float>(xs.rows() - 1)));

  const auto dy    = ys(1) - ys(0);
  const auto y_idx = static_cast<Eigen::Index>(std::clamp(
      std::round((y - ys(0)) / dy), static_cast<Float>(0), static_cast<Float>(ys.rows() - 1)));

  return u(y_idx, x_idx);
}

// -------------------------------------------------------------------------------------------------
template <typename FUNC>
[[nodiscard]] auto simpsons_rule_2d(FUNC f,
                                    Float x_min,
                                    Float x_max,
                                    Float y_min,
                                    Float y_max,
                                    size_t nx,
                                    size_t ny) noexcept -> Float {
  auto get_w = [](size_t idx, size_t end_idx) -> Float {
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

  Float dx = (x_max - x_min) / static_cast<Float>(nx);
  Float dy = (y_max - y_min) / static_cast<Float>(ny);

  Float res = 0;
#pragma omp parallel for reduction(+ : res)
  for (size_t i = 0; i < (nx + 1) * (ny + 1); ++i) {
    const size_t ix = i / (nx + 1);
    const size_t iy = i % (nx + 1);
    const Float x   = x_min + static_cast<Float>(ix) * dx;
    const Float y   = y_min + static_cast<Float>(iy) * dy;
    const Float wx  = get_w(ix, nx);
    const Float wy  = get_w(iy, ny);
    res += wx * wy * f(x, y);
  }
  return res * dx * dy / 9;
}

// -------------------------------------------------------------------------------------------------
auto compare(size_t nx,
             size_t ny,
             Float tend,
             const std::tuple<Zap::MatBased::Vector<Float>,
                              Zap::MatBased::Vector<Float>,
                              Zap::MatBased::Matrix<Float>>& hr_res,
             Float L1_hr,
             std::ostream& out) noexcept -> bool {
  const auto cell_res = run_cell_based(nx, ny, tend);
  const auto mat_res  = run_mat_based(nx, ny, tend);

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
  Float min_mat                    = mat_u(0, 0);
  Float max_mat                    = mat_u(0, 0);
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
    auto f = [&](Float x, Float y) {
      return std::abs(
          eval_mat_at(x, y, std::get<0>(*mat_res), std::get<1>(*mat_res), std::get<2>(*mat_res)));
    };
    const auto L1 = simpsons_rule_2d(f, X_MIN, X_MAX, Y_MIN, Y_MAX, 2'000, 2'000);
    out << "L1_mat = " << L1 << '\n';
  }

  // - L1 of cell-based ----------------------------------------------------------------------------
  {
    auto f        = [&](Float x, Float y) { return std::abs(cell_res->eval({x, y})(0)); };
    const auto L1 = simpsons_rule_2d(f, X_MIN, X_MAX, Y_MIN, Y_MAX, 2'000, 2'000);
    out << "L1_cell = " << L1 << '\n';
  }

  // - L1 of high-resolution -----------------------------------------------------------------------
  { out << "L1_hr = " << L1_hr << '\n'; }

  // - L1 error mat-based to cell-based ------------------------------------------------------------
  {
    auto f = [&](Float x, Float y) {
      return std::abs(
          cell_res->eval({x, y})(0) -
          eval_mat_at(x, y, std::get<0>(*mat_res), std::get<1>(*mat_res), std::get<2>(*mat_res)));
    };
    const auto L1_error = simpsons_rule_2d(f, X_MIN, X_MAX, Y_MIN, Y_MAX, 2'000, 2'000);
    out << "L1_error_mat_cell = " << L1_error << '\n';
  }

  // - L1 error high-resolution to cell-based ------------------------------------------------------
  {
    auto f = [&](Float x, Float y) {
      return std::abs(
          cell_res->eval({x, y})(0) -
          eval_mat_at(x, y, std::get<0>(hr_res), std::get<1>(hr_res), std::get<2>(hr_res)));
    };
    const auto L1_error = simpsons_rule_2d(f, X_MIN, X_MAX, Y_MIN, Y_MAX, 2'000, 2'000);
    out << "L1_error_hr_cell = " << L1_error << '\n';
  }

  // - L1 error high-resolution to mat-based -------------------------------------------------------
  {
    auto f = [&](Float x, Float y) {
      return std::abs(
          eval_mat_at(x, y, std::get<0>(*mat_res), std::get<1>(*mat_res), std::get<2>(*mat_res)) -
          eval_mat_at(x, y, std::get<0>(hr_res), std::get<1>(hr_res), std::get<2>(hr_res)));
    };
    const auto L1_error = simpsons_rule_2d(f, X_MIN, X_MAX, Y_MIN, Y_MAX, 2'000, 2'000);
    out << "L1_error_hr_mat = " << L1_error << '\n';
  }

  out << "--------------------------------------------------------------------------------\n";

  return true;
}

auto main(int argc, char** argv) -> int {
  if (argc < 2) {
    Igor::Warn("Usage: {} <tend>", *argv);
    return 1;
  }

  {
    std::error_code ec;
    std::filesystem::create_directories(OUTPUT_DIR, ec);
    if (ec) {
      Igor::Warn("Could not create directory `{}`: {}", OUTPUT_DIR, ec.message());
      return 1;
    }
  }

  const auto tend = static_cast<Float>(parse_double(argv[1]));  // NOLINT

  if (tend <= 0.0) {
    Igor::Warn("tend must be larger than 0, but is {}.", tend);
    return 1;
  }

  Igor::Info("tend = {}", tend);
  Igor::Info("Save results in `" OUTPUT_DIR "`.");

  constexpr auto output_file = "./static_cut_results.txt";
  std::ofstream out(output_file);
  if (!out) {
    Igor::Warn("Could not open output file `{}`: {}", output_file, std::strerror(errno));
  }

  bool all_success = true;
  // TODO: Find out why 71 does not work
  // TODO: Investigate why 161 has constant max. value 6.24992
  //        -> maybe because a subcell becomes very small
  constexpr std::array ns = {
      3UZ,  5UZ,  7UZ,   9UZ,   11UZ,  15UZ,  21UZ,  31UZ,  41UZ,  51UZ,  61UZ, /* 71UZ, */
      81UZ, 91UZ, 101UZ, 111UZ, 121UZ, 131UZ, 141UZ, 151UZ, 161UZ, 171UZ, 181UZ, 191UZ,
  };

  const auto hr_nx  = 800;
  const auto hr_ny  = 800;
  const auto hr_res = run_mat_based(hr_nx, hr_ny, tend);
  if (!hr_res.has_value()) {
    Igor::Warn("High-resultion ({}x{}) solver failed.", hr_nx, hr_ny);
    return 1;
  }

  auto f_hr = [&](Float x, Float y) {
    return std::abs(
        eval_mat_at(x, y, std::get<0>(*hr_res), std::get<1>(*hr_res), std::get<2>(*hr_res)));
  };
  const auto L1_hr = simpsons_rule_2d(f_hr, X_MIN, X_MAX, Y_MIN, Y_MAX, 2'000, 2'000);

  for (size_t n : ns) {
    const auto success = compare(n, n, tend, *hr_res, L1_hr, out);
    all_success        = all_success && success;
  }

  Igor::Info("Saved errors to `{}`.", output_file);

  return all_success ? 0 : 1;
}
