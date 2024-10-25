#include <filesystem>
#include <numbers>

// #define ZAP_TANGENTIAL_CORRECTION

#include "CellBased/EigenDecomp.hpp"
#include "CellBased/Solver.hpp"
#include "IO/IncCellWriter.hpp"
#include "IO/IncMatrixWriter.hpp"

#include "Igor/Logging.hpp"
#include "Igor/Macros.hpp"
#include "Igor/Timer.hpp"

#include "Common.hpp"

#define OUTPUT_DIR IGOR_STRINGIFY(ZAP_OUTPUT_DIR) "example_full_circle/"

// -------------------------------------------------------------------------------------------------
auto main(int argc, char** argv) -> int {
  using ActiveFloat    = double;
  using PassiveFloat   = double;
  constexpr size_t DIM = 1;

  if (argc < 4) {
    Igor::Warn("Usage: {} <nx> <ny> <tend>", *argv);
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

  const auto nx   = parse_size_t(argv[1]);  // NOLINT
  const auto ny   = parse_size_t(argv[2]);  // NOLINT
  const auto tend = parse_double(argv[3]);  // NOLINT
  assert(nx < std::numeric_limits<int>::max());
  assert(ny < std::numeric_limits<int>::max());

  if (tend <= 0.0) {
    Igor::Warn("tend must be larger than 0, but is {}.", tend);
    return 1;
  }

  Igor::Info("nx = {}", nx);
  Igor::Info("ny = {}", ny);
  Igor::Info("tend = {}", tend);

  const PassiveFloat x_min = 0.0;
  const PassiveFloat x_max = 5.0;
  const PassiveFloat y_min = 0.0;
  const PassiveFloat y_max = 5.0;

  Zap::CellBased::UniformGrid<ActiveFloat, PassiveFloat, DIM> grid(
      x_min, x_max, nx, y_min, y_max, ny);
  // grid.same_value_boundary();
  grid.periodic_boundary();

  const auto x_mid = (x_min + x_max) / 2;
  const auto y_mid = (y_min + y_max) / 2;
  const auto r     = (x_min + x_max + y_min + y_max) / 8;

  auto u0 = [=](ActiveFloat x, ActiveFloat y) -> ActiveFloat {
    // return r *
    return (std::pow(x - x_mid, 2) + std::pow(y - y_mid, 2)) *
           static_cast<ActiveFloat>((std::pow(x - x_mid, 2) + std::pow(y - y_mid, 2)) <=
                                    std::pow(r, 2));
  };

  auto init_shock = [=]<typename T>(T t) -> Zap::CellBased::SimCoord<T> {
    t /= 4;
    // assert(t >= 0 && t <= 1);
    return {
        r * std::cos(2 * std::numbers::pi_v<PassiveFloat> * t) + x_mid,
        r * std::sin(2 * std::numbers::pi_v<PassiveFloat> * t) + y_mid,
    };
  };

  const size_t N = (nx + ny) / 2;
  std::vector<Zap::CellBased::SimCoord<ActiveFloat>> init_shock_points(N);
  const auto step = 1 / static_cast<PassiveFloat>(N - 1);
  for (size_t i = 0; i < N; ++i) {
    init_shock_points[i] = init_shock(static_cast<PassiveFloat>(i) * step);
  }

  IGOR_TIME_SCOPE("Cutting the grid") {
    if (!grid.cut_piecewise_linear<Zap::CellBased::ExtendType::NEAREST>(init_shock_points)) {
      return 1;
    }
  }

  // grid.fill_center(u0);
  grid.fill_four_point(u0);

  // grid.dump_cells(std::cout);

  constexpr auto u_file = OUTPUT_DIR "u.grid";
  Zap::IO::IncCellWriter<ActiveFloat, PassiveFloat, DIM> grid_writer{u_file, grid};

  constexpr auto t_file = OUTPUT_DIR "t.mat";
  Zap::IO::IncMatrixWriter<PassiveFloat, 1, 1, 0> t_writer(t_file, 1, 1, 0);

// #define SAVE_ONLY_INITIAL_STATE
#ifdef SAVE_ONLY_INITIAL_STATE
  if (!grid_writer.write_data(grid)) { return 1; }
  if (!t_writer.write_data(Float{0})) { return 1; }
#else
  IGOR_TIME_SCOPE("Solver") {
    auto solver = Zap::CellBased::make_solver<Zap::CellBased::ExtendType::NEAREST>(
        Zap::CellBased::SingleEq::A{}, Zap::CellBased::SingleEq::B{});
    const auto res =
        solver.solve(grid, static_cast<PassiveFloat>(tend), grid_writer, t_writer, 0.25);
    if (!res.has_value()) {
      Igor::Warn("Solver failed.");
      return 1;
    }

    const auto final_shock = res->get_shock_curve();
    // Igor::Info("Final shock curve: {}", final_shock);
    const auto mean_shock_x = std::transform_reduce(std::cbegin(final_shock),
                                                    std::cend(final_shock),
                                                    ActiveFloat{0},
                                                    std::plus<>{},
                                                    [](const auto& p) { return p.x; }) /
                              static_cast<PassiveFloat>(final_shock.size());

    const auto std_dev_shock_x = std::sqrt(
        std::transform_reduce(std::cbegin(final_shock),
                              std::cend(final_shock),
                              ActiveFloat{0},
                              std::plus<>{},
                              [=](const auto& p) { return std::pow(p.x - mean_shock_x, 2); }) /
        static_cast<PassiveFloat>(final_shock.size()));

    Igor::Info("mean_shock_x = {}", mean_shock_x);
    Igor::Info("std_dev_shock_x = {}", std_dev_shock_x);
  }
  Igor::Info("Solver finished successfully.");
#endif
  Igor::Info("Saved grid to {}.", u_file);
  Igor::Info("Saved time steps to {}.", t_file);
}
