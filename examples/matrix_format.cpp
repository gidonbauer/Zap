#include <iostream>

#include <Eigen/Dense>

#include "Igor.hpp"
#include "ReadMatrixInc.hpp"
#include "WriteMatrixInc.hpp"

auto write() -> int {
  try {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajorBit> mat(5, 3);
    mat << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;

    constexpr auto filename = "output/test_row_major.dat";
    Zap::IncMatrixWriter writer(filename, mat);
    if (!writer.write_data(mat)) {
      return 1;
    }
    mat += Eigen::MatrixXd::Ones(5, 3);
    if (!writer.write_data(mat)) {
      return 1;
    }

    Igor::Info("Wrote `{}`", filename);
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }

  try {
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> mat(5, 3);
    mat << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;

    constexpr auto filename = "output/test_col_major.dat";
    Zap::IncMatrixWriter writer(filename, mat);
    if (!writer.write_data(mat)) {
      return 1;
    }

    if (!writer.write_data(mat * 5)) {
      return 1;
    }
    Igor::Info("Wrote `{}`", filename);
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }

  try {
    constexpr auto filename = "output/test_scalar.dat";
    Zap::IncMatrixWriter<double, 1, 1, 0> writer(filename, 1, 1, 0);
    if (!writer.write_data(4.2)) {
      return 1;
    }

    if (!writer.write_data(22.0 / 7.0)) {
      return 1;
    }

    if (!writer.write_data(1e-8)) {
      return 1;
    }

    if (!writer.write_data(-1.0)) {
      return 1;
    }
    Igor::Info("Wrote `{}`", filename);
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
  return 0;
}

auto read() -> int {
  try {
    constexpr auto filename = "output/test_row_major.dat";
    Igor::Info("Read `{}`", filename);
    Zap::IncMatrixReader<double> reader(filename);

    while (reader.read_next<false>()) {
      for (int row = 0; row < reader.rows(); ++row) {
        for (int col = 0; col < reader.cols(); ++col) {
          std::cout << reader(row, col) << ' ';
        }
        std::cout << '\n';
      }
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }

  try {
    constexpr auto filename = "output/test_col_major.dat";
    Igor::Info("Read `{}`", filename);
    Zap::IncMatrixReader<int> reader(filename);

    while (reader.read_next<false>()) {
      for (int row = 0; row < reader.rows(); ++row) {
        for (int col = 0; col < reader.cols(); ++col) {
          std::cout << reader(row, col) << ' ';
        }
        std::cout << '\n';
      }
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }

  try {
    constexpr auto filename = "output/test_scalar.dat";
    Igor::Info("Read `{}`", filename);
    Zap::IncMatrixReader<double> reader(filename);

    while (reader.read_next<false>()) {
      for (int row = 0; row < reader.rows(); ++row) {
        for (int col = 0; col < reader.cols(); ++col) {
          std::cout << reader(row, col) << ' ';
        }
        std::cout << '\n';
      }
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }

  return 0;
}

auto main(int argc, char** argv) -> int {
  if (argc < 2) {
    Igor::Warn("Usage: {} <write|read>", *argv);
    Igor::Warn("  Missing option, choose write or read.");
    return 1;
  }

  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  if (std::string{argv[1]} == std::string{"write"}) {
    return write();
  } else if (std::string{argv[1]} == std::string{"read"}) {
    return read();
  }
  Igor::Warn("Usage: {} <write|read>", *argv);
  Igor::Warn("  Invalid option `{}`, choose write or read.", argv[1]);
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  return 1;
}
