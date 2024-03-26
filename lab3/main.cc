// clang-format off
#include <sprout/math/exp.hpp>
#include <sprout/math/log.hpp>
#include <sprout/math/pow.hpp>
#include <sprout/math/sqrt.hpp>
// clang-format on

#include <array>
#include <cmath>
#include <iomanip>
#include <limits>

namespace {

using Deltas = std::array<double, 9>;

struct Parameters {
  double x_a, x_g, x_h;
  double y_a, y_g, y_h;
  double z_x_a, z_x_g, z_x_h;
  Deltas deltas;
};

struct Matrix {
  double a11, a12, a21, a22;
};

struct Vector {
  double b1, b2;
};

struct AugmentedMatrix {
  Matrix m;
  Vector v;
};

constexpr double ArithmeticMean(const double a, const double b) {
  return (a + b) / 2;
}

constexpr double GeometricMean(const double a, const double b) {
  return sprout::sqrt(a * b);
}

constexpr double HarmonicMean(const double a, const double b) {
  return 2 / (1 / a + 1 / b);
}

template <std::size_t N>
constexpr AugmentedMatrix BuildAugmentedMatrix(
    const std::array<double, N>& xs, const std::array<double, N>& ys) {
  Matrix m{};
  Vector v{};

  double ln_x = 0, ln_y = 0;
  for (std::size_t i = 0; i < N; ++i) {
    ln_x = sprout::log(xs[i]);
    ln_y = sprout::log(ys[i]);

    m.a11 += ln_x * ln_x;
    m.a12 += ln_x;
    m.a21 += ln_x;

    v.b1 += ln_y * ln_x;
    v.b2 += ln_y;
  }

  m.a22 = N;

  return {m, v};
}

constexpr Vector CalculateSolution(const AugmentedMatrix& am) {
  const auto det = am.m.a11 * am.m.a22 - am.m.a21 * am.m.a12;

  // Transposed matrix of algebraic complements
  const Matrix com{am.m.a22, -am.m.a12, -am.m.a21, am.m.a11};

  return {
      (com.a11 * am.v.b1 + com.a12 * am.v.b2) / det,
      (com.a21 * am.v.b1 + com.a22 * am.v.b2) / det,
  };
}

template <std::size_t N>
constexpr double CalculateStandartDeviation(const double b, const double ln_a,
                                            const std::array<double, N>& xs,
                                            const std::array<double, N>& ys) {
  double sum = 0;
  for (std::size_t i = 0; i < N; ++i) {
    sum += sprout::pow(ln_a + b * sprout::log(xs[i]) - sprout::log(ys[i]), 2);
  }

  return sprout::sqrt(sum / N);
}

template <std::size_t N>
constexpr Parameters CalculateParameters(const std::array<double, N>& xs,
                                         const std::array<double, N>& ys) {
  const auto x_a = ArithmeticMean(xs.front(), xs.back());
  const auto x_g = GeometricMean(xs.front(), xs.back());
  const auto x_h = HarmonicMean(xs.front(), xs.back());

  const auto y_a = ArithmeticMean(ys.front(), ys.back());
  const auto y_g = GeometricMean(ys.front(), ys.back());
  const auto y_h = HarmonicMean(ys.front(), ys.back());

  // Set according to the graph
  const auto kZXA = 4.2;
  const auto kZXG = 1.9;
  const auto kZXH = 0.9;

  const Deltas deltas{
      std::abs(kZXA - y_a), std::abs(kZXG - y_g), std::abs(kZXA - y_g),
      std::abs(kZXG - y_a), std::abs(kZXH - y_a), std::abs(kZXA - y_h),
      std::abs(kZXH - y_h), std::abs(kZXH - y_g), std::abs(kZXG - y_h),
  };

  return {x_a, x_g, x_h, y_a, y_g, y_h, kZXA, kZXG, kZXH, deltas};
}

template <std::size_t N>
void PrintCoordinates(const std::array<double, N>& xs,
                      const std::array<double, N>& ys) {
  static constexpr std::size_t kWidth = 6;

  std::cout << "Coordinates:\n\n";

  std::cout << std::setw(kWidth) << "i";
  for (std::size_t i = 0; i < N; ++i) {
    std::cout << std::setw(kWidth) << i;
  }
  std::cout << '\n';

  std::cout << std::setw(kWidth) << "x_i";
  for (std::size_t i = 0; i < N; ++i) {
    std::cout << std::setw(kWidth) << xs[i];
  }
  std::cout << '\n';

  std::cout << std::setw(kWidth) << "y_i";
  for (std::size_t i = 0; i < N; ++i) {
    std::cout << std::setw(kWidth) << ys[i];
  }
  std::cout << "\n\n";
}

void PrintParameters(const Parameters& ps) {
  std::cout << "Parameters:\n\n";

  std::cout << "x_a = " << ps.x_a << '\n';
  std::cout << "x_g = " << ps.x_g << '\n';
  std::cout << "x_h = " << ps.x_h << "\n\n";

  std::cout << "y_a = " << ps.y_a << '\n';
  std::cout << "y_g = " << ps.y_g << '\n';
  std::cout << "y_h = " << ps.y_h << "\n\n";

  std::cout << "z(x_a) = " << ps.z_x_a << '\n';
  std::cout << "z(x_g) = " << ps.z_x_g << '\n';
  std::cout << "z(x_h) = " << ps.z_x_h << "\n\n";

  for (std::size_t i = 0, end = ps.deltas.size(); i < end; ++i) {
    std::cout << "delta_" << i + 1 << " = " << ps.deltas[i] << '\n';
  }
  std::cout << '\n';
}

void PrintAugmentedMatrix(const AugmentedMatrix& am) {
  std::cout << "Augmented matrix of the equations system:\n\n";

  std::cout << "a11 = " << am.m.a11 << '\n';
  std::cout << "a21 = " << am.m.a21 << '\n';
  std::cout << "a12 = " << am.m.a12 << '\n';
  std::cout << "a22 = " << am.m.a22 << "\n\n";

  std::cout << "b1 = " << am.v.b1 << '\n';
  std::cout << "b2 = " << am.v.b2 << "\n\n";
}

void PrintSolution(const Vector& v) {
  std::cout << "Solution:\n\n";

  std::cout << "ln(a) = " << v.b2 << ", a = " << sprout::exp(v.b2) << '\n';
  std::cout << "b = " << v.b1 << "\n\n";
}

void PrintDeviation(const double d) {
  std::cout << "Standart deviation:\n\n";

  std::cout << "Δ = " << d << '\n';
}

}  // namespace

int main() {
  constexpr std::size_t kN = 9;
  constexpr std::array<double, kN> kXs{
      1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
  };
  constexpr std::array<double, kN> kYs{
      0.16, 0.68, 1.96, 2.79, 3.80, 6.81, 9.50, 15.60, 24.86,
  };
  PrintCoordinates(kXs, kYs);

  constexpr auto parameters = CalculateParameters(kXs, kYs);
  PrintParameters(parameters);

  constexpr auto augmented_matrix = BuildAugmentedMatrix(kXs, kYs);
  PrintAugmentedMatrix(augmented_matrix);

  constexpr auto solution = CalculateSolution(augmented_matrix);
  PrintSolution(solution);

  constexpr auto deviation =
      CalculateStandartDeviation(solution.b1, solution.b2, kXs, kYs);
  PrintDeviation(deviation);
}
