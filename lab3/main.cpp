#include <algorithm>
#include <array>
#include <iostream>
#include <limits>

namespace {

constexpr double SqrtNewtonRaphson(const double x, const double curr,
                                   const double prev) {
  return curr == prev ? curr
                      : SqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
}

// Constexpr version of C++ square root.
// Source: https://gist.github.com/alexshtf/eb5128b3e3e143187794
constexpr double Sqrt(const double x) {
  return x >= 0 && x < std::numeric_limits<double>::infinity()
             ? SqrtNewtonRaphson(x, x, 0)
             : std::numeric_limits<double>::quiet_NaN();
}

constexpr double ArithmeticMean(const double a, const double b) {
  return (a + b) / 2;
}

constexpr double GeometricMean(const double a, const double b) {
  return Sqrt(a * b);
}

constexpr double HarmonicMean(const double a, const double b) {
  return 2 / (1 / a + 1 / b);
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

  constexpr auto x_a = ArithmeticMean(kXs.front(), kXs.back());  // 3.00
  constexpr auto x_g = GeometricMean(kXs.front(), kXs.back());   // 2.24
  constexpr auto x_h = HarmonicMean(kXs.front(), kXs.back());    // 1.67

  constexpr auto y_a = ArithmeticMean(kYs.front(), kYs.back());
  constexpr auto y_g = GeometricMean(kYs.front(), kYs.back());
  constexpr auto y_h = HarmonicMean(kYs.front(), kYs.back());

  // Set according to the graph
  constexpr auto z_x_a = 4.2;
  constexpr auto z_x_g = 1.9;
  constexpr auto z_x_h = 0.9;

  constexpr std::size_t kDeltas = 9;
  constexpr std::array<double, kDeltas> deltas{
      std::abs(z_x_a - y_a), std::abs(z_x_g - y_g), std::abs(z_x_a - y_g),
      std::abs(z_x_g - y_a), std::abs(z_x_h - y_a), std::abs(z_x_a - y_h),
      std::abs(z_x_h - y_h), std::abs(z_x_h - y_g), std::abs(z_x_g - y_h),
  };

  const auto min_delta = std::min_element(deltas.cbegin(), deltas.cend());
}
