#include <array>
#include <cmath>
#include <limits>

namespace details {

constexpr double SqrtNewtonRaphson(const double x, const double curr,
                                   const double prev) {
  return curr == prev ? curr
                      : SqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
}

}  // namespace details

namespace {

// Constexpr version of C++ square root.
// Source: https://gist.github.com/alexshtf/eb5128b3e3e143187794
constexpr double Sqrt(const double x) {
  return x >= 0 && x < std::numeric_limits<double>::infinity()
             ? details::SqrtNewtonRaphson(x, x, 0)
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

constexpr std::size_t kPoints = 9;
constexpr std::array<double, kPoints> kXs{
    1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
};
constexpr std::array<double, kPoints> kYs{
    0.16, 0.68, 1.96, 2.79, 3.80, 6.81, 9.50, 15.60, 24.86,
};

constexpr auto kXA = ArithmeticMean(kXs.front(), kXs.back());  // 3.00
constexpr auto kXG = GeometricMean(kXs.front(), kXs.back());   // 2.24
constexpr auto kXH = HarmonicMean(kXs.front(), kXs.back());    // 1.67

constexpr auto kYA = ArithmeticMean(kYs.front(), kYs.back());
constexpr auto kYG = GeometricMean(kYs.front(), kYs.back());
constexpr auto kYH = HarmonicMean(kYs.front(), kYs.back());

// Set according to the graph
constexpr auto kZXA = 4.2;
constexpr auto kZXG = 1.9;
constexpr auto kZXH = 0.9;

constexpr std::size_t kDeltasSize = 9;
constexpr std::array<double, kDeltasSize> kDeltas{
    std::abs(kZXA - kYA), std::abs(kZXG - kYG), std::abs(kZXA - kYG),
    std::abs(kZXG - kYA), std::abs(kZXH - kYA), std::abs(kZXA - kYH),
    std::abs(kZXH - kYH), std::abs(kZXH - kYG), std::abs(kZXG - kYH),
};

}  // namespace

int main() {
}
