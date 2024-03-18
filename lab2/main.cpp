#include <cassert>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numbers>
#include <vector>

namespace {

constexpr auto kEps = 1e-3;
constexpr std::size_t kW = 15;

struct Integral {
  double a, b;
  std::function<double(const double x)> f;
  double val;
};

struct Approximation {
  std::size_t n;
  double val;
  double r;
};

class IntegrationMethod {
 public:
  virtual ~IntegrationMethod() = default;

  virtual double Accuracy() const = 0;

  virtual double Calculate(const Integral& i, const std::size_t n) const = 0;

  virtual std::string_view Name() const = 0;
};

class RectangleMethod final : public IntegrationMethod {
 public:
  double Accuracy() const noexcept override { return 2; }

  double Calculate(const Integral& i,
                   const std::size_t n) const noexcept override {
    const auto h = (i.b - i.a) / n;

    auto f_sum = 0.0;
    auto x = i.a + 0.5 * h;
    for (std::size_t j = 0; j < n; ++j) {
      f_sum += i.f(x);
      x += h;
    }

    return f_sum * h;
  }

  std::string_view Name() const noexcept override { return "Rectangle"; }
};

class TrapezoidMethod final : public IntegrationMethod {
 public:
  double Accuracy() const noexcept override { return 2; }

  double Calculate(const Integral& i,
                   const std::size_t n) const noexcept override {
    const auto h = (i.b - i.a) / n;

    auto f_sum = 0.0;
    auto x = i.a + h;
    for (std::size_t j = 1; j < n; ++j) {
      f_sum += i.f(x);
      x += h;
    }

    return h * (f_sum + 0.5 * (i.f(i.a) + i.f(i.b)));
  }

  std::string_view Name() const noexcept override { return "Trapezoid"; }
};

class SimpsonMethod final : public IntegrationMethod {
 public:
  double Accuracy() const noexcept override { return 4; }

  double Calculate(const Integral& i,
                   const std::size_t n) const noexcept override {
    const auto h = (i.b - i.a) / n;

    auto f_sum1 = 0.0;
    auto x = i.a + 0.5 * h;
    for (std::size_t j = 0; j < n; ++j) {
      f_sum1 += i.f(x);
      x += h;
    }

    auto f_sum2 = 0.0;
    x = i.a + h;
    for (std::size_t j = 1; j < n; ++j) {
      f_sum2 += i.f(x);
      x += h;
    }

    return h / 6 * (i.f(i.a) + i.f(i.b) + 4 * f_sum1 + 2 * f_sum2);
  }

  std::string_view Name() const noexcept override { return "Simpson"; }
};

void PrintTable(const std::vector<std::unique_ptr<IntegrationMethod>>& methods,
                const std::vector<Approximation>& approxs) {
  // clang-format off
  std::cout << std::setw(kW) << "Method" 
            << std::setw(kW) << "n"
            << std::setw(kW) << "I"
            << std::setw(kW) << "R"
            << std::setw(kW) << "I + R"
            << '\n';
  
  auto end = methods.size();
  assert(end == approxs.size());

  for (std::size_t i = 0; i < end; ++i) {
    const auto& method = *methods[i];
    const auto& approx = approxs[i];

    std::cout << std::setw(kW) << method.Name()
              << std::setw(kW) << approx.n
              << std::setw(kW) << approx.val
              << std::setw(kW) << approx.r
              << std::setw(kW) << approx.val + approx.r
              << '\n';
  }
  // clang-format on
}

double CalculateRichardson(const double val, const double val_freq,
                           const std::size_t k) {
  return (val_freq - val) / ((2 << k) - 1);
}

Approximation CalculateApproximation(const Integral& i,
                                     const IntegrationMethod& method,
                                     const double eps) noexcept {
  std::size_t n = 1;
  double val_freq = method.Calculate(i, n);
  double val, r;

  do {
    n <<= 1;
    val = val_freq;
    val_freq = method.Calculate(i, n);
    r = CalculateRichardson(val, val_freq, method.Accuracy());
  } while (std::abs(r) >= eps);

  return {n, val_freq, r};
}

}  // namespace

int main() {
  const Integral i{0, 1, [](const double x) { return std::exp(x); },
                   std::numbers::e - 1};

  std::vector<std::unique_ptr<IntegrationMethod>> methods;
  methods.push_back(std::make_unique<RectangleMethod>());
  methods.push_back(std::make_unique<TrapezoidMethod>());
  methods.push_back(std::make_unique<SimpsonMethod>());

  std::vector<Approximation> approxs;
  approxs.reserve(methods.size());

  for (const auto& method : methods) {
    approxs.push_back(CalculateApproximation(i, *method, kEps));
  }

  std::cout << "Actual value: " << i.val << '\n';
  PrintTable(methods, approxs);
}
