#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>

namespace {

using Segment = std::pair<double, double>;
using Function = std::function<double(const double x)>;

int Sgn(const double x) { return x < 0 ? -1 : x > 0; }

struct Result final {
  std::size_t k;
  double x;
};

class Method {
 public:
  virtual ~Method() = default;

  virtual Result FindRoot(const Segment segment, const double epsilon) = 0;
  virtual std::string_view get_name() const = 0;
};

class SplitMethod final : public Method {
 public:
  SplitMethod(Function f) : f_(std::move(f)) {}

  Result FindRoot(const Segment segment, const double epsilon) override {
    auto [a, b] = segment;
    auto k = std::size_t{};
    auto x = (a + b) / 2;

    while (b - a > 2 * epsilon) {
      if (f_(a) * f_(x) < 0) {
        b = x;
      } else {
        a = x;
      }

      x = (a + b) / 2;
      ++k;
    }

    return {k, x};
  }

  std::string_view get_name() const override { return kName; }

 private:
  static constexpr std::string_view kName = "Dichotomy";

  Function f_;
};

class NewtonMethod final : public Method {
 public:
  NewtonMethod(Function f, Function df_dx, Function d2f_dx2)
      : f_(std::move(f)),
        df_dx_(std::move(df_dx)),
        d2f_dx2_(std::move(d2f_dx2)) {}

  Result FindRoot(const Segment segment, const double epsilon) override {
    auto [a, b] = segment;
    auto k = std::size_t{};
    auto x_curr = (f_(a) * d2f_dx2_(a) > 0) ? a : b;
    auto x_prev = 0.0;

    do {
      x_prev = x_curr;
      x_curr = x_prev - f_(x_prev) / df_dx_(x_prev);
      ++k;
    } while (f_(x_curr) * f_(x_curr + Sgn(x_curr - x_prev) * epsilon) >= 0);

    return {k, x_curr};
  }

  std::string_view get_name() const override { return kName; }

 private:
  static constexpr std::string_view kName = "Newton";

  Function f_, df_dx_, d2f_dx2_;
};

}  // namespace

int main() {
  constexpr std::size_t kColW = 18;

  {
    std::cout << "Dichotomy and Newton nonlinear equation methods comparison.\n"
              << "* The columns store {iterations, root} on the segment.\n";

    constexpr auto kEpsilon = 1e-3;

    constexpr auto kF = [](const double x) {
      return 2 * std::pow(x, 3) + 9 * std::pow(x, 2) - 21;
    };
    constexpr auto kDfDx = [](const double x) {
      return 6 * std::pow(x, 2) + 18 * x;
    };
    constexpr auto kD2fDx2 = [](const double x) { return 12 * x + 18; };

    const auto segments =
        std::vector<Segment>{{-4.0, -3.5}, {-2.5, -2.0}, {1.0, 1.5}};

    auto methods = std::vector<std::unique_ptr<Method>>{};
    methods.push_back(std::make_unique<SplitMethod>(kF));
    methods.push_back(std::make_unique<NewtonMethod>(kF, kDfDx, kD2fDx2));

    std::cout << std::setw(kColW) << "Method";
    for (auto&& [a, b] : segments) {
      std::ostringstream oss;
      oss << "[" << a << ", " << b << "]";
      std::cout << std::setw(kColW) << oss.str();
    }
    std::cout << "\n";

    for (auto&& method : methods) {
      std::cout << std::setw(kColW) << method->get_name();
      for (auto&& segment : segments) {
        const auto [k, x] = method->FindRoot(segment, kEpsilon);
        std::ostringstream oss;
        oss << "{" << k << ", " << x << "}";
        std::cout << std::setw(kColW) << oss.str();
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }

  {
    std::cout << "Newton nonlinear equation system method.\n";
    constexpr auto kF = [](const double x, const double y) {
      return std::cos(x) + y - 1.5;
    };
    constexpr auto kG = [](const double x, const double y) {
      return 2 * x - std::sin(y - 0.5) - 1;
    };

    constexpr auto kDfDx = [](const double x, const double y) {
      return -std::sin(x);
    };
    constexpr auto kDfDy = [](const double x, const double y) { return 1; };
    constexpr auto kDgDx = [](const double x, const double y) { return 2; };
    constexpr auto kDgDy = [](const double x, const double y) {
      return -std::cos(y - 0.5);
    };

    constexpr auto JacDetRev = [=](const double x, const double y) {
      return 1 / (kDfDx(x, y) * kDgDy(x, y) - kDfDy(x, y) * kDgDx(x, y));
    };

    constexpr auto kEpsilon = 1e-2;
    constexpr auto kX0 = 0.582;
    constexpr auto kY0 = 0.665;

    auto x_c = 0.5, y_c = 0.6;
    double x_p, y_p, jdr, f, g;
    std::size_t k = 0;

    do {
      x_p = x_c;
      y_p = y_c;

      jdr = JacDetRev(x_p, y_p);
      f = kF(x_p, y_p);
      g = kG(x_p, y_p);

      x_c = x_p - jdr * (kDgDy(x_p, y_p) * f - kDfDy(x_p, y_p) * g);
      y_c = y_p - jdr * (-kDgDx(x_p, y_p) * f + kDfDx(x_p, y_p) * g);

      ++k;
    } while (std::max(std::abs(x_c - x_p), std::abs(y_c - y_p)) >= kEpsilon);

    std::cout << "Analytical solution: (" << kX0 << ", " << kY0 << ").\n";
    std::cout << "Newton method's solution: (" << x_c << ", " << y_c << ").\n";
    std::cout << "Iterations: " << k << ".\n";
    std::cout << "Error: (" << std::abs(kX0 - x_c) << ", "
              << std::abs(kY0 - y_c) << ").\n";
  }
}
