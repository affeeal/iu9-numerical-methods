#include <cmath>
#include <functional>
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

 private:
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

 private:
  Function f_, df_dx_, d2f_dx2_;
};

}  // namespace

int main() {
  constexpr auto kEpsilon = 1e-3;

  const auto f = [](const double x) {
    return 2 * std::pow(x, 3) + 9 * std::pow(x, 2) - 21;
  };
  const auto df_dx = [](const double x) { return 6 * std::pow(x, 2) + 18 * x; };
  const auto d2f_dx2 = [](const double x) { return 12 * x + 18; };

  const auto segments =
      std::vector<Segment>{{-4.0, -3.5}, {-2.5, -2.0}, {1.0, 1.5}};

  auto methods = std::vector<std::unique_ptr<Method>>{};
  methods.push_back(std::make_unique<SplitMethod>(f));
  methods.push_back(std::make_unique<NewtonMethod>(f, df_dx, d2f_dx2));

  for (auto&& method : methods) {
    for (auto&& segment : segments) {
      const auto [k, x] = method->FindRoot(segment, kEpsilon);
      std::cout << k << ", " << x << "\n";
    }
  }
}
