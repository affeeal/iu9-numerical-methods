#include <cassert>
#include <iomanip>
#include <cmath>
#include <ios>
#include <iostream>
#include <numbers>
#include <vector>

int main() {
  const auto y = [](const double x) -> double { return std::exp(x) + x * x; };
  constexpr auto p = [](const double x) -> double { return -1; };
  constexpr auto q = [](const double x) -> double { return 0; };
  constexpr auto f = [](const double x) -> double { return 2 * (1 - x); };

  constexpr double a = 0, b = 1;
  constexpr double A = 1, B = std::numbers::e + 1;
  constexpr std::size_t n = 10;

  constexpr auto h = (b - a) / n;
  constexpr auto h_sqr = h * h;
  constexpr auto h_half = h / 2;

  std::vector<double> xs;
  xs.reserve(n + 1);

  auto x = a;
  for (std::size_t i = 0; i <= n; ++i) {
    xs.push_back(x);
    x += h;
  }

  constexpr auto O_h = 1e-3;
  constexpr auto D0 = A + O_h;
  constexpr auto D1 = O_h;

  std::vector<double> y0s, y1s;
  y0s.reserve(n + 1);
  y1s.reserve(n + 1);

  y0s.push_back(A);
  y0s.push_back(D0);

  y1s.push_back(0);
  y1s.push_back(D1);

  for (std::size_t i = 1; i < n; ++i) {
    const auto pi = p(xs[i]), qi = q(xs[i]), fi = f(xs[i]);
    y0s.push_back((fi * h_sqr + (2 - qi * h_sqr) * y0s[i] -
                   (1 - pi * h_half) * y0s[i - 1]) /
                  (1 + pi * h_half));
    y1s.push_back(((2 - qi * h_sqr) * y1s[i] - (1 - pi * h_half) * y1s[i - 1]) /
                  (1 + pi * h_half));
  }

  const auto C1 = (B - y0s[n]) / y1s[n];

  std::vector<double> ys;
  ys.reserve(n + 1);

  for (std::size_t i = 0; i <= n; ++i) {
    ys.push_back(y0s[i] + C1 * y1s[i]);
  }

  const auto w = 9;

  std::cout << std::setw(w) << "i";
  for (std::size_t i = 0; i <= n; ++i) {
    std::cout << std::setw(w) << i;
  }
  std::cout << std::endl;

  std::cout << std::setw(w) << "x";
  for (auto&& x : xs) {
    std::cout << std::setw(w) << x;
  }
  std::cout << std::endl;

  std::cout << std::setw(w) << "y_an";
  for (std::size_t i = 0; i <= n; ++i) {
    std::cout << std::setw(w) << std::fixed << y(xs[i]);
  }
  std::cout << std::endl;

  std::cout << std::setw(w) << "y_nm";
  for (auto&& y : ys) {
    std::cout << std::setw(w) << std::fixed << y;
  }
  std::cout << std::endl;

  std::cout << std::setw(w) << "error";
  for (std::size_t i = 0; i <= n; ++i) {
    std::cout << std::setw(w) << std::fixed << std::abs(y(xs[i]) - ys[i]);
  }
  std::cout << std::endl;
}
