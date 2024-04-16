#include <cassert>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>

static constexpr std::size_t kWidth = 12;

namespace {

void PrintResult(const std::vector<double>& xs, const std::vector<double>& ys1,
                 const std::vector<double>& ys2,
                 const std::vector<double>& errors) {
  const auto n = xs.size();

  assert(n == ys1.size() && n == ys2.size() && n == errors.size());

  std::cout << "Comparsion\n"
            << std::setw(kWidth) << "i" << std::setw(kWidth) << "x_i"
            << std::setw(kWidth) << "actual y_i" << std::setw(kWidth)
            << "spline y_i" << std::setw(kWidth) << "error" << '\n';

  for (std::size_t i = 0; i < n; ++i) {
    std::cout << std::setw(kWidth) << i << std::setw(kWidth) << xs[i]
              << std::setw(kWidth) << ys1[i] << std::setw(kWidth) << ys2[i]
              << std::setw(kWidth) << errors[i] << '\n';
  }
}

void PrintCoeffs(const std::vector<double>& as, const std::vector<double>& bs,
                 const std::vector<double>& cs, const std::vector<double>& ds) {
  const auto n = as.size();

  assert(n == bs.size() && n == cs.size() && n == ds.size());

  std::cout << "Coefficient arrays a, b, c, d\n"
            << std::setw(kWidth) << "i" << std::setw(kWidth) << "a_i"
            << std::setw(kWidth) << "b_i" << std::setw(kWidth) << "c_i"
            << std::setw(kWidth) << "d_i" << '\n';

  for (std::size_t i = 0; i < n; ++i) {
    std::cout << std::setw(kWidth) << i << std::setw(kWidth) << as[i]
              << std::setw(kWidth) << bs[i] << std::setw(kWidth) << cs[i]
              << std::setw(kWidth) << ds[i] << '\n';
  }
}

std::pair<std::vector<double>, std::vector<double>> Tabulate(
    const std::function<double(const double)>& f, const std::size_t n,
    const double a, const double h) {
  std::vector<double> xs, ys;
  xs.reserve(n + 1);
  ys.reserve(n + 1);

  xs.push_back(a);
  ys.push_back(f(a));

  for (std::size_t i = 1; i <= n; ++i) {
    xs.push_back(xs[i - 1] + h);
    ys.push_back(f(xs[i]));
  }

  return {std::move(xs), std::move(ys)};
}

std::vector<std::vector<double>> GetMatrix(const std::size_t n) {
  assert(n >= 3);

  std::vector<std::vector<double>> m(n - 1, std::vector<double>(n - 1));

  m[0][0] = 4;
  m[0][1] = 1;

  for (std::size_t i = 1, end = n - 2; i < end; ++i) {
    m[i][i - 1] = 1;
    m[i][i] = 4;
    m[i][i + 1] = 1;
  }

  m[n - 2][n - 3] = 1;
  m[n - 2][n - 2] = 4;

  return m;
}

std::vector<double> GetFreeCoeffs(const std::vector<double>& ys,
                                  const std::size_t n, const std::size_t h) {
  std::vector<double> fc;
  fc.reserve(n - 1);

  for (std::size_t i = 1; i < n; ++i) {
    fc.push_back((ys[i + 1] - 2 * ys[i] + ys[i - 1]) / h / h);
  }

  return fc;
}

bool CheckDiagonalPredominanceConditions(
    const std::vector<std::vector<double>>& matrix) {
  const auto n = matrix.size();

  if (std::abs(matrix[0][1] / matrix[0][0]) > 1) {
    return false;
  }

  if (std::abs(matrix[n - 1][n - 2] / matrix[n - 1][n - 1]) > 1) {
    return false;
  }

  for (auto i = 1; i < n; i++) {
    if (std::abs(matrix[i][i]) <
        std::abs(matrix[i][i - 1]) + std::abs(matrix[i][i + 1])) {
      return false;
    }
  }

  return true;
}

std::vector<double> CalculateSolution(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& free_coeffs) {
  assert(CheckDiagonalPredominanceConditions(matrix));

  const auto n = free_coeffs.size();

  std::vector<double> alphas(n - 1);
  std::vector<double> betas(n - 1);

  assert(matrix[0][0] != 0);

  alphas[0] = -matrix[0][1] / matrix[0][0];
  betas[0] = free_coeffs[0] / matrix[0][0];

  for (auto i = 1; i < n - 1; i++) {
    alphas[i] =
        -matrix[i][i + 1] / (matrix[i][i - 1] * alphas[i - 1] + matrix[i][i]);
    betas[i] = (free_coeffs[i] - matrix[i][i - 1] * betas[i - 1]) /
               (matrix[i][i - 1] * alphas[i - 1] + matrix[i][i]);
  }

  std::vector<double> result(n);
  result[n - 1] = (free_coeffs[n - 1] - matrix[n - 1][n - 2] * betas[n - 2]) /
                  (matrix[n - 1][n - 2] * alphas[n - 2] + matrix[n - 1][n - 1]);

  for (int i = n - 2; i >= 0; i--) {
    result[i] = alphas[i] * result[i + 1] + betas[i];
  }

  return result;
}

}  // namespace

int main() {
  constexpr std::size_t n = 32;

  constexpr auto a = 0.0;
  constexpr auto b = std::numbers::pi;
  constexpr auto h = (b - a) / n;
  const auto f = [](const double x) { return 2 * x * std::cos(x / 2); };

  const auto [xs, ys] = Tabulate(f, n, a, h);

  const auto m = GetMatrix(n);
  const auto fc = GetFreeCoeffs(ys, n, h);
  const auto s = CalculateSolution(m, fc);  // cs[1..n-1]

  std::vector<double> as, bs, cs, ds;
  as.reserve(n);
  bs.reserve(n);
  cs.reserve(n);
  ds.reserve(n);

  for (std::size_t i = 0; i < n; ++i) {
    as.push_back(ys[i]);
  }

  cs.push_back(0);
  for (std::size_t i = 0, end = n - 1; i < end; ++i) {
    cs.push_back(ys[i]);
  }

  for (std::size_t i = 0, end = n - 1; i < end; ++i) {
    bs.push_back((ys[i + 1] - ys[i]) / 3 * (cs[i + 1] + 2 * cs[i]));
  }
  bs.push_back((ys[n] - ys[n - 1]) * 2 * cs[n - 1] / 3);

  for (std::size_t i = 0, end = n - 1; i < end; ++i) {
    ds.push_back((cs[i + 1] - cs[i]) / 3 / h);
  }
  ds.push_back(-cs[n] / 3 / h);

  PrintCoeffs(as, bs, cs, ds);

  std::vector<double> test_xs, f_ys, spline_ys, errors;
  test_xs.reserve(n);
  f_ys.reserve(n);
  spline_ys.reserve(n);
  errors.reserve(n);

  for (std::size_t i = 0; i < n; i++) {
    test_xs.push_back(a + (i + 0.5) * h);
    f_ys.push_back(f(test_xs[i]));

    const auto x = test_xs[i] - xs[i];
    spline_ys.push_back(as[i] + bs[i] * x + cs[i] * std::pow(x, 2) +
                        ds[i] * std::pow(x, 3));

    errors.push_back(std::abs(f_ys[i] - spline_ys[i]));
  }

  PrintResult(test_xs, f_ys, spline_ys, errors);
}
