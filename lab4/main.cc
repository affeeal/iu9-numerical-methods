// clang-format off
#include <sprout/math/exp.hpp>
// clang-format on

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <vector>

namespace {

bool CheckDiagonalPredominanceConditions(
    const std::vector<std::vector<double>>& mat) {
  const auto n = mat.size();

  if (std::abs(mat[0][1] / mat[0][0]) > 1) {
    return false;
  }

  if (std::abs(mat[n - 1][n - 2] / mat[n - 1][n - 1]) > 1) {
    return false;
  }

  for (std::size_t i = 1; i < n; i++) {
    if (std::abs(mat[i][i]) <
        std::abs(mat[i][i - 1]) + std::abs(mat[i][i + 1])) {
      return false;
    }
  }

  return true;
}

std::vector<double> RunThrough(const std::vector<std::vector<double>>& mat,
                               const std::vector<double>& cfs) {
  if (!CheckDiagonalPredominanceConditions(mat)) {
    throw std::runtime_error("Run-through predominance conditions failed");
  }

  const auto n = cfs.size();

  std::vector<double> as(n - 1);
  std::vector<double> bs(n - 1);

  as[0] = -mat[0][1] / mat[0][0];
  bs[0] = cfs[0] / mat[0][0];

  for (std::size_t i = 1, end = n - 1; i < end; i++) {
    as[i] = -mat[i][i + 1] / (mat[i][i - 1] * as[i - 1] + mat[i][i]);
    bs[i] = (cfs[i] - mat[i][i - 1] * bs[i - 1]) /
            (mat[i][i - 1] * as[i - 1] + mat[i][i]);
  }

  std::vector<double> res(n);
  res[n - 1] = (cfs[n - 1] - mat[n - 1][n - 2] * bs[n - 2]) /
               (mat[n - 1][n - 2] * as[n - 2] + mat[n - 1][n - 1]);

  for (int i = n - 2; i >= 0; i--) {
    res[i] = as[i] * res[i + 1] + bs[i];
  }

  return res;
}

void PrintResult(const std::vector<double>& xs, const std::vector<double>& ys,
                 const std::vector<double>& errs,
                 const std::function<double(const double)>& y) {
  static constexpr std::size_t kWidth = 10;

  const auto n = xs.size();
  assert(ys.size() == n);

  std::cout << std::setw(kWidth) << "i";
  for (std::size_t i = 0; i < n; i++) {
    std::cout << std::setw(kWidth) << i;
  }
  std::cout << "\n";

  std::cout << std::setw(kWidth) << "x_i";
  for (std::size_t i = 0; i < n; i++) {
    std::cout << std::setw(kWidth) << std::fixed << xs[i];
  }
  std::cout << "\n";

  std::cout << std::setw(kWidth) << "y(x_i)";
  for (std::size_t i = 0; i < n; i++) {
    std::cout << std::setw(kWidth) << std::fixed << y(xs[i]);
  }
  std::cout << "\n";

  std::cout << std::setw(kWidth) << "y_i";
  for (std::size_t i = 0; i < n; i++) {
    std::cout << std::setw(kWidth) << std::fixed << ys[i];
  }
  std::cout << "\n";

  std::cout << std::setw(kWidth) << "error";
  for (std::size_t i = 0; i < n; i++) {
    std::cout << std::setw(kWidth) << std::fixed << errs[i];
  }
  std::cout << "\n";

  const auto max_err = std::max_element(errs.begin(), errs.end());
  std::cout << "max error = " << *max_err << "\n";
}

}  // namespace

int main() {
  constexpr auto y = [](const double x) {
    return x * x + 9 * sprout::exp(x) - 8;
  };

  constexpr double a = 1;
  constexpr auto b = y(1);
  constexpr auto p = [](const double x) { return -1; };
  constexpr auto q = [](const double x) { return 0; };
  constexpr auto f = [](const double x) { return 2 * (1 - x); };

  constexpr std::size_t n = 10;
  constexpr auto h = 1 / static_cast<double>(n);

  // Fill xs, ps, qs, fs
  std::vector<double> xs, ps, qs, fs;

  xs.reserve(n + 1);
  ps.reserve(n + 1);
  qs.reserve(n + 1);
  fs.reserve(n + 1);

  double x = 0;
  for (std::size_t i = 0; i <= n; ++i) {
    xs.push_back(x);
    ps.push_back(p(x));
    qs.push_back(q(x));
    fs.push_back(f(x));

    x += h;
  }

  // Fill matrix, free coefficients
  std::vector<std::vector<double>> mat(n - 1, std::vector<double>(n - 1));

  std::vector<double> cfs;
  cfs.reserve(n - 1);

  const auto h_sqr = h * h;
  const auto h_hlf = h / 2;

  mat[0][0] = h_sqr * qs[1] - 2;
  mat[0][1] = 1 + h_hlf * ps[1];

  cfs.push_back(h_sqr * fs[1] - a * (1 - h_hlf * ps[1]));

  for (std::size_t i = 1, end = n - 2; i < end; ++i) {
    mat[i][i - 1] = 1 - h_hlf * ps[i + 1];
    mat[i][i] = h_sqr * qs[i + 1] - 2;
    mat[i][i + 1] = 1 + h_hlf * ps[i + 1];

    cfs.push_back(h_sqr * fs[i + 1]);
  }

  mat[n - 2][n - 3] = 1 - h_hlf * ps[n - 1];
  mat[n - 2][n - 2] = h_sqr * qs[n - 1] - 2;

  cfs.push_back(h_sqr * fs[n - 1] - b * (1 + h_hlf * ps[n - 1]));

  // Run-through
  const auto sol = RunThrough(mat, cfs);

  std::vector<double> ys;
  ys.reserve(n + 1);

  ys.push_back(a);

  for (std::size_t i = 0, end = n - 1; i < end; ++i) {
    ys.push_back(sol[i]);
  }

  ys.push_back(b);

  // Error
  std::vector<double> errs;
  errs.reserve(n + 1);

  for (std::size_t i = 0; i <= n; ++i) {
    errs.push_back(std::abs(y(xs[i]) - ys[i]));
  }

  PrintResult(xs, ys, errs, y);
}
