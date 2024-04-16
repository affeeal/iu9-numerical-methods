#include <cassert>
#include <cmath>
#include <vector>

namespace {

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
  const std::vector<std::vector<double>> matrix{
      {4, 1, 0, 0},
      {1, 4, 1, 0},
      {0, 1, 4, 1},
      {0, 0, 1, 4},
  };
  const std::vector<double> free_coeffs{5, 6, 6, 5};
  const std::vector<double> expected_result{1, 1, 1, 1};

  const auto result = CalculateSolution(matrix, free_coeffs);
  assert(result == expected_result);
}
