#include <cmath>
#include <iostream>

namespace {

double F(const double x, const double y) {
  return std::exp(x) + std::pow(x + y, 2);
};

double DfDx(const double x, const double y) {
  return std::exp(x) + 2 * (x + y);
};

double D2fDx2(const double x, const double y) { return std::exp(x) + 2; };

double DfDy(const double x, const double y) { return 2 * (x + y); };

double D2fDy2(const double x, const double y) { return 2; };

double D2fDxDy(const double x, const double y) { return 2; };

}  // namespace

int main() {
  constexpr auto kEpsilon = 1e-3;

  auto x = 1.0;
  auto y = 1.0;

  do {
    const auto df_dx = DfDx(x, y);
    const auto df_dy = DfDy(x, y);

    if (std::max(df_dx, df_dy) <= kEpsilon) {
      break;
    }

    const auto sqr_df_dx = std::pow(df_dx, 2);
    const auto sqr_df_dy = std::pow(df_dy, 2);

    const auto dphi_dt = -sqr_df_dx - sqr_df_dy;
    const auto d2phi_dt2 = D2fDx2(x, y) * sqr_df_dx +
                           2 * D2fDxDy(x, y) * df_dx * df_dy +
                           D2fDy2(x, y) * sqr_df_dy;
    const auto t = -dphi_dt / d2phi_dt2;

    x -= t * df_dx;
    y -= t * df_dy;
  } while (true);
  
  const auto f_min = F(x, y);
  std::cout << "Function: exp(x) + (x + y)^2, (x0, y0) = (1.0, 1.0);\n"
            << "Calculated point of minimum: (" << x << ", " << y << ");\n"
            << "Calculated minimum: " << f_min << ";\n"
            << "Analytical minimum: tends to 0;\n" 
            << "Error: " << f_min << ".\n";
}
