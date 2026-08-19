// Huber helpers for the robust segmentation cost. Not an NB engine.

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
double huber_mean_cpp(NumericVector y, double k = 1.345, int maxit = 50,
                      double tol = 1e-8) {
  const int n = y.size();
  if (n == 0) return NA_REAL;
  std::vector<double> v(y.begin(), y.end());
  std::nth_element(v.begin(), v.begin() + n / 2, v.end());
  double mu = v[n / 2];
  // MAD
  std::vector<double> dev(n);
  for (int i = 0; i < n; ++i) dev[i] = std::fabs(y[i] - mu);
  std::nth_element(dev.begin(), dev.begin() + n / 2, dev.end());
  double s = 1.4826 * dev[n / 2];
  if (!(s > 0.0) || !std::isfinite(s)) s = 1.0;
  for (int it = 0; it < maxit; ++it) {
    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; ++i) {
      const double r = (y[i] - mu) / s;
      const double w = (std::fabs(r) <= k) ? 1.0 : (k / std::max(std::fabs(r), 1e-12));
      num += w * y[i];
      den += w;
    }
    const double mu_new = (den > 0.0) ? (num / den) : mu;
    if (std::fabs(mu_new - mu) < tol) {
      mu = mu_new;
      break;
    }
    mu = mu_new;
  }
  return mu;
}

// [[Rcpp::export]]
double huber_cost_cpp(NumericVector y, double k = 1.345) {
  const double mu = huber_mean_cpp(y, k);
  if (!std::isfinite(mu)) return NA_REAL;
  double c = 0.0;
  for (int i = 0; i < y.size(); ++i) {
    const double r = y[i] - mu;
    const double ar = std::fabs(r);
    if (ar <= k) c += 0.5 * r * r;
    else c += k * ar - 0.5 * k * k;
  }
  return c;
}
