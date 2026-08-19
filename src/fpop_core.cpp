// Gaussian change-in-mean PELT (Killick et al.). Used by the MVP
// segmentation engine. This is NOT a Negative-Binomial FPOP engine.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <utility>

using namespace Rcpp;

// [[Rcpp::export]]
List pelt_gaussian_cpp(NumericVector y, double penalty, int min_width) {
  const int n = y.size();
  IntegerVector empty(0);
  if (n == 0) {
    return List::create(Named("changepoints") = empty, Named("cost") = 0.0);
  }
  if (min_width < 1) min_width = 1;

  std::vector<double> S(n + 1, 0.0), Q(n + 1, 0.0);
  for (int i = 1; i <= n; ++i) {
    S[i] = S[i - 1] + y[i - 1];
    Q[i] = Q[i - 1] + y[i - 1] * y[i - 1];
  }
  auto cost = [&](int s, int t) {
    const int m = t - s + 1;
    const double ss = S[t] - S[s - 1];
    const double qq = Q[t] - Q[s - 1];
    return qq - (ss * ss) / static_cast<double>(m);
  };

  const double INF = std::numeric_limits<double>::infinity();
  std::vector<double> F(n + 1, INF);
  std::vector<int> last(n + 1, 0);
  F[0] = -penalty;
  std::vector<int> R(n + 2, 0);
  std::vector<int> R2(n + 2, 0);
  int rlen = 1;
  R[0] = 0;
  int max_r = 1;

  for (int t = 1; t <= n; ++t) {
    if (t % 20000 == 0) Rcpp::checkUserInterrupt();
    double best = INF;
    int best_tau = 0;
    for (int i = 0; i < rlen; ++i) {
      const int tau = R[i];
      const int len = t - tau;
      if (len < min_width && t < n) continue;
      if (len < 1) continue;
      const double val = F[tau] + cost(tau + 1, t) + penalty;
      if (val < best) {
        best = val;
        best_tau = tau;
      }
    }
    if (!std::isfinite(best)) {
      best = F[0] + cost(1, t) + penalty;
      best_tau = 0;
    }
    F[t] = best;
    last[t] = best_tau;
    if (rlen > max_r) max_r = rlen;

    int r2 = 0;
    for (int i = 0; i < rlen; ++i) {
      const int tau = R[i];
      if (t - tau < 1) continue;
      if (F[tau] + cost(tau + 1, t) <= F[t]) R2[r2++] = tau;
    }
    R2[r2++] = t;
    R.swap(R2);
    rlen = r2;
  }

  std::vector<int> cps;
  int curr = n;
  int guard = 0;
  while (curr > 0 && guard++ <= n) {
    const int tau = last[curr];
    if (tau > 0) cps.push_back(tau);
    curr = tau;
  }
  std::reverse(cps.begin(), cps.end());
  IntegerVector cp(cps.begin(), cps.end());
  return List::create(Named("changepoints") = cp, Named("cost") = F[n],
                      Named("max_candidates") = max_r);
}

// 1-D Gaussian-mean FPOP (Maidstone / Rigaill / Fearnhead).
// Same RSS + penalty cost as pelt_gaussian_cpp. This is NOT NB-FPOP.

namespace {

struct FIval {
  double lo, hi;
};

struct FPiece {
  int tau;
  double a, b, c;
  double lo, hi;
};

const double F_INF = std::numeric_limits<double>::infinity();
const double F_EPS = 1e-12;

std::vector<FIval> where_positive(double A, double B, double C) {
  std::vector<FIval> out;
  if (std::abs(A) < F_EPS) {
    if (std::abs(B) < F_EPS) {
      if (C > 0) out.push_back({-F_INF, F_INF});
      return out;
    }
    const double root = -C / B;
    if (B > 0) out.push_back({root, F_INF});
    else out.push_back({-F_INF, root});
    return out;
  }
  double disc = B * B - 4.0 * A * C;
  if (disc < 0) {
    if (A > 0) out.push_back({-F_INF, F_INF});
    return out;
  }
  const double s = std::sqrt(std::max(disc, 0.0));
  double r1 = (-B - s) / (2.0 * A);
  double r2 = (-B + s) / (2.0 * A);
  if (r1 > r2) std::swap(r1, r2);
  if (A > 0) {
    out.push_back({-F_INF, r1});
    out.push_back({r2, F_INF});
  } else {
    out.push_back({r1, r2});
  }
  return out;
}

std::vector<FIval> clip_ivals(const std::vector<FIval>& xs, double lo, double hi) {
  std::vector<FIval> out;
  for (size_t i = 0; i < xs.size(); ++i) {
    const double L = std::max(xs[i].lo, lo);
    const double R = std::min(xs[i].hi, hi);
    if (R > L + 1e-15) out.push_back({L, R});
  }
  return out;
}

std::vector<FIval> complement_inside(double lo, double hi, std::vector<FIval> cut) {
  std::sort(cut.begin(), cut.end(), [](const FIval& u, const FIval& v) {
    return u.lo < v.lo;
  });
  std::vector<FIval> out;
  double cur = lo;
  for (size_t i = 0; i < cut.size(); ++i) {
    if (cut[i].lo > cur + 1e-15) out.push_back({cur, cut[i].lo});
    if (cut[i].hi > cur) cur = cut[i].hi;
  }
  if (hi > cur + 1e-15) out.push_back({cur, hi});
  return out;
}

void emit_piece(std::vector<FPiece>& out, const FPiece& proto, double lo, double hi) {
  if (!(hi > lo + 1e-15)) return;
  FPiece q = proto;
  q.lo = lo;
  q.hi = hi;
  if (!out.empty() && out.back().tau == q.tau &&
      std::abs(out.back().hi - lo) < 1e-10) {
    out.back().hi = hi;
    return;
  }
  out.push_back(q);
}

std::vector<FPiece> insert_piece(const std::vector<FPiece>& env, FPiece neu) {
  if (env.empty()) {
    neu.lo = -F_INF;
    neu.hi = F_INF;
    return std::vector<FPiece>(1, neu);
  }
  std::vector<FPiece> out;
  out.reserve(env.size() + 4);
  for (size_t i = 0; i < env.size(); ++i) {
    const FPiece& e = env[i];
    const double A = e.a - neu.a;
    const double B = e.b - neu.b;
    const double C = e.c - neu.c;
    std::vector<FIval> pwin = clip_ivals(where_positive(A, B, C), e.lo, e.hi);
    std::vector<FIval> ewin = complement_inside(e.lo, e.hi, pwin);
    std::vector<std::pair<double, int> > ev;
    for (size_t j = 0; j < ewin.size(); ++j) {
      ev.push_back(std::make_pair(ewin[j].lo, 0));
      ev.push_back(std::make_pair(ewin[j].hi, 0));
    }
    for (size_t j = 0; j < pwin.size(); ++j) {
      ev.push_back(std::make_pair(pwin[j].lo, 1));
      ev.push_back(std::make_pair(pwin[j].hi, 1));
    }
    std::sort(ev.begin(), ev.end());
    // Walk unique endpoints in order and assign ownership by midpoint test
    std::vector<double> xs;
    xs.push_back(e.lo);
    for (size_t j = 0; j < pwin.size(); ++j) {
      xs.push_back(pwin[j].lo);
      xs.push_back(pwin[j].hi);
    }
    xs.push_back(e.hi);
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    for (size_t j = 1; j < xs.size(); ++j) {
      const double L = xs[j - 1];
      const double R = xs[j];
      if (!(R > L + 1e-15)) continue;
      double mid = 0.5 * (L + R);
      if (!std::isfinite(mid)) {
        mid = std::isfinite(L) ? L + 1.0 : R - 1.0;
      }
      const double pv = neu.a * mid * mid + neu.b * mid + neu.c;
      const double evv = e.a * mid * mid + e.b * mid + e.c;
      if (pv < evv) emit_piece(out, neu, L, R);
      else emit_piece(out, e, L, R);
    }
  }
  return out;
}

double piece_min(const FPiece& p) {
  if (p.a > F_EPS) return p.c - (p.b * p.b) / (4.0 * p.a);
  return p.c;
}

}  // namespace

// [[Rcpp::export]]
List fpop_gaussian_cpp(NumericVector y, double penalty, int min_width) {
  const int n = y.size();
  IntegerVector empty(0);
  if (n == 0) {
    return List::create(Named("changepoints") = empty, Named("cost") = 0.0,
                        Named("max_candidates") = 0);
  }
  if (min_width < 1) min_width = 1;

  std::vector<double> F(n + 1, F_INF);
  std::vector<int> last(n + 1, 0);
  F[0] = -penalty;

  FPiece first;
  first.tau = 0;
  first.a = 0;
  first.b = 0;
  first.c = 0;
  first.lo = -F_INF;
  first.hi = F_INF;
  std::vector<FPiece> env(1, first);
  int max_k = 1;

  for (int t = 1; t <= n; ++t) {
    if (t % 100000 == 0) Rcpp::checkUserInterrupt();
    const double yt = y[t - 1];
    const double yt2 = yt * yt;
    for (size_t i = 0; i < env.size(); ++i) {
      env[i].a += 1.0;
      env[i].b -= 2.0 * yt;
      env[i].c += yt2;
    }
    double best = F_INF;
    int best_tau = 0;
    for (size_t i = 0; i < env.size(); ++i) {
      const int tau = env[i].tau;
      const int len = t - tau;
      if (len < min_width && t < n) continue;
      if (len < 1) continue;
      const double val = piece_min(env[i]);
      if (val < best) {
        best = val;
        best_tau = tau;
      }
    }
    if (!std::isfinite(best)) {
      best = 0.0;
      best_tau = 0;
    }
    F[t] = best;
    last[t] = best_tau;

    FPiece neu;
    neu.tau = t;
    neu.a = 0;
    neu.b = 0;
    neu.c = F[t] + penalty;
    neu.lo = -F_INF;
    neu.hi = F_INF;
    env = insert_piece(env, neu);
    if (static_cast<int>(env.size()) > max_k) max_k = static_cast<int>(env.size());
  }

  std::vector<int> cps;
  int curr = n;
  int guard = 0;
  while (curr > 0 && guard++ <= n) {
    const int tau = last[curr];
    if (tau > 0) cps.push_back(tau);
    curr = tau;
  }
  std::reverse(cps.begin(), cps.end());
  IntegerVector cp(cps.begin(), cps.end());
  return List::create(Named("changepoints") = cp, Named("cost") = F[n],
                      Named("max_candidates") = max_k);
}
