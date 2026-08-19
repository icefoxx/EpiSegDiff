#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector interval_overlap_cpp(IntegerVector as, IntegerVector ae,
                                   IntegerVector bs, IntegerVector be) {
  const int n = as.size();
  IntegerVector out(n);
  for (int i = 0; i < n; ++i) {
    const int left = std::max(as[i], bs[i]);
    const int right = std::min(ae[i], be[i]);
    out[i] = std::max(0, right - left);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector reciprocal_overlap_cpp(IntegerVector as, IntegerVector ae,
                                     IntegerVector bs, IntegerVector be) {
  const int n = as.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const int wa = ae[i] - as[i];
    const int wb = be[i] - bs[i];
    const int left = std::max(as[i], bs[i]);
    const int right = std::min(ae[i], be[i]);
    const int inter = std::max(0, right - left);
    if (wa <= 0 || wb <= 0) {
      out[i] = 0.0;
    } else {
      const double ra = static_cast<double>(inter) / static_cast<double>(wa);
      const double rb = static_cast<double>(inter) / static_cast<double>(wb);
      out[i] = std::min(ra, rb);
    }
  }
  return out;
}

// [[Rcpp::export]]
CharacterVector interval_relation_cpp(int as, int ae, int bs, int be) {
  if (as == bs && ae == be) return "equal";
  if (as >= bs && ae <= be) return "i_in_j";
  if (bs >= as && be <= ae) return "j_in_i";
  const int left = std::max(as, bs);
  const int right = std::min(ae, be);
  if (right > left) return "partial";
  return "disjoint";
}
