// cellMMCD cell step in C++ (block 4). Same semantics as cm_wstep in R/cellmmcd.R: sequential over positions c,
// for each position the kept observations' flags are reset to Wfix, Delta_ic - lambda_kl is computed with the
// current flagged set M_i of observation i (S = M_i + {c}, Schur-form conditional), the largest positive values are
// flagged, at most budget of them; then (shadow) the trimmed observations' cells are re-decided sequentially.
// Tf: p x q x n array Or (X_i - M) Oc; W, Wfix: logical p x q x n (vec index c = k + l p, 0-based here).
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

static double delta_cell(const double* Tf, const arma::mat& Or, const arma::mat& Oc, const std::vector<int>& M,
                         int c, int p) {
  std::vector<int> S; S.reserve(M.size() + 1);
  for (size_t a = 0; a < M.size(); ++a) if (M[a] != c) S.push_back(M[a]);
  S.push_back(c);
  int jc = S.size() - 1;
  if (jc == 0) {
    int k = c % p, l = c / p; double s2 = 1.0 / (Or(k, k) * Oc(l, l));
    return std::log(2.0 * M_PI * s2) + Tf[c] * Tf[c] * s2;
  }
  arma::mat A(jc + 1, jc + 1); arma::vec t(jc + 1);
  for (int a = 0; a <= jc; ++a) {
    int ka = S[a] % p, la = S[a] / p; t(a) = Tf[S[a]];
    for (int b = 0; b <= jc; ++b) A(a, b) = Or(ka, S[b] % p) * Oc(la, S[b] / p);
  }
  arma::mat inv = arma::inv_sympd(A); arma::vec g = inv * t; double s2 = inv(jc, jc);
  return std::log(2.0 * M_PI * s2) + g(jc) * g(jc) / s2;
}

// [[Rcpp::export]]
LogicalVector cm_wstep_cpp(NumericVector Tf, NumericMatrix Or_, NumericMatrix Oc_, LogicalVector W, LogicalVector Wfix,
                           NumericMatrix lambda, int h, LogicalVector v, bool shadow) {
  int p = Or_.nrow(), q = Oc_.nrow(), pq = p * q, n = v.size();
  arma::mat Or(Or_.begin(), p, p, false), Oc(Oc_.begin(), q, q, false);
  LogicalVector Wout = clone(W);
  std::vector<int> kept; for (int i = 0; i < n; ++i) if (v[i]) kept.push_back(i);
  int nv = kept.size();
  std::vector<int> M; std::vector<double> dec; std::vector<int> free_;
  for (int c = 0; c < pq; ++c) {
    int k = c % p, l = c / p;
    free_.clear();
    for (size_t a = 0; a < kept.size(); ++a) { int i = kept[a]; Wout[c + i * pq] = Wfix[c + i * pq]; if (Wfix[c + i * pq]) free_.push_back(i); }
    int budget = nv - h - (nv - (int) free_.size());
    if (budget <= 0) continue;
    dec.assign(free_.size(), 0.0);
    for (size_t a = 0; a < free_.size(); ++a) {
      int i = free_[a]; M.clear();
      for (int j = 0; j < pq; ++j) if (!Wout[j + i * pq]) M.push_back(j);
      dec[a] = delta_cell(&Tf[i * pq], Or, Oc, M, c, p) - lambda(k, l);
    }
    std::vector<int> pos; for (size_t a = 0; a < free_.size(); ++a) if (dec[a] > 0) pos.push_back(a);
    if ((int) pos.size() > budget) {
      std::stable_sort(pos.begin(), pos.end(), [&](int a, int b) { return dec[a] > dec[b]; });
      pos.resize(budget);
    }
    for (size_t a = 0; a < pos.size(); ++a) Wout[c + free_[pos[a]] * pq] = false;
  }
  if (shadow) for (int i = 0; i < n; ++i) {
    if (v[i]) continue;
    for (int j = 0; j < pq; ++j) Wout[j + i * pq] = Wout[j + i * pq] && Wfix[j + i * pq];
    for (int c = 0; c < pq; ++c) {
      if (!Wfix[c + i * pq]) continue;
      M.clear(); for (int j = 0; j < pq; ++j) if (!Wout[j + i * pq]) M.push_back(j);
      Wout[c + i * pq] = delta_cell(&Tf[i * pq], Or, Oc, M, c, p) - lambda(c % p, c / p) <= 0;
    }
  }
  Wout.attr("dim") = W.attr("dim");
  return Wout;
}

// Reweighting flags (cm_reweight): for every kept observation and every cell c, the squared standardized residual
// of c given the kept cells of the raw fit other than c (flagged set M_i of the raw fit); flag iff r2 > crw or forced.
// [[Rcpp::export]]
LogicalVector cm_reweight_cpp(NumericVector Tf, NumericMatrix Or_, NumericMatrix Oc_, LogicalVector W, LogicalVector Wfix,
                              double crw, LogicalVector v) {
  int p = Or_.nrow(), q = Oc_.nrow(), pq = p * q, n = v.size();
  arma::mat Or(Or_.begin(), p, p, false), Oc(Oc_.begin(), q, q, false);
  LogicalVector Wout = clone(W);
  std::vector<int> M, S;
  for (int i = 0; i < n; ++i) {
    if (!v[i]) continue;
    M.clear(); for (int j = 0; j < pq; ++j) if (!W[j + i * pq]) M.push_back(j);
    const double* T = &Tf[i * pq];
    for (int c = 0; c < pq; ++c) {
      S.clear(); for (size_t a = 0; a < M.size(); ++a) if (M[a] != c) S.push_back(M[a]);
      S.push_back(c); int jc = S.size() - 1; double r2;
      if (jc == 0) { int k = c % p, l = c / p; r2 = T[c] * T[c] / (Or(k, k) * Oc(l, l)); }
      else {
        arma::mat A(jc + 1, jc + 1); arma::vec t(jc + 1);
        for (int a = 0; a <= jc; ++a) { int ka = S[a] % p, la = S[a] / p; t(a) = T[S[a]];
          for (int b = 0; b <= jc; ++b) A(a, b) = Or(ka, S[b] % p) * Oc(la, S[b] / p); }
        arma::mat inv = arma::inv_sympd(A); arma::vec g = inv * t; r2 = g(jc) * g(jc) / inv(jc, jc);
      }
      Wout[c + i * pq] = Wfix[c + i * pq] && (r2 <= crw);
    }
  }
  Wout.attr("dim") = W.attr("dim");
  return Wout;
}
