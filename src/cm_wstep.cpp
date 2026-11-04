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
LogicalVector cm_wstep_cpp0(NumericVector Tf, NumericMatrix Or_, NumericMatrix Oc_, LogicalVector W, LogicalVector Wfix,
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

// ---- Cached-inverse version (block 6) -------------------------------------------------------------------
// Per observation i keep M_i (flagged set), inv_i = Omega_{M_i M_i}^{-1} and u_i = inv_i t_{M_i}. The delta of a
// candidate c not in M_i is then Schur: s2 = 1 / (Omega_cc - a' inv_i a), a = Omega_{M_i c}, g_c = s2 (t_c - a' u_i),
// O(m^2) instead of O(m^3). Flagging c borders inv_i (O(m^2)); un-flagging c downdates it (O(m^2)). The inverse is
// recomputed from scratch every REFRESH updates to keep the rounding error bounded.
struct Cache {
  std::vector<int> M; arma::mat inv; arma::vec u; int upd = 0;
  void rebuild(const double* Tf, const arma::mat& Or, const arma::mat& Oc, int p) {
    int m = M.size(); inv.set_size(m, m); u.set_size(m); upd = 0;
    if (m == 0) return;
    arma::mat A(m, m); arma::vec t(m);
    for (int a = 0; a < m; ++a) { int ka = M[a] % p, la = M[a] / p; t(a) = Tf[M[a]];
      for (int b = 0; b < m; ++b) A(a, b) = Or(ka, M[b] % p) * Oc(la, M[b] / p); }
    inv = arma::inv_sympd(A); u = inv * t;
  }
  // delta for candidate c (not in M); fills a, w = inv a, s2 for a possible add
  double delta(const double* Tf, const arma::mat& Or, const arma::mat& Oc, int c, int p, arma::vec& a, arma::vec& w, double& s2) const {
    int m = M.size(), k = c % p, l = c / p; a.set_size(m);
    for (int b = 0; b < m; ++b) a(b) = Or(k, M[b] % p) * Oc(l, M[b] / p);
    double tc = Tf[c], au = 0.0, schur = Or(k, k) * Oc(l, l);
    if (m) { w = inv * a; schur -= arma::dot(a, w); au = arma::dot(a, u); }
    s2 = 1.0 / schur; double gc = s2 * (tc - au);
    return std::log(2.0 * M_PI * s2) + gc * gc / s2;
  }
  void add(const double* Tf, const arma::mat& Or, const arma::mat& Oc, int c, int p, const arma::vec& w, double s2) {
    int m = M.size(); M.push_back(c);
    if (++upd > 64) { rebuild(Tf, Or, Oc, p); return; }
    arma::mat ni(m + 1, m + 1);
    if (m) { ni.submat(0, 0, m - 1, m - 1) = inv + s2 * w * w.t(); ni.submat(0, m, m - 1, m) = -s2 * w; ni.submat(m, 0, m, m - 1) = -s2 * w.t(); }
    ni(m, m) = s2; inv = ni;
    arma::vec t(m + 1); for (int b = 0; b <= m; ++b) t(b) = Tf[M[b]]; u = inv * t;
  }
  void remove(const double* Tf, const arma::mat& Or, const arma::mat& Oc, int c, int p) {
    int m = M.size(), j = -1; for (int b = 0; b < m; ++b) if (M[b] == c) { j = b; break; }
    if (j < 0) return;
    M.erase(M.begin() + j);
    if (++upd > 64 || m == 1) { rebuild(Tf, Or, Oc, p); return; }
    arma::uvec keep(m - 1); for (int b = 0, r = 0; b < m; ++b) if (b != j) keep(r++) = b;
    arma::vec col = inv.col(j); double d = inv(j, j);
    arma::mat ni = inv.submat(keep, keep) - (col(keep) * col(keep).t()) / d; inv = ni;
    arma::vec t(m - 1); for (int b = 0; b < m - 1; ++b) t(b) = Tf[M[b]]; u = inv * t;
  }
};

// [[Rcpp::export]]
LogicalVector cm_wstep_cpp(NumericVector Tf, NumericMatrix Or_, NumericMatrix Oc_, LogicalVector W, LogicalVector Wfix,
                           NumericMatrix lambda, int h, LogicalVector v, bool shadow) {
  int p = Or_.nrow(), q = Oc_.nrow(), pq = p * q, n = v.size();
  arma::mat Or(Or_.begin(), p, p, false), Oc(Oc_.begin(), q, q, false);
  LogicalVector Wout = clone(W);
  std::vector<int> kept; for (int i = 0; i < n; ++i) if (v[i]) kept.push_back(i);
  int nv = kept.size();
  std::vector<Cache> cache(n);
  for (int i = 0; i < n; ++i) { if (!v[i] && !shadow) continue;
    for (int j = 0; j < pq; ++j) if (!Wout[j + i * pq]) cache[i].M.push_back(j);
    cache[i].rebuild(&Tf[i * pq], Or, Oc, p); }
  std::vector<double> dec; std::vector<int> free_; std::vector<arma::vec> ws; std::vector<double> s2s; arma::vec a, w; double s2;
  for (int c = 0; c < pq; ++c) {
    int k = c % p, l = c / p;
    free_.clear();
    for (size_t a_ = 0; a_ < kept.size(); ++a_) { int i = kept[a_];
      if (Wfix[c + i * pq]) { if (!Wout[c + i * pq]) cache[i].remove(&Tf[i * pq], Or, Oc, c, p); Wout[c + i * pq] = true; free_.push_back(i); } }
    int budget = nv - h - (nv - (int) free_.size());
    if (budget <= 0) continue;
    dec.assign(free_.size(), 0.0); ws.resize(free_.size()); s2s.assign(free_.size(), 0.0);
    for (size_t a_ = 0; a_ < free_.size(); ++a_) { int i = free_[a_];
      dec[a_] = cache[i].delta(&Tf[i * pq], Or, Oc, c, p, a, ws[a_], s2s[a_]) - lambda(k, l); }
    std::vector<int> pos; for (size_t a_ = 0; a_ < free_.size(); ++a_) if (dec[a_] > 0) pos.push_back(a_);
    if ((int) pos.size() > budget) {
      std::stable_sort(pos.begin(), pos.end(), [&](int x, int y) { return dec[x] > dec[y]; });
      pos.resize(budget);
    }
    for (size_t a_ = 0; a_ < pos.size(); ++a_) { int i = free_[pos[a_]]; Wout[c + i * pq] = false; cache[i].add(&Tf[i * pq], Or, Oc, c, p, ws[pos[a_]], s2s[pos[a_]]); }
  }
  if (shadow) for (int i = 0; i < n; ++i) {
    if (v[i]) continue;
    Cache& ch = cache[i];
    for (int j = 0; j < pq; ++j) if (!Wfix[j + i * pq] && Wout[j + i * pq]) { Wout[j + i * pq] = false; ch.M.push_back(j); }
    ch.rebuild(&Tf[i * pq], Or, Oc, p);
    for (int c = 0; c < pq; ++c) {
      if (!Wfix[c + i * pq]) continue;
      if (!Wout[c + i * pq]) { ch.remove(&Tf[i * pq], Or, Oc, c, p); Wout[c + i * pq] = true; }
      double d = ch.delta(&Tf[i * pq], Or, Oc, c, p, a, w, s2) - lambda(c % p, c / p);
      if (d > 0) { Wout[c + i * pq] = false; ch.add(&Tf[i * pq], Or, Oc, c, p, w, s2); }
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
