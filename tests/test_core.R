# Unit tests for the core functions. Run: Rscript tests/test_core.R
source("R/load.R")
set.seed(1)
ok <- function(msg, cond) { cat(sprintf("[%s] %s\n", if (cond) "PASS" else "FAIL", msg)); invisible(cond) }
res <- logical()

## 1. Sampler: empirical cov of vec(X) -> Sc %x% Sr
p <- 3; q <- 4; S <- true_sigma(p, q, "full")
X <- rmatnorm(20000, Sr = S$Sr, Sc = S$Sc)
V <- cov(t(matrix(X, p * q)))
err <- max(abs(V - S$Sc %x% S$Sr))
res["sampler"] <- ok(sprintf("sampler: max |cov(vec X) - Sc x Sr| = %.3f (n = 20000)", err), err < 0.05)
if (has_robustmatrix) {
  X2 <- robustmatrix::rmatnorm(20000, matrix(0, p, q), S$Sr, S$Sc)
  err2 <- max(abs(cov(t(matrix(X2, p * q))) - S$Sc %x% S$Sr))
  res["sampler_rm"] <- ok(sprintf("robustmatrix::rmatnorm has the same convention (err %.3f)", err2), err2 < 0.05)
}

## 2. flipflop vs robustmatrix::mmle (both estimate M here)
X <- rmatnorm(30, Sr = S$Sr, Sc = S$Sc)
ff <- flipflop(X, center = TRUE)
if (has_robustmatrix) {
  # mmle stops at a loose internal tolerance (~1e-4), so compare at its iteration count
  ml <- robustmatrix::mmle(X, max_iter = 1000L)
  ffk <- flipflop(X, center = TRUE, maxit = ml$iterations, tol = 0)
  e <- max(abs(ffk$Sc %x% ffk$Sr - ml$cov_col %x% ml$cov_row))
  res["mmle"] <- ok(sprintf("flipflop == robustmatrix::mmle after %d iterations (max diff %.1e)",
                            ml$iterations, e), e < 1e-10)
}

## 3. root: A A' = S^{-1}
for (tp in c("sym", "chol", "randrot")) {
  A <- root(S$Sc, tp)
  res[paste0("root_", tp)] <- ok(paste("root", tp, ": A A' = S^-1"),
                                 max(abs(A %*% t(A) - solve(S$Sc))) < 1e-10)
}

## 4. Identity sum_j d_ij = sum_k e_ik = D_i^2 for any roots
Sr0 <- crossprod(matrix(rnorm(p * p), p)) + diag(p); Sc0 <- crossprod(matrix(rnorm(q * q), q)) + diag(q)
for (tp in c("sym", "chol", "randrot")) {
  D <- distances(X, Sr0, Sc0, tp, rot = list(c = rorth(q), r = rorth(p)))
  D2 <- sapply(1:30, function(i) sum(diag(solve(Sc0, t(X[, , i]) %*% solve(Sr0, X[, , i])))))
  e1 <- max(abs(rowSums(D$d) - rowSums(D$e))); e2 <- max(abs(rowSums(D$d) - D2))
  res[paste0("ident_", tp)] <- ok(sprintf("identity sum d = sum e = D^2 (%s): %.1e, %.1e", tp, e1, e2),
                                  e1 < 1e-9 && e2 < 1e-9)
}

## 5. No trimming: tff == flipflop (M = 0), every rule and root, both denominators
ff0 <- flipflop(X, scale_id = "det")
for (rl in c("none", "R0", "R1", "R2", "R3", "RJ")) for (tp in c("sym", "chol", "randrot")) {
  f <- tff(X, rule = rl, beta = 1, root_type = tp, init = "identity", maxit = 2000, tol = 1e-12)
  e <- max(abs(f$Sc %x% f$Sr - ff0$Sc %x% ff0$Sr))
  res[paste("noTrim", rl, tp)] <- ok(sprintf("beta = 1: tff(%s, %s) == flipflop (max diff %.1e, %s)",
                                             rl, tp, e, f$status), e < 1e-7)
}

## 6. kappa: kappa(1, m) = 1, and Monte Carlo check of the trimmed second moment
res["kappa1"] <- ok("kappa(1, m) = 1", kappa_mcd(1, 5) == 1)
x <- rchisq(1e6, 5); k <- kappa_mcd(0.75, 5)
mc <- k * mean(x[x <= quantile(x, 0.75)]) / 5
res["kappa_mc"] <- ok(sprintf("kappa(0.75, 5) * E[d | d <= q75] / p = %.4f (should be 1)", mc), abs(mc - 1) < 0.01)

## 7. RJ objective is monotone (a quick sanity run)
dat <- gen_data(25, S$Sr, S$Sc, "col", 0.1)
f <- tff(dat$X, rule = "RJ", beta = 0.75, init = "identity")
res["rj_mono"] <- ok(sprintf("RJ: J* non-increasing (max increase %.2e)", max(diff(f$J))),
                     max(diff(f$J)) <= 1e-8 * abs(f$J[1]))

## 8. The two RJ implementations minimise the same J*: tff(rule = "RJ") vs fit_rj() in exp/RJ_simple.R
# They are not step-for-step identical: tff minimises the A block exactly (rj_update_A, up to 20 coordinate
# sweeps), fit_rj does one sweep per iteration, so they may stop at different local minima of J*.
# Test: fit_rj's end point is a fixed point of tff's exact block updates (C-step, Sr step, A step).
rjs <- new.env()
for (e in parse("exp/RJ_simple.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function")))
    eval(e, rjs)
}
dat <- gen_data(40, S$Sr, S$Sc, "col", 0.1)
f1 <- rjs$fit_rj(dat$X, beta = 0.75, tol = 1e-14, maxit = 5000)
Xc <- sweep(dat$X, 1:2, f1$M)
n1 <- dim(Xc)[3]; h1 <- sum(f1$keep[, 1]); q1 <- dim(Xc)[2]
d1 <- rjs$rj_dist(Xc, f1$A, f1$Sr)
keep_c <- apply(d1, 2, function(dj) rank(dj, ties.method = "first") <= h1)
Z1 <- zmat(Xc, f1$A)
Sr_c <- tcrossprod(Z1[, as.vector(f1$keep), drop = FALSE]) / (h1 * q1)
V1 <- array(0, c(q1, q1, q1))
for (i in seq_len(n1)) for (j in which(f1$keep[i, ])) V1[, , j] <- V1[, , j] + crossprod(Xc[, , i], solve(f1$Sr, Xc[, , i]))
A_c <- rj_update_A(f1$A, V1, h1 * q1, dim(Xc)[1], q1, sweeps = 200)
e <- max(abs(Sr_c - f1$Sr)) + max(abs(A_c - f1$A)) / max(abs(f1$A))
f2 <- tff(Xc, rule = "RJ", beta = 0.75, init = "identity", maxit = 500, tol = 1e-12, trace = FALSE)
J1 <- rjs$rj_crit(Xc, f1$A, f1$Sr, f1$keep)
J2 <- rjs$rj_crit(Xc, f2$A, f2$Sr / kappa_mcd(0.75, dim(Xc)[1]), f2$a)
res["rj_two_impl"] <- ok(sprintf(paste("RJ: fit_rj end point is a fixed point of tff's block updates",
                                       "(same C-step: %s, step size %.1e); J* fit_rj %.3f vs tff %.3f"),
                                 identical(keep_c, f1$keep), e, J1, J2),
                         identical(keep_c, f1$keep) && e < 1e-5)

cat(sprintf("\n%d / %d tests passed\n", sum(res), length(res)))
if (!all(res)) quit(status = 1)
