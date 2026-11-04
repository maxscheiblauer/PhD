# Block 8, step 24: (a) the counterexample to Conjecture C (JOURNAL block 8, result 1): p = q = 2, n = 1, kept cells
# (2,1) and (1,2) (no partial column or row of any size, so Prop. A is silent), F along the explicit path
# Sigma_OO = yy' + tI realised in the Kronecker cone, and by direct minimisation (BFGS over log-Cholesky factors);
# (b) the same pattern with n = 2: bounded, the minimum equals the unstructured 2-D MLE value 2 logdet S + 4;
# (c) location-free counting (STATE-7 item 4): does the joint descent with M free reach an unbounded configuration at
# 3x3, n = 6 (Prop. A remark iii: sum_l (N_l - 1)_+ <= p - 1 makes the joint objective unbounded)? Seed 24.
# Output: output/cm_step24_crossmode.log
source("R/load.R"); set.seed(24)
sink("output/cm_step24_crossmode.log"); cat("cellMMCD step 24 (block 8): cross-mode counterexample and location-free counting\n")
chol_par <- function(th, m) { L <- matrix(0, m, m); L[lower.tri(L, TRUE)] <- th; diag(L) <- exp(diag(L)); tcrossprod(L) }
fobj <- function(th, Xc, W, pp, qq) { k <- pp * (pp + 1) / 2; Sr <- chol_par(th[1:k], pp); Sc <- chol_par(th[-(1:k)], qq)
  tryCatch(cm_obj(Xc, Sr, Sc, W, matrix(0, pp, qq)), error = function(e) 1e10) }
minimise <- function(Xc, W, p, q, maxit = 5000) { th <- rep(0, p * (p + 1) / 2 + q * (q + 1) / 2); path <- c()
  for (r in 1:10) { o <- optim(th, fobj, Xc = Xc, W = W, pp = p, qq = q, method = "BFGS", control = list(maxit = maxit, reltol = 1e-14)); th <- o$par; path <- c(path, o$value) }
  k <- p * (p + 1) / 2; Sr <- chol_par(th[1:k], p); Sc <- chol_par(th[-(1:k)], q)
  list(path = path, evr_r = min(eigen(Sr)$values) / max(eigen(Sr)$values), evr_c = min(eigen(Sc)$values) / max(eigen(Sc)$values)) }
## (a) n = 1
y <- c(1.3, -0.7); W1 <- array(c(FALSE, TRUE, TRUE, FALSE), c(2, 2, 1)); X1 <- array(0, c(2, 2, 1)); X1[2, 1, 1] <- y[1]; X1[1, 2, 1] <- y[2]
cat("\n(a) n = 1, kept cells (2,1), (1,2), y =", y, "; cm_check_pattern:", unlist(cm_check_pattern(W1)[1:5]), "\n")
for (t in 10^-(1:8)) { a <- y[1]^2 + t; d <- y[2]^2 + t; cc <- prod(y); rr <- (1 + abs(cc) / sqrt(a * d)) / 2
  Sr <- matrix(c(1, rr, rr, 1), 2); Sc <- matrix(c(a, cc / rr, cc / rr, d), 2)
  cat(sprintf("  t = %.0e: F = %9.3f  (log t = %7.3f; rho_r = %.6f, rho_c = %.6f)\n", t, cm_obj(X1, Sr, Sc, W1, matrix(0, 2, 2)), log(t), rr, (cc / rr) / sqrt(a * d))) }
m1 <- minimise(X1, W1, 2, 2); cat("  BFGS restarts, F:", round(m1$path, 2), "; eigenvalue ratios", signif(m1$evr_r, 3), signif(m1$evr_c, 3), "\n")
## (b) n = 2
X2 <- array(0, c(2, 2, 2)); Y2 <- matrix(rnorm(4), 2); X2[2, 1, ] <- Y2[1, ]; X2[1, 2, ] <- Y2[2, ]; W2 <- array(c(FALSE, TRUE, TRUE, FALSE), c(2, 2, 2))
S2 <- tcrossprod(Y2) / 2; m2 <- minimise(X2, W2, 2, 2); rS <- S2[1, 2] / sqrt(S2[1, 1] * S2[2, 2]); rr <- sqrt(abs(rS))
FS <- cm_obj(X2, matrix(c(1, rr, rr, 1), 2), matrix(c(S2[1, 1], sign(rS) * rr * sqrt(S2[1, 1] * S2[2, 2]), sign(rS) * rr * sqrt(S2[1, 1] * S2[2, 2]), S2[2, 2]), 2), W2, matrix(0, 2, 2))
cat(sprintf("\n(b) n = 2: S has correlation %.4f, logdet %.3f. F at the Kronecker realisation of Sigma_OO = S: %.4f = unstructured minimum 2 logdet S + 4 + 4 log 2 pi = %.4f (bounded: the cone covers every PD 2x2 Sigma_OO on this pattern).\n    BFGS from the identity (log-Cholesky, 10 restarts) stalls at %s (eigenvalue ratios %.3g %.3g): a generic optimiser does not find the coupled boundary direction; same for (a), where it stalls at %.2f.\n",
  rS, logdet(S2), FS, 2 * logdet(S2) + 4 + 4 * log(2 * pi), paste(unique(round(m2$path, 4)), collapse = " "), m2$evr_r, m2$evr_c, m1$path[10]))
## (c) location-free counting, 3x3, n = 6
cat("\n(c) 3x3, n = 6, clean, default cellmmcd (loc = TRUE), maxit 300 vs 1500: F_end, dF between the two, eigenvalue ratios, N_l, sum (N_l-1)_+, viol\n")
S <- true_sigma(3, 3, "full")
for (s in 1:6) { set.seed(2400 + s); X <- rmatnorm(6, Sr = S$Sr, Sc = S$Sc)
  f1 <- tryCatch(suppressWarnings(cellmmcd(X, maxit = 300)), error = function(e) NULL); f2 <- tryCatch(suppressWarnings(cellmmcd(X, maxit = 1500)), error = function(e) NULL)
  if (is.null(f1) || is.null(f2)) { cat(sprintf("  seed %d: error\n", s)); next }
  Nl <- colSums(apply(f1$W[, , f1$v, drop = FALSE], c(2, 3), all)); Nk <- rowSums(apply(f1$W[, , f1$v, drop = FALSE], c(1, 3), all)); ck <- cm_check_pattern(f1)
  cat(sprintf("  seed %d: F %9.2f, dF(1500-300) %9.2e, evr %.2e %.2e, N_l %s (sum(N_l-1)+ = %d vs p-1 = 2), N_k %s (%d), c_full %d d_full %d viol %d, trimmed %d, flags %d\n", s,
    tail(f1$obj, 1), tail(f2$obj, 1) - tail(f1$obj, 1), min(eigen(f1$Sr)$values) / max(eigen(f1$Sr)$values), min(eigen(f1$Sc)$values) / max(eigen(f1$Sc)$values),
    paste(Nl, collapse = ","), sum(pmax(Nl - 1, 0)), paste(Nk, collapse = ","), sum(pmax(Nk - 1, 0)), ck$c_full, ck$d_full, ck$viol, sum(!f1$v), sum(!f1$W[, , f1$v]))) }
sink(); cat("done\n")
