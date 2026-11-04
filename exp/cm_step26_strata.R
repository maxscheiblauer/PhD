# Block 9, step 26: the stratum criterion for the existence problem (JOURNAL block 9, results 1-2). (a) Prop. D (dim N_i on a
# stratum) against the nullity of G_0 on random patterns; (b) the slope of F along the (a,b)-path (Prop. C') on the 2x2 and the
# step-25 patterns; (c) Conjecture C'' on random 3x3 patterns: strata that the count flags, (a,b) found by least squares on the
# orthogonality system, slope of F along the path vs the predicted nu; (d) algorithm-produced patterns (3x3 n=6, 4x4 n=12,
# 5x5 n=15 and n=100 with 20 % cells): worst stratum and slack; (e) the random-contamination regime at 20x20, n = 100. Seed 26.
# Output: output/cm_step26_strata.log
source("R/load.R"); set.seed(26)
Fobs <- function(X, W, Sr, Sc) { Sig <- Sc %x% Sr; tot <- 0; for (i in seq_len(dim(X)[3])) { o <- which(W[, , i]); x <- as.vector(X[, , i])[o]
  So <- Sig[o, o, drop = FALSE]; tot <- tot + logdet(So) + sum(x * solve(So, x)) }; tot }
Fpath <- function(X, W, a, b, d) Fobs(X, W, diag(length(a)) - (1 - d) * tcrossprod(a), if (is.null(b)) diag(dim(X)[2]) else diag(length(b)) - (1 - d) * tcrossprod(b))
slope <- function(X, W, a, b) (Fpath(X, W, a, b, 1e-6) - Fpath(X, W, a, b, 1e-3)) / log(1e-3)
nu_of <- function(W, T, U) sum(sapply(seq_len(dim(W)[3]), function(i) cm_null_dim(W[, , i], T, U) + as.numeric(length(U) > 0 && all(W[T, U, i]))))
find_ab <- function(X, W, T, U, nstart = 8) { p <- dim(W)[1]; q <- dim(W)[2]; n <- dim(W)[3]
  unpack <- function(th) { a <- numeric(p); b <- numeric(q); a[T] <- th[seq_along(T)]; b[U] <- th[length(T) + seq_along(U)]
    a <- a / sqrt(sum(a^2)); if (length(U)) b <- b / sqrt(sum(b^2)); list(a = a, b = b) }
  res <- function(th) { ab <- unpack(th); tot <- 0; for (i in 1:n) { B <- cm_null_basis(W[, , i], T, U, ab$a, ab$b)
    if (!is.null(B)) tot <- tot + sum(crossprod(B, as.vector(X[, , i]))^2) }; tot }
  best <- NULL; for (s in 1:nstart) { o <- optim(rnorm(length(T) + length(U)), res, method = "BFGS", control = list(maxit = 500, reltol = 1e-14))
    if (is.null(best) || o$value < best$value) best <- o }
  ab <- unpack(best$par); list(a = ab$a, b = if (length(U)) ab$b else NULL, resid = best$value, full = min(abs(c(ab$a[T], ab$b[U]))) > 1e-3) }
sink("output/cm_step26_strata.log"); cat("cellMMCD step 26 (block 9): stratum criterion, Prop. C' slopes, Conjecture C'' on random patterns\n")
cat("\n(a) Prop. D vs nullity of G_0 on 40 random patterns (p,q in 2..4, kept prob .6, random strata incl. one-factor):\n")
agree <- 0; for (r in 1:40) { p <- sample(2:4, 1); q <- sample(2:4, 1); Wi <- matrix(runif(p * q) < 0.6, p, q); Wi[1, 1] <- TRUE
  T <- if (r %% 7 == 0) integer() else sort(sample(p, sample(p, 1))); U <- if (r %% 5 == 0 && length(T)) integer() else sort(sample(q, sample(q, 1)))
  a <- numeric(p); if (length(T)) { a[T] <- rnorm(length(T)); a <- a / sqrt(sum(a^2)) }; b <- numeric(q); if (length(U)) { b[U] <- rnorm(length(U)); b <- b / sqrt(sum(b^2)) }
  P0 <- (if (length(U)) diag(q) - tcrossprod(b) else diag(q)) %x% (if (length(T)) diag(p) - tcrossprod(a) else diag(p)); o <- which(Wi)
  nul <- sum(eigen(P0[o, o, drop = FALSE], symmetric = TRUE, only.values = TRUE)$values < 1e-9); d <- cm_null_dim(Wi, T, U); agree <- agree + (nul == d) }
cat("  agreement:", agree, "/ 40\n")
cat("\n(b) slopes of F along the (a,b)-path (predicted nu = sum_i dim N_i + #delta^2 modes):\n")
W24 <- array(c(FALSE, TRUE, TRUE, FALSE), c(2, 2, 1)); X24 <- array(rnorm(4), c(2, 2, 1)); f24 <- find_ab(X24, W24, 1:2, 1:2)
cat(sprintf("  2x2 cross (step 24): resid %.1e, full support %s, slope %.3f, nu %d\n", f24$resid, f24$full, slope(X24, W24, f24$a, f24$b), nu_of(W24, 1:2, 1:2)))
W25 <- array(FALSE, c(3, 3, 3)); W25[2:3, 1, 1] <- TRUE; W25[1, 2:3, 1] <- TRUE; W25[2:3, 1, 2] <- TRUE; W25[1, 2:3, 3] <- TRUE; X25 <- array(rnorm(27), c(3, 3, 3))
f25 <- find_ab(X25, W25, 1:3, 1:3); cat(sprintf("  3x3 step-25 pattern: resid %.1e, full %s, slope %.3f, nu %d\n", f25$resid, f25$full, slope(X25, W25, f25$a, f25$b), nu_of(W25, 1:3, 1:3)))
a <- c(0.6, 0.8); b <- c(-0.8, 0.6); Xd <- array(tcrossprod(c(-0.8, 0.6), c(0.6, 0.8)), c(2, 2, 1)); Wd <- array(TRUE, c(2, 2, 1))
cat(sprintf("  2x2 all kept, y = a_perp b_perp' (delta^2 mode): slope %.3f, nu %d\n", (Fpath(Xd, Wd, a, b, 1e-5) - Fpath(Xd, Wd, a, b, 1e-3)) / log(1e-2), nu_of(Wd, 1:2, 1:2)))
cat("\n(c) Conjecture C'' on random 3x3 patterns (kept prob .55; 40 patterns per n):\n")
for (n in 2:4) { cnt <- c(A = 0, cross = 0, silent = 0); shown <- 0
  for (r in 1:40) { W <- array(runif(9 * n) < 0.55, c(3, 3, n)); X <- array(rnorm(9 * n), c(3, 3, n)); cs <- cm_check_strata(W)
    cls <- if (cs$viol_A == 1) "A" else if (cs$viol_cross == 1) "cross" else "silent"; cnt[cls] <- cnt[cls] + 1
    if (cls == "cross" && shown < 6) { shown <- shown + 1; w <- cs$worst; f <- find_ab(X, W, w$T, w$U); nu <- nu_of(W, w$T, w$U)
      ctrl <- { a2 <- rnorm(3); a2 <- a2 / sqrt(sum(a2^2)); b2 <- rnorm(3); b2 <- b2 / sqrt(sum(b2^2)); slope(X, W, a2, b2) }
      cat(sprintf("  n=%d pattern %2d: worst stratum T={%s} U={%s} count %d need %d; (a,b) resid %.1e full %s; slope %.3f (nu %d; control slope %.0f)\n",
        n, r, paste(w$T, collapse = ""), paste(w$U, collapse = ""), w$count, w$need, f$resid, f$full, slope(X, W, f$a, f$b), nu, ctrl)) } }
  cat(sprintf("  n=%d: Prop. A fires %d, cross-only %d, silent %d\n", n, cnt["A"], cnt["cross"], cnt["silent"])) }
cat("\n(d) algorithm-produced patterns (default cellmmcd), worst stratum over all (T,U):\n")
run_d <- function(p, q, n, eps, seed) { set.seed(seed); S <- true_sigma(p, q); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
  if (eps > 0) { idx <- sample(p * q * n, round(eps * p * q * n)); ai <- arrayInd(idx, dim(X)); X[idx] <- X[idx] + 6 * sample(c(-1, 1), length(idx), TRUE) * sqrt(diag(S$Sr)[ai[, 1]] * diag(S$Sc)[ai[, 2]]) }
  fit <- suppressWarnings(cellmmcd(X)); t0 <- proc.time()[3]; cs <- cm_check_strata(fit); w <- cs$worst
  cat(sprintf("  %dx%d n=%3d eps=%.2f seed %d: kept matrices %d, c_full %d, d_full %d, viol_A %d, viol_cross %d, min slack %s (T={%s} U={%s} count %d need %d), %.1f s\n",
    p, q, n, eps, seed, sum(fit$v), cs$A$c_full, cs$A$d_full, cs$viol_A, cs$viol_cross, format(cs$slack), paste(w$T, collapse = ""), paste(w$U, collapse = ""), w$count, w$need, proc.time()[3] - t0)) }
for (s in 1:3) run_d(3, 3, 6, 0, 2400 + s); for (s in 1:2) run_d(4, 4, 12, 0, 2600 + s); for (s in 1:2) run_d(5, 5, 15, 0, 2610 + s); for (s in 1:2) run_d(5, 5, 100, 0.2, 2620 + s)
cat("\n(e) random-contamination regime: expected fully kept columns n q (1-eps)^p vs p (n = 100), and fitted 20x20 patterns:\n")
for (p in c(5, 8, 20, 40)) cat(sprintf("  p=q=%2d: eps .05 %7.1f | .10 %7.1f | .20 %7.1f | .25 %7.1f   (needs 0 or >= %d)\n", p, 100 * p * .95^p, 100 * p * .9^p, 100 * p * .8^p, 100 * p * .75^p, p))
for (eps in c(0.2, 0.25)) { set.seed(2630); p <- q <- 20; n <- 100; S <- true_sigma(p, q); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
  idx <- sample(p * q * n, round(eps * p * q * n)); ai <- arrayInd(idx, dim(X)); X[idx] <- X[idx] + 6 * sample(c(-1, 1), length(idx), TRUE) * sqrt(diag(S$Sr)[ai[, 1]] * diag(S$Sc)[ai[, 2]])
  t0 <- proc.time()[3]; fit <- suppressWarnings(cellmmcd(X)); A <- cm_check_pattern(fit)
  cat(sprintf("  20x20 n=100 eps=%.2f: c_full %d, d_full %d (expected before false flags %.1f); Prop. A fires on T=[p] or U=[q]: %s; KL %.3f; kept matrices %d; %.1f s\n",
    eps, A$c_full, A$d_full, n * q * (1 - eps)^p, (A$c_full >= 1 && A$c_full <= p - 1) || (A$d_full >= 1 && A$d_full <= q - 1), kl_mn(fit$Sr, fit$Sc, S$Sr, S$Sc), sum(fit$v), proc.time()[3] - t0)) }
sink(); cat("done\n")
