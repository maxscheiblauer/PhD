# cellMMCD step 1: benchmark under cellwise contamination (shift outliers, 6 sd, random sign).
# Methods: FF (MLE), kpair, cellMMCD, MMCD (robustmatrix), RJ (tff_multi, beta .75),
# cellMCD on vec(X) + nearest Kronecker projection. Seed 1. Output: output/cm_step1.log, output/cm_step1.rds
source("R/load.R"); library(parallel); set.seed(1)
nearest_kron <- function(Sig, p, q) {           # Van Loan rearrangement, rank-1 SVD
  A <- matrix(0, q * q, p * p)
  for (l in 1:q) for (l2 in 1:q) A[l + (l2 - 1) * q, ] <- as.vector(Sig[(l - 1) * p + 1:p, (l2 - 1) * p + 1:p])
  s <- svd(A, nu = 1, nv = 1); Sc <- matrix(s$u[, 1], q, q) * sqrt(s$d[1]); Sr <- matrix(s$v[, 1], p, p) * sqrt(s$d[1])
  if (Sc[1, 1] < 0) { Sc <- -Sc; Sr <- -Sr }
  list(Sr = (Sr + t(Sr)) / 2, Sc = (Sc + t(Sc)) / 2)
}
one_rep <- function(seed, p, q, n, eps) {
  set.seed(seed); S <- true_sigma(p, q, "full")
  X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); truth <- array(FALSE, c(p, q, n))
  m <- round(eps * p * q * n)
  if (m > 0) { bad <- sample(p * q * n, m); truth[bad] <- TRUE
    sdc <- array(sqrt(outer(diag(S$Sr), diag(S$Sc))), c(p, q, n))
    X[bad] <- X[bad] + 6 * sdc[bad] * sample(c(-1, 1), m, TRUE) }
  ev <- function(f) if (is.null(f)) c(kl = NA, skl = NA, tr = NA) else c(kl = kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), skl = kl_shape(f$Sr, f$Sc, S$Sr, S$Sc), tr = tr_ratio(f$Sr, f$Sc, S$Sr, S$Sc))
  tm <- function(expr) { t0 <- proc.time()[3]; r <- tryCatch(expr, error = function(e) { message("fit error: ", conditionMessage(e)); NULL }); list(r = r, t = proc.time()[3] - t0) }
  out <- list()
  a <- tm(flipflop(X)); out$FF <- c(ev(a$r), t = a$t)
  a <- tm(kpair(X)); out$kpair <- c(ev(a$r), t = a$t)
  a <- tm(cellmmcd(X)); out$cellMMCD <- if (is.null(a$r)) c(kl = NA, skl = NA, tr = NA, t = a$t, recall = NA, fpr = NA, iter = NA, flagged = NA) else c(ev(a$r), t = a$t, recall = if (m > 0) mean(!a$r$W[truth]) else NA, fpr = mean(!a$r$W[!truth]), iter = a$r$iter, flagged = a$r$flagged)
  a <- tm({ f <- robustmatrix::mmcd(X, alpha = 0.75); list(Sr = f$cov_row, Sc = f$cov_col) }); out$MMCD <- c(ev(a$r), t = a$t)
  a <- tm({ f <- tff_multi(X, rule = "RJ", beta = 0.75); f$fits[[f$best]] }); out$RJ <- if (is.null(a$r)) c(kl = NA, skl = NA, tr = NA, t = a$t) else c(ev(a$r), t = a$t)
  a <- tm({ f <- cellmcd_quiet(t(matrix(X, p * q)), alpha = 0.75); nearest_kron(f$S, p, q) })
  out$cellMCD_vec <- if (is.null(a$r)) c(kl = NA, skl = NA, tr = NA, t = a$t) else c(ev(a$r), t = a$t)
  out
}
grid <- expand.grid(pq = c("3x4", "5x5"), eps = c(0, 0.05, 0.10, 0.20), stringsAsFactors = FALSE)
n <- 100; R <- 20; results <- list()
sink("output/cm_step1.log"); cellmcd_quiet <- function(...) { tmp <- tempfile(); capture.output(r <- cellWise::cellMCD(...), file = tmp); r }
cat("cellMMCD step 1: n = 100, 20 reps, cellwise shift outliers 6 sd. Medians (IQR in rds).\n")
for (g in seq_len(nrow(grid))) {
  p <- as.integer(substr(grid$pq[g], 1, 1)); q <- as.integer(substr(grid$pq[g], 3, 3)); eps <- grid$eps[g]
  reps <- mclapply(1:R, function(r) one_rep(1000 * g + r, p, q, n, eps), mc.cores = 4)
  meths <- names(reps[[1]])
  tab <- t(sapply(meths, function(m) { v <- do.call(rbind, lapply(reps, `[[`, m)); apply(v, 2, median, na.rm = TRUE)[c("kl", "skl", "tr", "t.elapsed")] }))
  cm <- do.call(rbind, lapply(reps, `[[`, "cellMMCD"))
  cat(sprintf("\n== %s, eps = %.2f ==  cellMMCD: recall %.2f, fpr %.3f, flagged %.3f, iter %.0f (medians)\n", grid$pq[g], eps,
              median(cm[, "recall"]), median(cm[, "fpr"]), median(cm[, "flagged"]), median(cm[, "iter"])))
  print(round(tab, 3))
  results[[g]] <- list(grid = grid[g, ], reps = reps)
}
cat("\n== timing: p = q = 10, n = 200, eps = .10, one run ==\n")
p <- q <- 10; S <- true_sigma(p, q, "full"); set.seed(99); X <- rmatnorm(200, Sr = S$Sr, Sc = S$Sc)
truth <- array(FALSE, c(p, q, 200)); bad <- sample(p * q * 200, 2000); truth[bad] <- TRUE
sdc <- array(sqrt(outer(diag(S$Sr), diag(S$Sc))), c(p, q, 200)); X[bad] <- X[bad] + 6 * sdc[bad] * sample(c(-1, 1), 2000, TRUE)
t0 <- proc.time()[3]; f <- cellmmcd(X); t1 <- proc.time()[3] - t0
cat(sprintf("cellMMCD: %.1f s, %d iterations, KL %.2f (MLE %.2f), recall %.2f, fpr %.3f\n", t1, f$iter,
            kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), kl_mn(flipflop(X)$Sr, flipflop(X)$Sc, S$Sr, S$Sc), mean(!f$W[truth]), mean(!f$W[!truth])))
sink(); saveRDS(results, "output/cm_step1.rds"); cat("done\n")
