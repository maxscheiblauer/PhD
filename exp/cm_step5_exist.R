# cellMMCD step 5: boundedness of the objective under the quota, numerically. Seed 5. Clean data,
# small n; records whether the iteration runs off to a degenerate Sigma (min eigenvalue ratio < 1e-6
# or F dropping below F(MLE start) - 50) within maxit = 300. Also the MLE existence condition
# n >= max(p/q, q/p) + 1 (Dutilleul) for reference. Output: output/cm_step5_exist.log
source("R/load.R"); library(parallel); set.seed(5)
one <- function(p, q, n, alpha, seed) {
  set.seed(seed); S <- true_sigma(p, q, "full"); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
  f <- tryCatch(suppressWarnings(cellmmcd(X, alpha = alpha, maxit = 300)), error = function(e) NULL)
  if (is.null(f)) return(c(err = 1, degen = NA, dF = NA, it = NA, fl = NA))
  ev <- function(A) { e <- eigen(A, symmetric = TRUE, only.values = TRUE)$values; min(e) / max(e) }
  c(err = 0, degen = as.numeric(min(ev(f$Sr), ev(f$Sc)) < 1e-6 || f$obj[length(f$obj)] < f$obj[1] - 50 * n),
    dF = f$obj[1] - f$obj[length(f$obj)], it = f$iter, fl = f$flagged)
}
grid <- expand.grid(pq = c("3x3", "5x5", "8x8", "3x8"), n = c(2, 3, 4, 6, 9, 12), alpha = c(0.5, 0.75), stringsAsFactors = FALSE)
sink("output/cm_step5_exist.log"); cat("cellMMCD step 5: existence grid, 4 seeds per cell. degen = fraction of runs ending degenerate; err = fit errors.\n")
cat(sprintf("%-5s %3s %5s %4s  %5s %5s %7s %6s %5s\n", "pq", "n", "alpha", "h", "degen", "err", "dF_med", "it", "flag"))
for (g in seq_len(nrow(grid))) {
  p <- as.integer(substr(grid$pq[g], 1, 1)); q <- as.integer(substr(grid$pq[g], 3, 3)); n <- grid$n[g]; a <- grid$alpha[g]
  r <- do.call(rbind, mclapply(1:4, function(s) one(p, q, n, a, 5000 + 10 * g + s), mc.cores = 4))
  cat(sprintf("%-5s %3d %5.2f %4d  %5.2f %5.2f %7.1f %6.0f %5.3f   MLE exists: %s\n", grid$pq[g], n, a, ceiling(a * n),
              mean(r[, "degen"], na.rm = TRUE), mean(r[, "err"]), median(r[, "dF"], na.rm = TRUE), median(r[, "it"], na.rm = TRUE), median(r[, "fl"], na.rm = TRUE),
              if (n >= max(p / q, q / p) + 1) "yes" else "no"))
}
sink()
