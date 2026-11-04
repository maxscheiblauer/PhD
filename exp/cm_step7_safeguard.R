# Block 3: does the full-column safeguard (cm_fullcol) prevent the degenerate descents of the existence grid
# (exp/cm_step5_exist.R), and what does it cost? Clean data, alpha = .75, 4 seeds; safeguard off vs on.
# degen as in step 5; unfl = cells unflagged by the safeguard; nonmono = fraction of runs where F increased
# at some iteration by more than 1e-6 |F|; KL vs the MLE at n = 20 and 40. Output: output/cm_step7_safeguard.log
source("R/load.R"); library(parallel); set.seed(37)
one <- function(p, q, n, seed, sg) {
  set.seed(seed); S <- true_sigma(p, q, "full"); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
  f <- tryCatch(suppressWarnings(cellmmcd(X, maxit = 300, safeguard = sg)), error = function(e) NULL)
  if (is.null(f)) return(c(err = 1, degen = NA, unfl = NA, nonmono = NA, kl = NA, klmle = NA))
  ev <- function(A) { e <- eigen(A, symmetric = TRUE, only.values = TRUE)$values; min(e) / max(e) }
  ff <- tryCatch(flipflop(X), error = function(e) NULL)
  c(err = 0, degen = as.numeric(min(ev(f$Sr), ev(f$Sc)) < 1e-6 || f$obj[length(f$obj)] < f$obj[1] - 50 * n),
    unfl = f$unflagged, nonmono = as.numeric(any(diff(f$obj) > 1e-6 * abs(f$obj[-1]))),
    kl = kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), klmle = if (is.null(ff)) NA else kl_mn(ff$Sr, ff$Sc, S$Sr, S$Sc))
}
grid <- expand.grid(pq = c("3x3", "5x5", "8x8", "3x8"), n = c(6, 9, 12, 20, 40), stringsAsFactors = FALSE)
sink("output/cm_step7_safeguard.log"); cat("cellMMCD step 7: full-column safeguard on the existence grid, alpha = .75, 4 seeds per cell.\n")
cat(sprintf("%-5s %3s | %5s %7s | %5s %7s %7s %7s %7s\n", "pq", "n", "degen", "KL_off", "degen", "unfl", "nonmono", "KL_on", "KL_MLE"))
for (g in seq_len(nrow(grid))) {
  p <- as.integer(substr(grid$pq[g], 1, 1)); q <- as.integer(substr(grid$pq[g], 3, 3)); n <- grid$n[g]
  r0 <- do.call(rbind, mclapply(1:4, function(s) one(p, q, n, 7000 + 10 * g + s, FALSE), mc.cores = 4))
  r1 <- do.call(rbind, mclapply(1:4, function(s) one(p, q, n, 7000 + 10 * g + s, TRUE), mc.cores = 4))
  m <- function(r, c) mean(r[, c], na.rm = TRUE); md <- function(r, c) median(r[, c], na.rm = TRUE)
  cat(sprintf("%-5s %3d | %5.2f %7.2f | %5.2f %7.1f %7.2f %7.2f %7.2f\n", grid$pq[g], n, m(r0, "degen"), md(r0, "kl"), m(r1, "degen"), m(r1, "unfl"), m(r1, "nonmono"), md(r1, "kl"), md(r1, "klmle")))
}
sink()
