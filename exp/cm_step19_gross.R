# Block 6: gross cells beyond the per-position quota. The marginal cut (3 MAD per position) force-flags every
# gross cell as long as fewer than half of the cells at the position are contaminated, whatever the quota
# (JOURNAL block 6, result 6). Does the estimator survive 30-45 % gross cells per position? Variants: default
# (with the observation step) and cells only (alpha_obs = 1). 3x4, n = 100, 5 reps, cells at 1000 sd placed
# position-wise (exactly round(eps * n) per position). Output: output/cm_step19_gross.log
source("R/load.R"); library(parallel); set.seed(19)
one <- function(seed, eps, p = 3, q = 4, n = 100) {
  set.seed(seed); S <- true_sigma(p, q, "full"); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); sdc <- sqrt(outer(diag(S$Sr), diag(S$Sc)))
  bad <- array(FALSE, c(p, q, n)); m <- round(eps * n)
  for (k in 1:p) for (l in 1:q) { idx <- sample(n, m); bad[k, l, idx] <- TRUE; X[k, l, idx] <- X[k, l, idx] + 1000 * sdc[k, l] * sample(c(-1, 1), m, TRUE) }
  ev <- function(f) c(kl = kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), tr = tr_ratio(f$Sr, f$Sc, S$Sr, S$Sc), trimmed = f$trimmed, recall = mean(!f$W[bad]), fp = mean(!f$W[!bad]))
  list(default = ev(cellmmcd(X)), cells_only = ev(cellmmcd(X, alpha_obs = 1)), hard_rw = ev(cellmmcd(X, quant_obs = NULL)))
}
sink("output/cm_step19_gross.log"); cat("Gross cells (1000 sd) at exactly eps*n per position, 3x4, n = 100, 5 reps. Medians.\n")
for (eps in c(0.2, 0.3, 0.4, 0.45)) { reps <- mclapply(1:5, function(r) one(19000 + 100 * eps * 100 + r, eps), mc.cores = 4)
  cat(sprintf("\n== eps = %.2f ==\n", eps)); print(round(t(sapply(names(reps[[1]]), function(m) apply(do.call(rbind, lapply(reps, `[[`, m)), 2, median))), 3)) }
sink()
