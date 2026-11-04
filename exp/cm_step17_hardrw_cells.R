# Block 6: reweighted hard quota vs the soft default on the cellwise grid (the interplay / lock-in cases of
# blocks 3-4 are the risk: the raw hard quota always trims 25 % of the matrices; the reweighting must re-admit
# them under 10-20 % cellwise contamination). p = q = 5 and 3x4, n = 100, 10 reps, cells at 6 sd (and 3 / 4 sd
# at 10 %). Output: output/cm_step17_hardrw_cells.log / .rds
source("R/load.R"); library(parallel); set.seed(17)
one_rep <- function(seed, p, q, n, eps, k) {
  set.seed(seed); S <- true_sigma(p, q, "full"); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
  bad <- array(FALSE, c(p, q, n)); m <- round(eps * p * q * n)
  if (m > 0) { idx <- sample(p * q * n, m); bad[idx] <- TRUE
    sdc <- array(sqrt(outer(diag(S$Sr), diag(S$Sc))), c(p, q, n)); X[idx] <- X[idx] + k * sdc[idx] * sample(c(-1, 1), m, TRUE) }
  ev <- function(f) c(kl = kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), tr = tr_ratio(f$Sr, f$Sc, S$Sr, S$Sc), trimmed = f$trimmed,
                      recall = if (m > 0) mean(!f$W[bad]) else NA, fp = mean(!f$W[!bad & array(rep(f$v, each = p * q), c(p, q, n))]))
  list(soft = ev(cellmmcd(X)), hard_rw = ev(cellmmcd(X, quant_obs = NULL)))
}
grid <- rbind(data.frame(p = 5, q = 5, eps = c(0, 0.05, 0.1, 0.2), k = 6), data.frame(p = 5, q = 5, eps = 0.1, k = c(3, 4)),
              data.frame(p = 3, q = 4, eps = c(0, 0.1, 0.2), k = 6))
n <- 100; R <- 10; results <- list()
sink("output/cm_step17_hardrw_cells.log"); cat("cellMMCD step 17 (block 6): reweighted hard quota on the cellwise grid, n = 100, 10 reps. Medians.\n")
for (g in seq_len(nrow(grid))) {
  reps <- mclapply(1:R, function(r) one_rep(17000 + 100 * g + r, grid$p[g], grid$q[g], n, grid$eps[g], grid$k[g]), mc.cores = 4)
  cat(sprintf("\n== %dx%d, cells %.2f at %d sd ==\n", grid$p[g], grid$q[g], grid$eps[g], grid$k[g]))
  tab <- t(sapply(names(reps[[1]]), function(m) apply(do.call(rbind, lapply(reps, `[[`, m)), 2, median, na.rm = TRUE)))
  print(round(tab, 3)); results[[g]] <- list(grid = grid[g, ], reps = reps)
}
saveRDS(results, "output/cm_step17_hardrw_cells.rds"); sink()
