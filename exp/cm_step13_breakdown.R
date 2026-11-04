# cellMMCD step 13 (block 5): numerical breakdown under cellwise and casewise contamination. Seed 13.
# p = q = 5, n = 100, 5 reps. Cellwise: fraction eps of cells (uniform positions) shifted by 4 / 10 / 1000 sd with random sign,
# eps from 0.10 to 0.50. Casewise: fraction of matrices shifted by 5 sd in every cell, 0.10 to 0.45. Default estimator
# (alpha = alpha_obs = 0.75: per-position cell quota 25 %, observation quota 25 %; forced marginal flags at |x|/MAD > 3).
# Reported: median and max KL, median flagged fraction, recall of bad cells, trimmed matrices; the start kpair and the MLE for reference.
# Output: output/cm_step13_breakdown.log, output/cm_step13_breakdown.rds
source("R/load.R"); library(parallel); set.seed(13)
p <- q <- 5; n <- 100; R <- 5; S <- true_sigma(p, q, "full"); sd_cell <- sqrt(outer(diag(S$Sr), diag(S$Sc)))
one <- function(seed, mode, eps, mag) { set.seed(seed); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); bad <- array(FALSE, c(p, q, n))
  if (mode == "cell") { bad <- array(runif(p * q * n) < eps, c(p, q, n)); X[bad] <- X[bad] + mag * sd_cell[array(1:(p * q), c(p, q, n))[bad]] * sign(rnorm(sum(bad))) }
  else { ib <- seq_len(round(eps * n)); X[, , ib] <- X[, , ib] + mag * sd_cell[rep(TRUE, p * q)]; bad[, , ib] <- TRUE }
  f <- tryCatch(cellmmcd(X), error = function(e) NULL); st <- kpair(X); ff <- flipflop(X)
  if (is.null(f)) return(c(kl = NA, flagged = NA, rec = NA, tr = NA, kl_start = kl_mn(st$Sr, st$Sc, S$Sr, S$Sc), kl_mle = kl_mn(ff$Sr, ff$Sc, S$Sr, S$Sc), err = 1))
  fl <- !f$W | array(rep(!f$v, each = p * q), c(p, q, n))
  c(kl = kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), flagged = mean(!f$W[, , f$v]), rec = mean(fl[bad]), tr = mean(!f$v),
    kl_start = kl_mn(st$Sr, st$Sc, S$Sr, S$Sc), kl_mle = kl_mn(ff$Sr, ff$Sc, S$Sr, S$Sc), err = 0) }
grid <- rbind(expand.grid(mode = "cell", eps = c(0.10, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50), mag = c(4, 10, 1000)),
              expand.grid(mode = "case", eps = c(0.10, 0.20, 0.25, 0.30, 0.35, 0.45), mag = 5))
res <- list(); sink("output/cm_step13_breakdown.log")
cat("cellMMCD step 13: p = q = 5, n = 100, 5 reps. KL median [max]; flagged = fraction of cells flagged among kept matrices; rec = recall of bad cells (flagged or trimmed).\n")
for (g in seq_len(nrow(grid))) {
  r <- do.call(rbind, mclapply(1:R, function(s) one(13000 + 10 * g + s, as.character(grid$mode[g]), grid$eps[g], grid$mag[g]), mc.cores = 4)); res[[g]] <- list(grid = grid[g, ], r = r)
  cat(sprintf("%s eps %.2f mag %4g: cellMMCD KL %7.3f [%7.3f] flagged %.3f rec %.2f trimmed %.2f errors %d | kpair KL %7.3f | MLE KL %8.2f\n", grid$mode[g], grid$eps[g], grid$mag[g],
              median(r[, "kl"], na.rm = TRUE), max(r[, "kl"], na.rm = TRUE), median(r[, "flagged"], na.rm = TRUE), median(r[, "rec"], na.rm = TRUE), median(r[, "tr"], na.rm = TRUE), sum(r[, "err"]),
              median(r[, "kl_start"]), median(r[, "kl_mle"]))) }
sink(); saveRDS(res, "output/cm_step13_breakdown.rds"); cat("done\n")
