# Block 3: is the efficiency gap at the Gaussian bias or variance? Clean 3x4 ("full" shape), n = 500, 2000,
# 8000, 4 seeds. If n * KL grows with n, the fixed point is biased (the trimming is not at random: a flagged
# clean cell is in the tail of its conditional, so treating it as missing shrinks the conditional variances,
# Sigma - (F1 - F3) diag(1 / Omega_cc) to first order, which is not a uniform scale). Variants: default
# (quant .99), q999, tail (censored parameter step, no t1), rw (reweighted at .999). Output: cm_step6_bias.log
source("R/load.R"); library(parallel); set.seed(36)
sink("output/cm_step6_bias.log")
p <- 3; q <- 4; S <- true_sigma(p, q, "full")
vars <- list(default = list(), q999 = list(quant = 0.999), tail = list(tail = TRUE), rw = list(reweight = TRUE))
grid <- expand.grid(seed = 1:4, n = c(500, 2000, 8000))
r <- mclapply(seq_len(nrow(grid)), function(g) { n <- grid$n[g]; set.seed(3600 + grid$seed[g])
  X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); ff <- flipflop(X)
  out <- c(n = n, MLE = kl_mn(ff$Sr, ff$Sc, S$Sr, S$Sc))
  for (v in names(vars)) { t0 <- proc.time()[3]; f <- do.call(cellmmcd, c(list(X = X), vars[[v]]))
    out[v] <- kl_mn(f$Sr, f$Sc, S$Sr, S$Sc); out[paste0(v, "_sh")] <- kl_shape(f$Sr, f$Sc, S$Sr, S$Sc)
    out[paste0(v, "_fl")] <- f$flagged; out[paste0(v, "_tr")] <- tr_ratio(f$Sr, f$Sc, S$Sr, S$Sc); out[paste0(v, "_t")] <- proc.time()[3] - t0 }
  out }, mc.cores = 4)
R <- do.call(rbind, r); saveRDS(R, "output/cm_step6_bias.rds")
cat("Clean 3x4 full shape, 4 seeds. Mean KL, n*KL, shape KL, flagged fraction, tr ratio, seconds per fit.\n")
for (n in unique(grid$n)) { v <- colMeans(R[R[, "n"] == n, , drop = FALSE])
  cat(sprintf("n = %5d  MLE: KL %.5f  n*KL %.2f\n", n, v["MLE"], n * v["MLE"]))
  for (nm in names(vars)) cat(sprintf("   %-7s KL %.5f (ratio %.2f) n*KL %6.2f  shape %.5f  flagged %.4f  tr %.3f  %.0fs\n", nm, v[nm], v[nm] / v["MLE"], n * v[nm], v[paste0(nm, "_sh")], v[paste0(nm, "_fl")], v[paste0(nm, "_tr")], v[paste0(nm, "_t")])) }
sink()
