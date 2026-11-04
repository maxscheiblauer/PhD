# Block 6: consistency of the block-5 default (robust rate, count rule, shadow step, joint location,
# reweighting) at large n. Clean 3x4 full shape, n = 500, 2000, 8000, 4 seeds (as cm_step6_bias, block 3).
# Variants: default, loc = FALSE, hard quota (quant_obs = NULL; block 6: with observation reweighting),
# hard quota without observation reweighting. Criterion: n*KL bounded; ratio to the MLE of order 1-1.5.
# Output: output/cm_step15_consistency.log / .rds
source("R/load.R"); library(parallel); set.seed(36)
sink("output/cm_step15_consistency.log")
p <- 3; q <- 4; S <- true_sigma(p, q, "full")
vars <- list(default = list(), nloc = list(loc = FALSE), hard_rw = list(quant_obs = NULL), hard_raw = list(quant_obs = NULL, quant_rw_obs = FALSE))
grid <- expand.grid(seed = 1:4, n = c(500, 2000, 8000))
r <- mclapply(seq_len(nrow(grid)), function(g) { n <- grid$n[g]; set.seed(3600 + grid$seed[g])
  X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); ff <- flipflop(X)
  out <- c(n = n, MLE = kl_mn(ff$Sr, ff$Sc, S$Sr, S$Sc))
  for (v in names(vars)) { t0 <- proc.time()[3]; f <- do.call(cellmmcd, c(list(X = X), vars[[v]]))
    out[v] <- kl_mn(f$Sr, f$Sc, S$Sr, S$Sc); out[paste0(v, "_sh")] <- kl_shape(f$Sr, f$Sc, S$Sr, S$Sc)
    out[paste0(v, "_fl")] <- f$flagged; out[paste0(v, "_trim")] <- f$trimmed; out[paste0(v, "_tr")] <- tr_ratio(f$Sr, f$Sc, S$Sr, S$Sc)
    out[paste0(v, "_M")] <- sum(f$M^2); out[paste0(v, "_t")] <- proc.time()[3] - t0 }
  out }, mc.cores = 4)
R <- do.call(rbind, r); saveRDS(R, "output/cm_step15_consistency.rds")
cat("Clean 3x4 full shape, 4 seeds, block-6 code. Mean KL, n*KL, shape KL, flagged, trimmed, tr ratio, n*|M|^2, seconds per fit.\n")
for (n in unique(grid$n)) { v <- colMeans(R[R[, "n"] == n, , drop = FALSE])
  cat(sprintf("n = %5d  MLE: KL %.5f  n*KL %.2f\n", n, v["MLE"], n * v["MLE"]))
  for (nm in names(vars)) cat(sprintf("   %-8s KL %.5f (ratio %.2f) n*KL %6.2f  shape %.5f  flagged %.4f  trimmed %.4f  tr %.3f  n|M|^2 %.2f  %.0fs\n", nm, v[nm], v[nm] / v["MLE"], n * v[nm], v[paste0(nm, "_sh")], v[paste0(nm, "_fl")], v[paste0(nm, "_trim")], v[paste0(nm, "_tr")], n * v[paste0(nm, "_M")], v[paste0(nm, "_t")])) }
sink()
