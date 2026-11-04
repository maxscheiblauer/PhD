# Efficiency diagnostic (block 2): clean 3x4, n = 1000, 8 seeds. Which ingredient costs the efficiency
# relative to the MLE? Variants: default (quant .99, marg_cut 3, obs trim .99); no forced marginal flags
# (marg_cut = Inf); quant = .999; no observation trimming; all three relaxed. Output: output/cm_step3c_eff.log
source("R/load.R"); library(parallel); set.seed(34)
sink("output/cm_step3c_eff.log")
p <- 3; q <- 4; S <- true_sigma(p, q, "full"); n <- 1000
vars <- list(default = list(), no_marg = list(marg_cut = Inf), q999 = list(quant = 0.999), no_obs = list(alpha_obs = 1),
             relaxed = list(marg_cut = Inf, quant = 0.999, alpha_obs = 1))
r <- mclapply(1:8, function(k) { set.seed(3400 + k); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); ff <- flipflop(X)
  out <- c(MLE = kl_mn(ff$Sr, ff$Sc, S$Sr, S$Sc))
  for (v in names(vars)) { f <- do.call(cellmmcd, c(list(X = X), vars[[v]])); out[v] <- kl_mn(f$Sr, f$Sc, S$Sr, S$Sc); out[paste0(v, "_fl")] <- f$flagged; out[paste0(v, "_tr")] <- tr_ratio(f$Sr, f$Sc, S$Sr, S$Sc) }
  out }, mc.cores = 4)
v <- colMeans(do.call(rbind, r))
cat(sprintf("3x4 n=1000 mean KL: MLE %.4f\n", v["MLE"]))
for (nm in names(vars)) cat(sprintf("  %-8s KL %.4f (ratio to MLE %.2f), flagged %.4f, tr ratio %.3f\n", nm, v[nm], v[nm] / v["MLE"], v[paste0(nm, "_fl")], v[paste0(nm, "_tr")]))
sink()
