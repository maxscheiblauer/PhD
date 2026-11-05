# Block 11, step 27c: at eps = .25 (20x20, n = 100, step-27 seeds 2716-2718) the default trims 7-16 % of the matrices
# while the oracle with the true flags (step 27b) is fine: is the count rule (flag_quant) or observation trimming the cause?
# Variants: default; flag_quant = NULL (max_flag 1/2 only); alpha_obs = 1 (no observation trimming at all).
# Output: output/cm_step27c_countrule.log
source("R/load.R"); library(parallel); p <- q <- 20; n <- 100; eps <- 0.25
one <- function(k) { set.seed(2700 + k); S <- true_sigma(p, q, "full"); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
  idx <- sample(p * q * n, round(eps * p * q * n)); ai <- arrayInd(idx, dim(X)); X[idx] <- X[idx] + 6 * sample(c(-1, 1), length(idx), TRUE) * sqrt(diag(S$Sr)[ai[, 1]] * diag(S$Sc)[ai[, 2]])
  out <- c(); for (v in c("default", "noflagquant", "notrim")) { f <- switch(v, default = cellmmcd(X), noflagquant = cellmmcd(X, flag_quant = NULL), notrim = cellmmcd(X, alpha_obs = 1))
    out <- c(out, sprintf("%s: KL %.3f trimmed %.2f eps_bar %.3f count_max %s", v, kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), f$trimmed, f$raw$eps_bar, format(f$raw$count_max))) }
  sprintf("  seed %d: %s", 2700 + k, paste(out, collapse = " | ")) }
res <- mclapply(16:18, one, mc.cores = 3); sink("output/cm_step27c_countrule.log"); cat("step 27c (block 11): count rule at eps = .25, 20x20, n = 100\n"); cat(unlist(res), sep = "\n"); sink()
cat(readLines("output/cm_step27c_countrule.log"), sep = "\n")
