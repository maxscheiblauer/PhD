# Step 10 (Q8 for RJ): when does RJ fail to converge at small n, and why?
# Run from project root: Rscript exp/step10_rj_existence.R      (~4 min on 15 cores)
# Step 5 found RJ hitting maxit at n = 10, beta = .5, p = q = 5 (h = 5 = p kept matrices per column index).
# Here: clean data, Sigma^c full, p in {3, 5}, q in {3, 5, 8}, beta in {.5, .75}, n chosen so that the
# per-column quota h = round(beta n) runs from min(p, q) - 1 to max(p, q) + 4. Two starts (ff, identity),
# maxit 500, 15 reps. A failed linear solve is recorded as status "singular" (the same degeneracy, caught
# by R before tff's own checks).
# Per run: status, iterations, and diagnostics to tell the failure modes apart:
#   drift   = J*(end) - J*(end - 50): large negative => J* heads to -Inf (no minimiser exists)
#   cond_Sr / cond_Sc = condition numbers of the final factors (blow-up => degenerate solution)
#   step    = relative change of the product in the last iteration (small => slow, not divergent)
source("R/load.R")
library(parallel)
set.seed(1010)
nrep <- 15; betas <- c(0.5, 0.75); ps <- c(3, 5); qs <- c(3, 5, 8)

grid <- do.call(rbind, lapply(ps, function(p) do.call(rbind, lapply(qs, function(q)
  do.call(rbind, lapply(betas, function(be) {
    n <- 2:60; h <- round(be * n)
    keep <- !duplicated(h) & h >= min(p, q) - 1 & h <= max(p, q) + 4
    data.frame(p = p, q = q, beta = be, n = n[keep], h = h[keep])
  }))))))
grid <- grid[rep(seq_len(nrow(grid)), each = nrep), ]
grid$rep <- ave(seq_len(nrow(grid)), grid$p, grid$q, grid$beta, grid$n, FUN = seq_along)
grid$seed <- sample.int(1e8, nrow(grid))

cond <- function(S) { ev <- eigen(S, TRUE, TRUE)$values; max(ev) / min(ev) }
one <- function(g) {
  set.seed(g$seed)
  q <- g$q; S <- true_sigma(g$p, q, "full")
  X <- rmatnorm(g$n, Sr = S$Sr, Sc = S$Sc)
  do.call(rbind, lapply(c("ff", "identity"), function(st) {
    f <- tryCatch(tff(X, rule = "RJ", beta = g$beta, init = st, maxit = 500, trace = FALSE),
                  error = function(e) NULL)
    if (is.null(f)) return(data.frame(g[c("p", "q", "beta", "n", "h", "rep")], start = st, status = "singular",
      iter = NA, drift = NA, cond_Sr = NA, cond_Sc = NA, step = NA, kl_shape = NA))
    k <- length(f$J)
    data.frame(g[c("p", "q", "beta", "n", "h", "rep")], start = st, status = f$status, iter = f$iter,
      drift = if (k > 50) f$J[k] - f$J[k - 50] else NA,
      cond_Sr = cond(f$Sr), cond_Sc = cond(f$Sc), step = tail(f$delta, 1),
      kl_shape = kl_shape(f$Sr, f$Sc, S$Sr, S$Sc))
  }))
}

t0 <- proc.time()[3]
cl <- makeCluster(min(15, detectCores() - 1))
clusterExport(cl, c("cond", "one"))
invisible(clusterEvalQ(cl, { source("R/load.R"); NULL }))
res <- do.call(rbind, parLapplyLB(cl, split(grid, seq_len(nrow(grid))), one))
stopCluster(cl)
elapsed <- proc.time()[3] - t0
saveRDS(list(res = res, elapsed = elapsed), "output/step10.rds")

options(width = 200)
cat(sprintf("Step 10 done in %.1f s; %d runs\n", elapsed, nrow(res)))
res$ok <- res$status == "converged"
cat("\nStatus counts (all settings, both starts):\n")
print(table(res$status))
cat("\nShare converged by (p, q, beta) and h (kept matrices per column index):\n")
print(ftable(xtabs(ok ~ p + q + beta + h, res) / xtabs(~ p + q + beta + h, res)), digits = 2)
# smallest h from which every run converges (for all larger h too)
hmin <- do.call(rbind, lapply(split(aggregate(ok ~ p + q + beta + h, res, mean), ~ p + q + beta), function(d) {
  d <- d[order(d$h), ]; bad <- which(d$ok < 1)
  data.frame(p = d$p[1], q = d$q[1], beta = d$beta[1],
             h_all_conv = if (length(bad)) d$h[max(bad)] + 1 else min(d$h))
}))
hmin$p_plus_1 <- hmin$p + 1; hmin$q_plus_1 <- hmin$q + 1
cat("\nSmallest h with 100 % convergence from there on, vs p + 1 and q + 1:\n")
print(hmin, row.names = FALSE)
cat("\nNon-converged runs: diagnostics (median) by status\n")
nc <- res[!res$ok, ]
print(aggregate(cbind(iter, drift, l10_step = log10(step)) ~ status, nc, median, na.action = na.pass), digits = 3)
cat("\nConverged runs: median iterations and shape-KL by h - q\n")
cv <- res[res$ok, ]; cv$hq <- cv$h - cv$q
print(aggregate(cbind(iter, kl_shape) ~ hq, cv, median), digits = 3)
