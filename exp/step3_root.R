# Step 3 (Q3): dependence of the estimate on the choice of root (sym, chol, 5 random rotations).
# Run from project root: Rscript exp/step3_root.R
source("R/load.R")
set.seed(303)
p <- q <- 5; n <- 25; S <- true_sigma(p, q, "full")   # non-diagonal Sigma^c, otherwise sym = chol
conds <- c("clean", "col", "wcol"); rules <- c("R0", "R1", "R3", "R2", "RJ")
betas <- c(1, 0.75); nds <- 8
variants <- c("sym", "chol", rep("randrot", 5))

t0 <- proc.time()[3]
res <- list()
for (cd in conds) for (ds in seq_len(nds)) {
  X <- gen_data(n, S$Sr, S$Sc, cd, eps = if (cd == "clean") 0 else 0.1)$X
  st <- make_init(X, "ff")
  rots <- lapply(1:5, function(k) list(c = rorth(q), r = rorth(p)))
  for (rl in rules) for (be in betas) {
    fits <- lapply(seq_along(variants), function(v)
      tff(X, rule = rl, beta = be, root_type = variants[v], init = st, trace = FALSE,
          rot = if (v > 2) rots[[v - 2]] else NULL, tol = 1e-10))
    pk <- combn(length(fits), 2)
    sk <- apply(pk, 2, function(ij) skl(fits[[ij[1]]], fits[[ij[2]]]))
    kl <- sapply(fits, function(f) kl_mn(f$Sr, f$Sc, S$Sr, S$Sc))
    res[[length(res) + 1]] <- data.frame(
      cond = cd, ds = ds, rule = rl, beta = be,
      sKL_max = max(sk), sKL_med = median(sk),
      sKL_sym_chol = sk[1],                         # pair (sym, chol)
      n_distinct_a = length(unique(sapply(fits, function(f) paste(as.integer(f$a), collapse = "")))),
      KL_truth_med = median(kl), KL_truth_range = diff(range(kl)),
      n_conv = sum(sapply(fits, `[[`, "converged")))
  }
}
res <- do.call(rbind, res)
elapsed <- proc.time()[3] - t0
saveRDS(list(res = res, elapsed = elapsed), "output/step3.rds")

summ <- aggregate(cbind(sKL_max, sKL_med, sKL_sym_chol, n_distinct_a, KL_truth_med, KL_truth_range) ~
                    beta + rule + cond, res, median)
cat(sprintf("Step 3 done in %.1f s\n", elapsed))
print(summ[order(summ$beta, summ$cond, summ$rule), ], digits = 3, row.names = FALSE)
cat("\nRatio median(max pairwise sKL across roots) / median(KL to truth), beta = 0.75:\n")
r75 <- summ[summ$beta == 0.75, ]
print(xtabs(I(sKL_max / KL_truth_med) ~ rule + cond, r75), digits = 2)

png("output/step3_root_spread.png", 1300, 450, res = 110)
par(mfrow = c(1, 3), mar = c(4, 4.5, 2.5, 0.5))
for (cd in conds) {
  r <- res[res$cond == cd & res$beta == 0.75, ]
  r$rule <- factor(r$rule, rules)
  boxplot(log10(sKL_max) ~ rule, r, ylim = c(-3, 1.6), col = "grey85",
          ylab = "log10 sKL", main = sprintf("%s, beta = .75", cd))
  boxplot(log10(KL_truth_med) ~ rule, r, add = TRUE, at = 1:5 + 0.3, boxwex = 0.3,
          col = "tomato", axes = FALSE)
  text(4, -2.5, "R2: spread ~1e-15\n(off scale)", cex = 0.75)
  legend("bottomleft", c("max pairwise sKL across 7 roots", "KL to truth (median over roots)"),
         fill = c("grey85", "tomato"), bty = "n", cex = 0.8)
}
dev.off()
