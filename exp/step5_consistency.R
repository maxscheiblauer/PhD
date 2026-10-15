# Step 5 (Q6): is kappa(beta, p) right on clean data at small n? Oracle vs. plug-in.
# Run from project root: Rscript exp/step5_consistency.R
source("R/load.R")
set.seed(505)
ns <- c(10, 25, 50); betas <- c(0.5, 0.75, 0.9); nrep <- 12
methods <- c("oracle", "R1", "RJ")   # oracle: true Sigma^c fixed, C-steps on d only

t0 <- proc.time()[3]
res <- list(); qq <- list()
for (pq in c(3, 5)) {
  p <- q <- pq; S <- true_sigma(p, q, "full")
  for (n in ns) for (rep in seq_len(nrep)) {
    X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
    st <- make_init(X, "ff")
    for (be in betas) for (m in methods) for (kp in c(TRUE, FALSE)) {
      f <- tff(X, rule = if (m == "oracle") "R1" else m, beta = be, kappa = kp, init = st,
               Sc_fixed = if (m == "oracle") S$Sc, trace = FALSE)
      res[[length(res) + 1]] <- data.frame(p = p, n = n, beta = be, method = m, kappa = kp,
        rep = rep, status = f$status, ratio = tr_ratio(f$Sr, f$Sc, S$Sr, S$Sc))
      if (pq == 5 && n == 50 && be == 0.75 && kp && m != "oracle")
        qq[[m]] <- c(qq[[m]], as.vector(distances(X, f$Sr, f$Sc)$d))
    }
    if (pq == 5 && n == 50) qq$truth <- c(qq$truth, as.vector(distances(X, S$Sr, S$Sc)$d))
  }
}
res <- do.call(rbind, res)
elapsed <- proc.time()[3] - t0
summ <- aggregate(ratio ~ p + n + beta + method + kappa, res, mean)
summ$se <- aggregate(ratio ~ p + n + beta + method + kappa, res, function(x) sd(x) / sqrt(length(x)))$ratio
saveRDS(list(res = res, summ = summ, qq = qq, elapsed = elapsed), "output/step5.rds")

cat(sprintf("Step 5 done in %.1f s; non-converged runs: %d\n", elapsed, sum(res$status != "converged")))
cat("Mean of tr(Sigma^-1 Sigma_hat)/(pq) (1 = unbiased scale); kappa on:\n")
print(ftable(xtabs(ratio ~ p + method + beta + n, summ[summ$kappa, ])), digits = 3)
cat("\nkappa off:\n")
print(ftable(xtabs(ratio ~ p + method + beta + n, summ[!summ$kappa, ])), digits = 3)
cat("\nUntrimmed MLE reference (beta = 1) has E ratio = 1 exactly up to O(1/n) terms.\n")

png("output/step5_ratio.png", 1200, 450, res = 110)
par(mfrow = c(1, 2), mar = c(4, 4.5, 2.5, 0.5))
for (pp in c(3, 5)) {
  s <- summ[summ$p == pp & summ$kappa, ]
  plot(NA, xlim = c(8, 52), ylim = range(c(0.5, 1.2, s$ratio)), log = "x", xlab = "n",
       ylab = "mean tr(Sigma^-1 Sigma_hat)/(pq)", main = sprintf("p = q = %d, kappa on", pp))
  abline(h = 1, lty = 2)
  cols <- c(oracle = "black", R1 = "tomato", RJ = "steelblue"); ltys <- c(1, 2, 3)
  for (m in methods) for (bi in seq_along(betas)) {
    r <- s[s$method == m & s$beta == betas[bi], ]
    lines(r$n, r$ratio, col = cols[m], lty = ltys[bi], type = "b", pch = 19, cex = 0.6)
  }
  legend("bottomright", c(methods, paste("beta =", betas)), col = c(cols, rep("grey40", 3)),
         lty = c(1, 1, 1, ltys), bty = "n", cex = 0.75)
}
dev.off()

png("output/step5_qq.png", 1200, 420, res = 110)
par(mfrow = c(1, 3), mar = c(4, 4.5, 2.5, 0.5))
pr <- ppoints(500)
for (m in c("truth", "R1", "RJ")) {
  plot(qchisq(pr, 5), quantile(qq[[m]], pr), pch = 19, cex = 0.4,
       xlab = "chi^2_5 quantiles", ylab = "empirical d_ij quantiles",
       main = sprintf("%s (p = q = 5, n = 50, beta = .75)", if (m == "truth") "true parameters" else m))
  abline(0, 1, col = "tomato")
}
dev.off()
