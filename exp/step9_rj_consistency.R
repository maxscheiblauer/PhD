# Step 9 (Q6 for RJ): is kappa(beta, p) on the product the right consistency factor for RJ?
# Run from project root: Rscript exp/step9_rj_consistency.R      (~3 min on 15 cores)
# Clean data, Sigma^c full. Per data set: RJ from 3 starts (ff, identity, elemental), keep lowest J*,
# and the RJ oracle (true Sigma^c fixed as A = Sc^(-1/2), only C-steps + Sr step).
# Theory (notes/rj_objective.md section 3): kappa(beta, p) Sigma_hat -> Sigma at the population level.
# Scale is measured two ways: ratio = tr(Sigma^-1 Sigma_hat)/(pq) (mean, the target of kappa) and the
# geometric ratio (|Sigma_hat| / |Sigma|)^(1/pq). c_emp = 1 / mean raw ratio is the factor the data ask for;
# compare it with kappa(beta, p) (RJ) and kappa(beta, pq) (MMCD-type, wrong unit dimension).
source("R/load.R")
library(parallel)
set.seed(909)
dims <- list(c(3, 3), c(5, 5), c(3, 6))
ns <- c(25, 50, 100, 200, 400, 1000); betas <- c(0.5, 0.75, 0.9); nrep <- 50

grid <- expand.grid(rep = seq_len(nrep), n = ns, d = seq_along(dims))
grid$seed <- sample.int(1e8, nrow(grid))

one <- function(g) {
  set.seed(g$seed)
  p <- dims[[g$d]][1]; q <- dims[[g$d]][2]; S <- true_sigma(p, q, "full")
  X <- rmatnorm(g$n, Sr = S$Sr, Sc = S$Sc)
  out <- list()
  for (be in betas) {
    m <- tff_multi(X, starts = c("ff", "identity", "elemental"), rule = "RJ", beta = be, kappa = FALSE)
    f <- m$fits[[m$best]]
    o <- tff(X, rule = "RJ", beta = be, kappa = FALSE, Sc_fixed = S$Sc, init = "ff", trace = FALSE)
    for (fit in list(list(method = "RJ", f = f), list(method = "oracle", f = o))) {
      ff <- fit$f
      out[[length(out) + 1]] <- data.frame(p = p, q = q, n = g$n, beta = be, method = fit$method,
        rep = g$rep, status = ff$status, n_conv = if (fit$method == "RJ") sum(sapply(m$fits, function(z)
          !is.null(z) && z$converged)) else NA,
        raw = tr_ratio(ff$Sr, ff$Sc, S$Sr, S$Sc),
        graw = exp((q * logdet(solve(S$Sr, ff$Sr)) + p * logdet(solve(S$Sc, ff$Sc))) / (p * q)),
        kl_shape = kl_shape(ff$Sr, ff$Sc, S$Sr, S$Sc))
    }
  }
  do.call(rbind, out)
}

t0 <- proc.time()[3]
cl <- makeCluster(min(15, detectCores() - 1))
clusterExport(cl, c("dims", "betas", "one"))
invisible(clusterEvalQ(cl, { source("R/load.R"); NULL }))
res <- parLapplyLB(cl, split(grid, seq_len(nrow(grid))), function(g)
  tryCatch(one(g), error = function(e) conditionMessage(e)))
stopCluster(cl)
err <- !sapply(res, is.data.frame)
if (any(err)) print(table(unlist(res[err])))
res <- do.call(rbind, res[!err])
elapsed <- proc.time()[3] - t0

res$kp <- mapply(kappa_mcd, res$beta, res$p)
res$kpq <- mapply(kappa_mcd, res$beta, res$p * res$q)
res$ratio <- res$raw * res$kp                 # RJ as implemented (kappa(beta, p) on the product)
res$gratio <- res$graw * res$kp
res$ratio_pq <- res$raw * res$kpq             # control: MMCD-type factor
summ <- aggregate(cbind(ratio, gratio, ratio_pq, raw, kl_shape) ~ method + p + q + beta + n, res, mean)
summ$se <- aggregate(ratio ~ method + p + q + beta + n, res, function(x) sd(x) / sqrt(length(x)))$ratio
summ$c_emp <- 1 / summ$raw
summ$kp <- mapply(kappa_mcd, summ$beta, summ$p)
summ$kpq <- mapply(kappa_mcd, summ$beta, summ$p * summ$q)
summ$nonconv <- aggregate(status ~ method + p + q + beta + n, res, function(s) mean(s != "converged"))$status
saveRDS(list(res = res, summ = summ, elapsed = elapsed), "output/step9.rds")

options(width = 200)
cat(sprintf("Step 9 done in %.1f s; %d fits; non-converged chosen fits: %d\n", elapsed, nrow(res),
            sum(res$status != "converged")))
cat("\nMean tr(Sigma^-1 Sigma_hat)/(pq) with kappa(beta, p) [SE]; 1 = consistent scale\n")
s <- summ[order(summ$method, summ$p, summ$q, summ$beta, summ$n), ]
s$cell <- sprintf("%.3f [%.3f]", s$ratio, s$se)
print(reshape(s[, c("method", "p", "q", "beta", "n", "cell")], idvar = c("method", "p", "q", "beta"),
              timevar = "n", direction = "wide"), row.names = FALSE)
cat("\nFactor the data ask for (c_emp = 1 / mean raw ratio) vs kappa(beta, p) and kappa(beta, pq), n = 1000:\n")
print(s[s$n == 1000, c("method", "p", "q", "beta", "c_emp", "kp", "kpq", "gratio", "kl_shape")],
      row.names = FALSE, digits = 3)
cat("\nShare of non-converged chosen fits by (p, q, beta, n), RJ only:\n")
print(ftable(xtabs(nonconv ~ p + q + beta + n, s[s$method == "RJ", ])), digits = 2)

png("output/step9_ratio.png", 1500, 450, res = 110)
par(mfrow = c(1, 3), mar = c(4, 4.5, 2.5, 0.5))
for (dd in dims) {
  s1 <- s[s$p == dd[1] & s$q == dd[2], ]
  plot(NA, xlim = range(ns), ylim = range(c(0.85, 1.1, s1$ratio)), log = "x", xlab = "n",
       ylab = "mean tr(Sigma^-1 kappa Sigma_hat)/(pq)", main = sprintf("RJ, p = %d, q = %d", dd[1], dd[2]))
  abline(h = 1, lty = 2)
  for (bi in seq_along(betas)) for (m in c("RJ", "oracle")) {
    r <- s1[s1$method == m & s1$beta == betas[bi], ]
    lines(r$n, r$ratio, col = bi, lty = if (m == "RJ") 1 else 3, type = "b", pch = 19, cex = 0.6)
  }
  legend("bottomright", c(paste("beta =", betas), "RJ", "oracle (true Sc)"),
         col = c(seq_along(betas), "grey40", "grey40"), lty = c(1, 1, 1, 1, 3), bty = "n", cex = 0.75)
}
dev.off()
