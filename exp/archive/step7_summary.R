# Step 7 summary: tables and plots from output/step7.rds.
# Run from project root: Rscript exp/archive/step7_summary.R
source("R/load.R")
res <- readRDS("output/step7.rds")$res
res$nonconv <- res$n_cycle + res$n_maxit
pct <- function(x) sprintf("%.1f", 100 * mean(x, na.rm = TRUE))
cat(sprintf("%d data sets, %d R1 runs\n\n", nrow(res), sum(res$n_starts)))

## 1. convergence
cat("Share of R1 RUNS that converge / cycle / hit maxit (1000) / singular, by setting [%]:\n")
conv <- function(d) c(runs = sum(d$n_starts), conv = 100 * sum(d$n_conv) / sum(d$n_starts),
  cycle = 100 * sum(d$n_cycle) / sum(d$n_starts), maxit = 100 * sum(d$n_maxit) / sum(d$n_starts),
  sing = 100 * sum(d$n_sing) / sum(d$n_starts),
  ds_any_nonconv = 100 * mean(d$nonconv > 0), ds_best_nonconv = 100 * mean(d$best_status != "converged"),
  iter_med = median(d$iter_med))
tab_by <- function(f) { s <- split(res, f, drop = TRUE); round(t(sapply(s, conv)), 1) }
print(tab_by(res$beta)); cat("\n"); print(tab_by(res$pq)); cat("\n"); print(tab_by(res$n)); cat("\n")
print(tab_by(res$type)); cat("\n")
tab_conv <- tab_by(list(res$pq, res$n, res$beta))
cat("By (p x q, n, beta), only settings with any non-convergence:\n")
print(tab_conv[tab_conv[, "conv"] < 100, ])

## 2. accuracy vs MLE and MMCD (KL_shape: scale-free)
cat("\nMedian KL_shape to truth by type x beta (all settings pooled):\n")
acc <- aggregate(cbind(MLE = KLs_MLE, MMCD = KLs_MMCD, R1 = KLs_R1) ~ type + beta, res, median)
print(acc, digits = 3)
cat("\nShare of data sets where R1 (best start) beats MMCD / MLE in KL_shape, by type [%]:\n")
win <- aggregate(cbind(R1_vs_MMCD = KLs_R1 < KLs_MMCD, R1_vs_MLE = KLs_R1 < KLs_MLE) ~ type + beta, res,
                 function(x) round(100 * mean(x)))
print(reshape(win, idvar = "type", timevar = "beta", direction = "wide"), row.names = FALSE)
cat("\nMedian KL_shape, p x q x n slices, beta = 0.75:\n")
print(aggregate(cbind(MLE = KLs_MLE, MMCD = KLs_MMCD, R1 = KLs_R1) ~ type + pq + n,
                res[res$beta == 0.75 & res$type %in% c("clean", "col", "row", "obs"), ], median), digits = 3)
cat("\nDoes non-convergence hurt? Median KL_shape ratio R1 / MMCD by status of the best start:\n")
print(aggregate(cbind(ratio = KLs_R1 / KLs_MMCD) ~ best_status, res, median), digits = 3)
print(table(res$best_status))
cat(sprintf("Non-converged best start: median step size over the last 50 iterations %.3g (shape-KL), median KL_shape(R1, truth) %.3g\n",
            median(res$amp_best[res$best_status != "converged"]), median(res$KLs_R1[res$best_status != "converged"])))

## 3. Q4 diagnostics (ff start)
cat("\nLocal contraction rate at converged fixed points (ff start): quantiles of rho (R1) and rho_ff (MLE flip-flop)\n")
print(rbind(R1 = quantile(res$rho, c(.1, .5, .9, .99, 1), na.rm = TRUE),
            FF = quantile(res$rho_ff, c(.1, .5, .9, .99, 1), na.rm = TRUE)), digits = 3)
cat(sprintf("rho > rho_ff in %s %% of fixed points; rho >= 0.9 in %s %%\n",
            pct(res$rho > res$rho_ff), pct(res$rho >= 0.9)))
print(aggregate(cbind(rho = rho, rho_ff = rho_ff) ~ beta, res, median), digits = 3)
cat("\nShare of ff-start runs in which the criterion increases at least once [%]:\n")
print(c(Lr = pct(res$Lr_inc > 0), Lc = pct(res$Lc_inc > 0), J = pct(res$J_inc > 0), gap = pct(res$gap_inc > 0)))
cat("\nLate iterations (step size 1e-9..1e-4): median |cross pull| / |own gain|\n")
print(rbind(Lr = quantile(res$late_ratio_r, c(.1, .5, .9), na.rm = TRUE),
            Lc = quantile(res$late_ratio_c, c(.1, .5, .9), na.rm = TRUE)), digits = 3)
cat(sprintf("share of runs whose late ratio > 1: Lr %s %%, Lc %s %%; late half-steps with J~ increasing: %s %%\n",
            pct(res$late_ratio_r > 1), pct(res$late_ratio_c > 1), pct(res$late_Jinc)))

saveRDS(list(conv_beta = tab_by(res$beta), conv_type = tab_by(res$type), conv_grid = tab_conv,
             acc = acc, win = win), "output/step7_summary.rds")

## plots
okabe <- c(MLE = "#999999", MMCD = "#0072B2", R1 = "#D55E00")
png("output/step7_convergence.png", 1500, 480, res = 110)
par(mfrow = c(1, 3), mar = c(5, 4.2, 2.5, 0.8))
for (v in c("beta", "pq", "n")) {
  s <- tab_by(res[[v]])
  barplot(t(s[, c("cycle", "maxit")]), beside = FALSE, col = c("#D55E00", "#E69F00"), border = NA,
          ylab = "% of R1 runs not converged", xlab = v, main = paste("non-convergence by", v),
          ylim = c(0, max(s[, "cycle"] + s[, "maxit"], 1) * 1.25))
  legend("topright", c("exact cycle", "maxit (1000)"), fill = c("#D55E00", "#E69F00"), border = NA, bty = "n")
}
dev.off()

png("output/step7_accuracy.png", 1700, 480, res = 110)
types <- c("clean", "t3", "obs", "col", "row", "cell", "wcol")
par(mfrow = c(1, 3), mar = c(5, 4.2, 2.5, 0.8))
for (be in c(0.5, 0.75, 0.9)) {
  d <- res[res$beta == be, ]
  m <- sapply(types, function(tp) sapply(c("KLs_MLE", "KLs_MMCD", "KLs_R1"), function(v) median(d[d$type == tp, v], na.rm = TRUE)))
  barplot(log10(m), beside = TRUE, col = okabe, border = NA, names.arg = types, las = 2,
          ylab = "log10 median KL_shape to truth", main = sprintf("beta = %.2f (all p, q, n pooled)", be))
  legend("topleft", names(okabe), fill = okabe, border = NA, bty = "n")
}
dev.off()

png("output/step7_rho.png", 1200, 520, res = 110)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 0.8))
cols <- c("0.5" = "#D55E00", "0.75" = "#0072B2", "0.9" = "#009E73")
plot(res$rho_ff, res$rho, pch = 20, cex = 0.5, col = adjustcolor(cols[as.character(res$beta)], 0.5),
     xlim = c(0, 1), ylim = c(0, 1.05), xlab = "rho, untrimmed flip-flop at the MLE",
     ylab = "rho, R1 at its fixed point", main = "local contraction (< 1: attracting)")
abline(0, 1, lty = 2, col = "grey50"); abline(h = 1, col = "grey50")
legend("bottomright", paste("beta =", names(cols)), col = cols, pch = 20, bty = "n")
r <- res$late_ratio_r[is.finite(res$late_ratio_r)]; cc <- res$late_ratio_c[is.finite(res$late_ratio_c)]
br <- seq(floor(log10(min(r, cc))), ceiling(log10(max(r, cc))), by = 0.5)
hist(log10(cc), breaks = br, col = adjustcolor("#D55E00", 0.6), border = "white",
     main = "late iterations: |cross pull| / |own gain|", xlab = "log10 ratio (> 0: the pull dominates)")
hist(log10(r), breaks = br, col = adjustcolor("#0072B2", 0.6), border = "white", add = TRUE)
abline(v = 0, lty = 2)
legend("topright", c("Lr~ (pulled by column step)", "Lc~ (pulled by row step)"),
       fill = adjustcolor(c("#0072B2", "#D55E00"), 0.6), border = NA, bty = "n")
dev.off()
