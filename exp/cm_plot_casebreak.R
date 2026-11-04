# Block 8: figure for the casewise breakdown of the three observation rules (block 7, step 21; STATE-7 item 3).
# Reads output/cm_step21_casebreak.rds (rerun in block 8, same seeds), writes output/cm_fig_casebreak.png.
res <- readRDS("output/cm_step21_casebreak.rds"); rules <- c(soft = "soft threshold (default)", hard_raw = "hard quota, raw", hard_rw = "hard quota, reweighted")
png("output/cm_fig_casebreak.png", width = 1600, height = 700, res = 160); par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
for (ty in c("shift", "struct")) { sel <- Filter(function(z) z$grid$type == ty, res); eps <- sapply(sel, function(z) z$grid$eps)
  plot(NA, xlim = range(eps), ylim = c(0.05, 20), log = "y", xlab = "fraction of casewise outlier matrices", ylab = "median KL (5 reps)",
       main = sprintf("%s matrices, 5x5, n = 100", if (ty == "shift") "shifted" else "structural"))
  abline(v = 0.25, lty = 3, col = "grey50"); for (j in seq_along(rules)) { v <- names(rules)[j]
    kl <- sapply(sel, function(z) median(z$r[, paste0(v, ".kl")], na.rm = TRUE)); lines(eps, kl, type = "b", pch = c(19, 17, 15)[j], col = c("black", "firebrick", "steelblue")[j], lwd = 2) }
  if (ty == "shift") legend("topleft", legend = rules, col = c("black", "firebrick", "steelblue"), pch = c(19, 17, 15), lwd = 2, bty = "n", cex = 0.85) }
dev.off(); cat("written output/cm_fig_casebreak.png\n")
