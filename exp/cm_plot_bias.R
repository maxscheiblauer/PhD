# Block 3: plot of exp/cm_step6_bias.R (n * KL against n: bounded for a consistent estimator). Output: output/cm_step6_bias.png
R <- readRDS("output/cm_step6_bias.rds"); ns <- sort(unique(R[, "n"])); vars <- c("MLE", "default", "q999", "tail", "rw")
lab <- c("MLE", "raw cellMMCD, quant .99", "quant .999", "censored (tail)", "reweighted .99 -> .999 (default)")
m <- sapply(vars, function(v) sapply(ns, function(n) n * mean(R[R[, "n"] == n, v])))
png("output/cm_step6_bias.png", 800, 550, res = 110); par(mar = c(4.5, 4.5, 2, 1))
matplot(ns, m, type = "b", pch = 16, lty = 1, log = "x", xlab = "n (log scale)", ylab = "n x mean KL", main = "Clean 3x4, Gaussian: n x KL (bounded iff consistent)", col = 1:5)
legend("topleft", lab, col = 1:5, pch = 16, lty = 1, bty = "n", cex = 0.85); dev.off()
