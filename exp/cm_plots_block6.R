# Block 6 figure: sensitivity curves in the Kronecker model (cm_step18_sc.rds). Output: output/cm_fig_sc.png
source("R/load.R")
pal <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100", magenta = "#e87ba4", green = "#008300")
ink <- "#1a1a19"; grid_col <- "#e4e3df"
R <- readRDS("output/cm_step18_sc.rds"); xs <- sort(unique(R[, "x"]))
agg <- function(col, ty) sapply(xs, function(x) mean(R[R[, "x"] == x & R[, "type"] == ty, col]))
png("output/cm_fig_sc.png", 1000, 480, res = 120); par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.5, 1), col.axis = ink, col.lab = ink, col.main = ink, cex.main = 1, font.main = 1, las = 1)
for (ty in 0:1) {
  plot(NA, xlim = c(0, 10), ylim = c(1, 1000), log = "y", xlab = if (ty == 0) "x: one cell z[1,1] = x sd" else "x: all cells shifted by x sd", ylab = "|| n (Sigma(X + z) - Sigma(X)) ||_F",
       main = if (ty == 0) "Sensitivity curve, one outlying cell (3x3, n = 500)" else "Sensitivity curve, one shifted matrix", bty = "l", yaxt = "n")
  yt <- c(1, 3, 10, 30, 100, 300, 1000); axis(2, yt, yt, col = grid_col); abline(h = yt, col = grid_col, lwd = 0.8)
  for (s in list(list("MLE_norm", pal["yellow"], "MLE"), list("cellMMCD_norm", pal["blue"], "cellMMCD (default)"), list("hard_rw_norm", pal["orange"], "reweighted hard quota"))) {
    y <- agg(s[[1]], ty); lines(xs[xs <= 10], y[xs <= 10], col = s[[2]], lwd = 2); points(xs[xs <= 10], y[xs <= 10], col = s[[2]], pch = 16, cex = 0.9) }
  if (ty == 0) legend("topleft", c("MLE (flip-flop): ~ x^2", "cellMMCD (default)", "reweighted hard quota"), col = pal[c("yellow", "blue", "orange")], lwd = 2, bty = "n", cex = 0.85, text.col = ink)
}
dev.off()
