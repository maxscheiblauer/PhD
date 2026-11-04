# Block 5: figures for the write-up. Reads the saved experiment results (no new simulation) and the proxy IF.
# Output: output/cm_fig_breakdown.png, output/cm_fig_timing.png, output/cm_fig_if.png, output/cm_fig_mixed.png
source("R/load.R")
pal <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100", magenta = "#e87ba4", green = "#008300")
ink <- "#1a1a19"; grid_col <- "#e4e3df"
setup <- function(file, w = 900, h = 560) { png(file, w, h, res = 120); par(mar = c(4.2, 4.5, 2.5, 1), col.axis = ink, col.lab = ink, col.main = ink, cex.main = 1, font.main = 1, las = 1) }
hgrid <- function(at) abline(h = at, col = grid_col, lwd = 0.8)

## 1. Breakdown (cm_step13): median KL against the contamination fraction, log scale
bd <- readRDS("output/cm_step13_breakdown.rds")
tab <- do.call(rbind, lapply(bd, function(b) data.frame(b$grid, kl = median(b$r[, "kl"], na.rm = TRUE))))
series <- list(list(mode = "cell", mag = 4, lab = "cells at 4 sd", col = pal["blue"]), list(mode = "cell", mag = 10, lab = "cells at 10 sd", col = pal["orange"]),
               list(mode = "cell", mag = 1000, lab = "cells at 1000 sd", col = pal["aqua"]), list(mode = "case", mag = 5, lab = "matrices shifted 5 sd", col = pal["yellow"]))
setup("output/cm_fig_breakdown.png")
plot(NA, xlim = c(0.1, 0.5), ylim = c(0.1, 1e7), log = "y", xlab = "contaminated fraction (cells, or matrices)", ylab = "median KL to the truth (log scale)",
     main = "cellMMCD breakdown, p = q = 5, n = 100 (cm_step13): quotas at 25 %", yaxt = "n", bty = "l")
yt <- 10^(-1:7); axis(2, yt, c("0.1", "1", "10", "100", "1e3", "1e4", "1e5", "1e6", "1e7"), col = grid_col); hgrid(yt); abline(v = 0.25, col = grid_col, lwd = 2)
for (s in series) { d <- tab[tab$mode == s$mod & tab$mag == s$mag, ]; lines(d$eps, d$kl, col = s$col, lwd = 2); points(d$eps, d$kl, col = s$col, pch = 16, cex = 1.1) }
text(0.25, 0.13, "quota 1 - alpha = 0.25", pos = 4, cex = 0.8, col = ink)
legend("topleft", sapply(series, `[[`, "lab"), col = sapply(series, `[[`, "col"), lwd = 2, pch = 16, bty = "n", cex = 0.85, text.col = ink)
dev.off()

## 2. Timing (cm_step10_timing.log and _Rengine.log, n = 100, 5 % cells at 6 sd, one core)
pq <- c(5, 10, 15, 20); t_cpp <- c(0.3, 0.4, 2.2, 8.6); t_R <- c(0.8, 2.7, 12.4, 38.8); t_mmcd <- c(0.2, 0.8, 2.2, 4.7); t_mle <- c(0.03, 0.04, 0.10, 0.21)
setup("output/cm_fig_timing.png")
plot(NA, xlim = c(5, 20), ylim = c(0.02, 50), log = "y", xlab = "p = q", ylab = "seconds (one core, n = 100)", main = "Run time (cm_step10): cellMMCD C++ cell step vs R loops, MMCD, MLE", bty = "l", yaxt = "n")
yt <- c(0.02, 0.1, 0.5, 2, 10, 50); axis(2, yt, yt, col = grid_col); hgrid(yt)
for (s in list(list(t_cpp, "cellMMCD (C++ cell step)", pal["blue"]), list(t_R, "cellMMCD (R loops)", pal["orange"]), list(t_mmcd, "MMCD", pal["aqua"]), list(t_mle, "MLE (flip-flop)", pal["yellow"]))) {
  lines(pq, s[[1]], col = s[[3]], lwd = 2); points(pq, s[[1]], col = s[[3]], pch = 16, cex = 1.1); text(20, s[[1]][4], s[[2]], pos = 2, cex = 0.8, col = ink, offset = 0.7) }
dev.off()

## 3. Influence functions in the 1-D proxy (METHOD 7a): raw, reweighted, MLE
x <- seq(-5, 5, by = 0.01); ifs <- sapply(x, cm_proxy_if)
setup("output/cm_fig_if.png")
plot(NA, xlim = c(-5, 5), ylim = c(-2, 12), xlab = "x (standardized cell residual)", ylab = "IF of the scale functional", main = "Influence functions in the 1-D proxy (METHOD 7a): bounded, ~0 beyond the cutoff", bty = "l")
hgrid(seq(-2, 12, 2)); abline(h = 0, col = ink, lwd = 0.8)
lines(x, x^2 - 1, col = pal["yellow"], lwd = 2); lines(x, ifs["raw", ], col = pal["blue"], lwd = 2); lines(x, ifs["rw", ], col = pal["orange"], lwd = 2)
legend("top", c("MLE: x^2 - 1 (unbounded)", "raw cellMMCD (quant .99)", "reweighted (quant .999, default)"), col = pal[c("yellow", "blue", "orange")], lwd = 2, bty = "n", cex = 0.85, text.col = ink)
dev.off()

## 4. Casewise and mixed benchmark (cm_step4_mixed.rds): default vs hard quota vs MMCD, log scale
mx <- readRDS("output/cm_step4_mixed.rds")
lab <- sapply(mx, function(m) { g <- m$grid; paste0(if (g$type == "none") "" else paste0(g$type, " ", 100 * g$eps_case, "%"), if (g$eps_cell > 0) paste0(if (g$type != "none") " + " else "", 100 * g$eps_cell, "% cells") else "", if (g$type == "none" && g$eps_cell == 0) "clean" else "") })
kl <- sapply(mx, function(m) sapply(c("cellMMCD_soft", "cellMMCD_hard", "MMCD"), function(v) median(sapply(m$reps, function(r) r[[v]]["kl"]), na.rm = TRUE)))
setup("output/cm_fig_mixed.png", 1000, 560); par(mar = c(7.5, 4.5, 2.5, 1))
bp <- barplot(kl, beside = TRUE, log = "y", ylim = c(0.1, 100), col = pal[c("blue", "orange", "aqua")], border = NA, names.arg = rep("", ncol(kl)), ylab = "median KL (log scale)",
              main = "Casewise and mixed contamination, p = q = 5, n = 100, 10 reps (cm_step4_mixed)", yaxt = "n", space = c(0.15, 1))
yt <- c(0.1, 0.3, 1, 3, 10, 30, 100); axis(2, yt, yt, col = grid_col); hgrid(yt)
text(colMeans(bp), 0.075, lab, srt = 35, adj = 1, xpd = TRUE, cex = 0.75, col = ink)
legend("topleft", c("cellMMCD default (soft observation step, count rule)", "cellMMCD hard quota", "MMCD"), fill = pal[c("blue", "orange", "aqua")], border = NA, bty = "n", cex = 0.85, text.col = ink)
dev.off(); cat("done\n")
