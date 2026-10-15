# Step 6 (Q7): column-trimmed rules vs. MMCD under five contamination types.
# Run from project root:
#   Rscript exp/step6_mmcd.R full   # Sigma^c = AR(0.7): main KL comparison
#   Rscript exp/step6_mmcd.R diag   # Sigma^c diagonal: precision/recall of flagged columns is meaningful
source("R/load.R")
sc <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(sc)) sc <- "full"
set.seed(if (sc == "full") 606 else 616)
p <- q <- 5; n <- 25; beta <- 0.75; nrep <- 8       # matched breakdown: beta = alpha = 0.75
S <- true_sigma(p, q, sc)
types <- c("obs", "col", "row", "cell", "wcol"); epss <- c(0.1, 0.2)
rules <- c("R0", "R1", "R2", "R3", "RJ")
starts <- c("identity", "ff", "elemental")
conds <- rbind(data.frame(type = "clean", eps = 0), expand.grid(type = types, eps = epss, stringsAsFactors = FALSE))

t0 <- proc.time()[3]
res <- list()
for (ci in seq_len(nrow(conds))) for (rep in seq_len(nrep)) {
  dat <- gen_data(n, S$Sr, S$Sc, if (conds$eps[ci] == 0) "clean" else conds$type[ci], conds$eps[ci])
  X <- dat$X
  add <- function(method, Sr, Sc, a) {
    pr <- prec_recall(a, dat$truth)
    res[[length(res) + 1]] <<- data.frame(type = conds$type[ci], eps = conds$eps[ci], rep = rep,
      method = method, KL = kl_mn(Sr, Sc, S$Sr, S$Sc), KL_shape = kl_shape(Sr, Sc, S$Sr, S$Sc),
      precision = pr[["precision"]], recall = pr[["recall"]])
  }
  ff <- flipflop(X); add("FF", ff$Sr, ff$Sc, matrix(TRUE, n, q))
  for (rl in rules) {
    m <- tff_multi(X, starts = starts, rule = rl, beta = beta)
    f <- m$fits[[m$best]]
    add(rl, f$Sr, f$Sc, f$a)
  }
  if (has_robustmatrix) {
    mm <- robustmatrix::mmcd(X, alpha = beta, adapt_alpha = FALSE, reweight = FALSE)
    keep <- rank(mm$md, ties.method = "first") <= floor(beta * n)
    add("MMCD", mm$cov_row, mm$cov_col, matrix(keep, n, q))
  }
}
res <- do.call(rbind, res)
elapsed <- proc.time()[3] - t0
methods <- c("FF", rules, if (has_robustmatrix) "MMCD")
res$method <- factor(res$method, methods)

fmt <- function(x) sprintf("%.2f (%.2f)", median(x), IQR(x))
tab <- aggregate(KL ~ type + eps + method, res, fmt)
tab_w <- reshape(tab, idvar = c("type", "eps"), timevar = "method", direction = "wide")
names(tab_w) <- sub("KL\\.", "", names(tab_w))
tab_w <- tab_w[order(match(tab_w$type, c("clean", types)), tab_w$eps), ]
pr <- aggregate(cbind(precision, recall) ~ type + eps + method, res[res$eps > 0, ],
                function(x) mean(x, na.rm = TRUE), na.action = na.pass)
saveRDS(list(res = res, tab = tab_w, pr = pr, sc = sc, elapsed = elapsed), sprintf("output/step6_%s.rds", sc))

png(sprintf("output/step6_%s_boxplots.png", sc), 1800, 420, res = 110)
par(mfrow = c(1, 5), mar = c(5.5, 4, 2.5, 0.5))
cols <- c("grey60", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#0072B2", "black")[seq_along(methods)]
for (tp in types) {
  r <- res[res$type %in% c(tp, "clean"), ]
  boxplot(log10(KL) ~ eps + method, r, col = rep(cols, each = 3), las = 2, cex.axis = 0.7,
          ylab = "log10 KL", xlab = "", main = sprintf("%s (Sigma^c %s)", tp, sc), pch = 20, cex = 0.5,
          names = rep(c("0", "", ""), length(methods)), at = rep(seq_along(methods) * 4, each = 3) + 0:2 - 4)
  axis(1, at = seq_along(methods) * 4 - 3, labels = methods, las = 2, tick = FALSE, line = 1, cex.axis = 0.75)
  mtext("eps = 0 / .1 / .2 within each method", 3, 0.1, cex = 0.55)
}
dev.off()

cat(sprintf("Step 6 (%s) done in %.1f s\n", sc, elapsed))
cat("Median KL (IQR), beta = alpha = 0.75, n = 25, p = q = 5:\n")
print(tab_w, row.names = FALSE)
if (sc == "diag") {
  cat("\nMean precision / recall of flagged columns (only meaningful for diagonal Sigma^c):\n")
  print(reshape(transform(pr, pr = sprintf("%.2f/%.2f", precision, recall))[, c("type", "eps", "method", "pr")],
                idvar = c("type", "eps"), timevar = "method", direction = "wide"), row.names = FALSE)
}
