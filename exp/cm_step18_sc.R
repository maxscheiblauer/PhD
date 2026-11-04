# Block 6: sensitivity curves (finite-sample influence functions) in the Kronecker model. Clean 3x3 full shape,
# n = 500, 3 seeds; one extra matrix z added: (a) one cell: z = 0 except z[1,1] = x * sd_11; (b) all cells
# shifted by x * sd (a casewise outlier). SC(x) = n * (Sigma_hat(X + z) - Sigma_hat(X)) with Sigma = Sc (x) Sr
# (scale-free comparison: det-identified), reported as the Frobenius norm and the (1,1) entry relative to
# Sigma_11. Methods: flip-flop MLE (SC ~ x^2, unbounded), cellMMCD default, reweighted hard quota.
# Output: output/cm_step18_sc.log / .rds, figure by exp/cm_plots_block6.R
source("R/load.R"); library(parallel); set.seed(18)
p <- q <- 3; n <- 500; S <- true_sigma(p, q, "full"); sdc <- sqrt(outer(diag(S$Sr), diag(S$Sc)))
xs <- c(0, 1, 2, 3, 4, 5, 6, 8, 10, 30)
kron <- function(f) { sc <- scale_ident(f$Sr, f$Sc, "det"); kronecker(sc$Sc, sc$Sr) }
fits <- function(X) list(MLE = kron(flipflop(X)), cellMMCD = kron(cellmmcd(X)), hard_rw = kron(cellmmcd(X, quant_obs = NULL)))
grid <- expand.grid(x = xs, type = c("cell", "matrix"), seed = 1:3, stringsAsFactors = FALSE)
base <- lapply(1:3, function(s) { set.seed(1800 + s); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); list(X = X, f = fits(X)) })
r <- mclapply(seq_len(nrow(grid)), function(g) { b <- base[[grid$seed[g]]]; z <- matrix(0, p, q)
  if (grid$type[g] == "cell") z[1, 1] <- grid$x[g] * sdc[1, 1] else z <- grid$x[g] * sdc
  X1 <- array(c(b$X, z), c(p, q, n + 1)); f1 <- fits(X1)
  out <- c(x = grid$x[g], seed = grid$seed[g], type = as.numeric(grid$type[g] == "matrix"))
  for (m in names(f1)) { D <- n * (f1[[m]] - b$f[[m]]); out[paste0(m, "_norm")] <- sqrt(sum(D^2)); out[paste0(m, "_11")] <- D[1, 1] / b$f[[m]][1, 1] }
  out }, mc.cores = 4)
R <- do.call(rbind, r); saveRDS(R, "output/cm_step18_sc.rds")
sink("output/cm_step18_sc.log")
cat("Sensitivity curves, 3x3 full shape, n = 500, 3 seeds (mean). n * (Sigma(X + z) - Sigma(X)): Frobenius norm | relative (1,1) entry.\n")
for (ty in 0:1) { cat(sprintf("\n-- %s --\n   x    MLE norm   cellMMCD   hard_rw  |  MLE (1,1)  cellMMCD   hard_rw\n", c("one cell z[1,1] = x sd", "whole matrix shifted by x sd")[ty + 1]))
  for (x in xs) { v <- colMeans(R[R[, "x"] == x & R[, "type"] == ty, , drop = FALSE])
    cat(sprintf("%4.0f  %9.2f  %9.2f  %9.2f  |  %9.2f  %9.2f  %9.2f\n", x, v["MLE_norm"], v["cellMMCD_norm"], v["hard_rw_norm"], v["MLE_11"], v["cellMMCD_11"], v["hard_rw_11"])) } }
sink()
