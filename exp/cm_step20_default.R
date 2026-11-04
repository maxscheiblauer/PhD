# Block 7, step 20: choice of the default, soft threshold vs reweighted hard quota, on a larger grid than cm_step16.
# p = q in {5, 8}, n in {100, 200}, 10 reps per cell of the grid, same contamination code as cm_step16 (6-sd cells,
# shift +5, scale x3, struct = AR(-0.5) x AR(-0.7)). MMCD and the MLE for reference. Seed 20. Medians of KL, shape KL,
# trace ratio, trimmed bad / good matrices, cell recall. Output: output/cm_step20_default.log / .rds
source("R/load.R"); library(parallel); suppressMessages(library(robustmatrix)); set.seed(20)
contaminate <- function(X, S, type, eps_case = 0.1, eps_cell = 0) {
  p <- dim(X)[1]; q <- dim(X)[2]; n <- dim(X)[3]
  bad_case <- rep(FALSE, n); bad_cell <- array(FALSE, c(p, q, n))
  if (eps_cell > 0) { m <- round(eps_cell * p * q * n); idx <- sample(p * q * n, m); bad_cell[idx] <- TRUE
    sdc <- array(sqrt(outer(diag(S$Sr), diag(S$Sc))), c(p, q, n)); X[idx] <- X[idx] + 6 * sdc[idx] * sample(c(-1, 1), m, TRUE) }
  if (eps_case > 0) { idx <- sample(n, round(eps_case * n)); bad_case[idx] <- TRUE
    if (type == "shift") X[, , idx] <- X[, , idx] + 5
    if (type == "scale") X[, , idx] <- X[, , idx] * 3
    if (type == "struct") X[, , idx] <- rmatnorm(length(idx), Sr = ar_cov(p, -0.5), Sc = ar_cov(q, -0.7)) }
  list(X = X, bad_case = bad_case, bad_cell = bad_cell)
}
one_rep <- function(seed, p, q, n, type, eps_case, eps_cell) {
  set.seed(seed); S <- true_sigma(p, q, "full"); d <- contaminate(rmatnorm(n, Sr = S$Sr, Sc = S$Sc), S, type, eps_case, eps_cell)
  X <- d$X; ev <- function(f) c(kl = kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), skl = kl_shape(f$Sr, f$Sc, S$Sr, S$Sc), tr = tr_ratio(f$Sr, f$Sc, S$Sr, S$Sc))
  out <- list()
  for (variant in c("soft", "hard_rw")) {
    f <- switch(variant, soft = cellmmcd(X), hard_rw = cellmmcd(X, quant_obs = NULL))
    out[[variant]] <- c(ev(f), trim_bad = if (any(d$bad_case)) mean(!f$v[d$bad_case]) else NA, trim_good = mean(!f$v[!d$bad_case]),
                    cell_recall = if (any(d$bad_cell)) mean(!f$W[d$bad_cell]) else NA, Merr = sum(f$M^2))
  }
  m <- tryCatch(mmcd(X, alpha = 0.75), error = function(e) NULL)
  out$mmcd <- c(kl = if (is.null(m)) NA else kl_mn(m$cov_row, m$cov_col, S$Sr, S$Sc), skl = NA, tr = NA, trim_bad = NA, trim_good = NA, cell_recall = NA, Merr = NA)
  ml <- flipflop(X, center = TRUE); out$mle <- c(kl = kl_mn(ml$Sr, ml$Sc, S$Sr, S$Sc), skl = NA, tr = NA, trim_bad = NA, trim_good = NA, cell_recall = NA, Merr = NA)
  out
}
cases <- rbind(data.frame(type = "none", eps_case = 0, eps_cell = 0), data.frame(type = "none", eps_case = 0, eps_cell = 0.10),
               data.frame(type = c("shift", "scale", "struct"), eps_case = 0.2, eps_cell = 0),
               data.frame(type = "shift", eps_case = 0.1, eps_cell = 0.05), data.frame(type = "struct", eps_case = 0.2, eps_cell = 0.05))
grid <- merge(expand.grid(p = c(5, 8), n = c(100, 200)), cases); R <- 10; results <- list()
sink("output/cm_step20_default.log"); cat("cellMMCD step 20 (block 7): soft default vs reweighted hard quota, 10 reps, medians. KL / shape KL / trace ratio / trimmed bad / trimmed good / cell recall / |M|^2.\n")
for (g in seq_len(nrow(grid))) {
  p <- q <- grid$p[g]; n <- grid$n[g]
  reps <- mclapply(1:R, function(r) one_rep(20000 + 100 * g + r, p, q, n, grid$type[g], grid$eps_case[g], grid$eps_cell[g]), mc.cores = 4)
  cat(sprintf("\n== %dx%d n = %d, case %s %.2f, cells %.2f ==\n", p, q, n, grid$type[g], grid$eps_case[g], grid$eps_cell[g]))
  tab <- t(sapply(names(reps[[1]]), function(m) apply(do.call(rbind, lapply(reps, `[[`, m)), 2, median, na.rm = TRUE)))
  print(round(tab, 3)); results[[g]] <- list(grid = grid[g, ], reps = reps)
}
saveRDS(results, "output/cm_step20_default.rds"); sink(); cat("done\n")
