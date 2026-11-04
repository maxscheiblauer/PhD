# cellMMCD step 4: casewise and mixed contamination (brief: "ideally casewise too"). Seed 4.
# Casewise types on 10 % of the matrices: shift (all cells + 5), scale (X * 3), struct (drawn from
# the Kronecker model with the correlation signs flipped: not marginally outlying, only relationally).
# Mixed: 5 % cells at 6 sd plus 10 % shifted matrices. Methods: FF (MLE), cellMMCD (calibrated), MMCD.
# Output: output/cm_step4_mixed.log, output/cm_step4_mixed.rds
source("R/load.R"); library(parallel); set.seed(4)
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
  out <- list(FF = ev(flipflop(X)))
  f <- cellmmcd(X); fl <- apply(!f$W, 3, mean)        # fraction of flagged cells per matrix
  out$cellMMCD <- c(ev(f), flag_bad = if (any(d$bad_case)) mean(fl[d$bad_case]) else NA, flag_good = mean(fl[!d$bad_case & !apply(d$bad_cell, 3, any)]),
                    cell_recall = if (any(d$bad_cell)) mean(!f$W[d$bad_cell]) else NA, iter = f$iter)
  m <- tryCatch(robustmatrix::mmcd(X, alpha = 0.75), error = function(e) NULL)
  out$MMCD <- if (is.null(m)) c(kl = NA, skl = NA, tr = NA) else ev(list(Sr = m$cov_row, Sc = m$cov_col))
  out
}
grid <- rbind(data.frame(type = c("shift", "scale", "struct"), eps_case = 0.1, eps_cell = 0),
              data.frame(type = c("shift", "scale", "struct"), eps_case = 0.2, eps_cell = 0),
              data.frame(type = "shift", eps_case = 0.1, eps_cell = 0.05), data.frame(type = "none", eps_case = 0, eps_cell = 0.05))
p <- q <- 5; n <- 100; R <- 10; results <- list()
sink("output/cm_step4_mixed.log"); cat("cellMMCD step 4: p = q = 5, n = 100, 10 reps. Medians.\n")
for (g in seq_len(nrow(grid))) {
  reps <- mclapply(1:R, function(r) one_rep(4000 + 100 * g + r, p, q, n, grid$type[g], grid$eps_case[g], grid$eps_cell[g]), mc.cores = 4)
  cm <- do.call(rbind, lapply(reps, `[[`, "cellMMCD"))
  cat(sprintf("\n== case %s %.2f, cells %.2f ==  cellMMCD flagged: %.2f of cells in bad matrices, %.3f in clean matrices; cell recall %.2f; iter %.0f\n",
              grid$type[g], grid$eps_case[g], grid$eps_cell[g], median(cm[, "flag_bad"]), median(cm[, "flag_good"]), median(cm[, "cell_recall"]), median(cm[, "iter"])))
  tab <- t(sapply(c("FF", "cellMMCD", "MMCD"), function(m) { v <- do.call(rbind, lapply(reps, `[[`, m)); apply(v[, c("kl", "skl", "tr"), drop = FALSE], 2, median, na.rm = TRUE) }))
  print(round(tab, 3)); results[[g]] <- list(grid = grid[g, ], reps = reps)
}
sink(); saveRDS(results, "output/cm_step4_mixed.rds"); cat("done\n")
