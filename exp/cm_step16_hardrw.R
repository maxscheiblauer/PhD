# Block 6: the reweighted hard quota (quant_rw_obs, default for quant_obs = NULL) against the raw hard quota
# and the soft default, on the casewise / mixed grid of cm_step4_mixed (same seeds, same contamination code),
# p = q = 5, n = 100, 10 reps. Output: output/cm_step16_hardrw.log / .rds
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
  out <- list()
  for (variant in c("soft", "hard_raw", "hard_rw")) {
    f <- switch(variant, soft = cellmmcd(X), hard_raw = cellmmcd(X, quant_obs = NULL, quant_rw_obs = FALSE), hard_rw = cellmmcd(X, quant_obs = NULL))
    fl <- apply(!f$W, 3, mean)
    out[[variant]] <- c(ev(f), trim_bad = if (any(d$bad_case)) mean(!f$v[d$bad_case]) else NA, trim_good = mean(!f$v[!d$bad_case]),
                    cell_recall = if (any(d$bad_cell)) mean(!f$W[d$bad_cell]) else NA, flag_good = mean(fl[!d$bad_case & !apply(d$bad_cell, 3, any)]), Merr = sum(f$M^2))
  }
  out
}
grid <- rbind(data.frame(type = c("shift", "scale", "struct"), eps_case = 0.1, eps_cell = 0),
              data.frame(type = c("shift", "scale", "struct"), eps_case = 0.2, eps_cell = 0),
              data.frame(type = "shift", eps_case = 0.1, eps_cell = 0.05), data.frame(type = "none", eps_case = 0, eps_cell = 0.05), data.frame(type = "none", eps_case = 0, eps_cell = 0),
              data.frame(type = "struct", eps_case = 0.2, eps_cell = 0.05))
p <- q <- 5; n <- 100; R <- 10; results <- list()
sink("output/cm_step16_hardrw.log"); cat("cellMMCD step 16 (block 6): reweighted hard quota. p = q = 5, n = 100, 10 reps. Medians.\n")
for (g in seq_len(nrow(grid))) {
  reps <- mclapply(1:R, function(r) one_rep(4000 + 100 * g + r, p, q, n, grid$type[g], grid$eps_case[g], grid$eps_cell[g]), mc.cores = 4)
  cat(sprintf("\n== case %s %.2f, cells %.2f ==\n", grid$type[g], grid$eps_case[g], grid$eps_cell[g]))
  tab <- t(sapply(names(reps[[1]]), function(m) apply(do.call(rbind, lapply(reps, `[[`, m)), 2, median, na.rm = TRUE)))
  print(round(tab, 3)); results[[g]] <- list(grid = grid[g, ], reps = reps)
}
saveRDS(results, "output/cm_step16_hardrw.rds"); sink()
