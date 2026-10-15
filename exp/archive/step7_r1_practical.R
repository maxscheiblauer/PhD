# Step 7 (R1 convergence, practical): how often does R1 converge, and how good is it vs. MLE and MMCD?
# Run from project root: Rscript exp/archive/step7_r1_practical.R      (~15-25 min on 15 cores)
# R1 = r1_fit (= tff R1 with kappa off, det scale, sym roots), 6 starts per data set, maxit 1000.
# Grid: (p, q) x n x beta x data type, 8 reps. Contamination eps = .1, delta = 5 (R/data.R);
# "t3": elliptical matrix-t with 3 df (X_i / sqrt(chi2_3 / 3)), same shape as the Gaussian.
# Accuracy: KL_shape (scale-free) of the product to the truth; KL with kappa(beta, q) on the product.
# MMCD: robustmatrix::mmcd(alpha = beta, raw, no reweighting); it estimates M, R1 and MLE know M = 0.
# Per data set, the start with the smallest Lr~ + Lc~ is "R1 best". For the ff start also:
#   rho      local contraction rate (r1_jacobian) at the fixed point if converged, and rho of the
#            untrimmed flip-flop at the MLE for comparison;
#   monotonicity of Lr~, Lc~, J~ = Lr~ + Lc~ and of the residual gap along the run;
#   own vs. cross: in the late iterations (step KL in [1e-9, 1e-4]) the median of
#            |cross pull| / |own gain| of each criterion (theory: grows like 1 / step size).
source("R/load.R")
library(parallel)

grid <- expand.grid(pq = c("2x2", "3x3", "5x5", "2x6", "6x2"), n = c(10, 30, 100),
                    beta = c(0.5, 0.75, 0.9),
                    type = c("clean", "t3", "obs", "col", "row", "cell", "wcol"),
                    rep = 1:8, stringsAsFactors = FALSE)
grid$seed <- 7000 + seq_len(nrow(grid))

one <- function(g) {
  set.seed(g$seed)
  p <- as.integer(substr(g$pq, 1, 1)); q <- as.integer(substr(g$pq, 3, 3)); n <- g$n; beta <- g$beta
  S <- true_sigma(p, q, "full")
  if (g$type == "t3") {
    X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
    X <- sweep(X, 3, sqrt(rchisq(n, 3) / 3), "/")
  } else X <- gen_data(n, S$Sr, S$Sc, g$type, if (g$type == "clean") 0 else 0.1)$X
  hr <- round(beta * n * q); hc <- round(beta * n * p)
  out <- data.frame(g[rep(1, 1), ], p = p, q = q)

  ff <- flipflop(X, scale_id = "det")
  starts <- c(list(identity = list(Sr = diag(p), Sc = diag(q)), ff = ff[c("Sr", "Sc")]),
              lapply(1:4, function(k) tryCatch(init_elemental(X), error = function(e) NULL)))
  starts <- starts[!sapply(starts, is.null)]
  fits <- lapply(seq_along(starts), function(s)
    tryCatch(r1_fit(X, beta, starts[[s]]$Sr, starts[[s]]$Sc, maxit = 1000, diag = s == 2),
             error = function(e) NULL))
  ok <- !sapply(fits, is.null)
  st <- sapply(fits, function(f) if (is.null(f)) "error" else f$status)
  Jt <- sapply(fits, function(f) if (is.null(f) || f$status == "singular") Inf else
    sum(r1_crit(X, f$Sr, f$Sc, f$hr, f$hc)))
  best <- which.min(Jt); fb <- fits[[best]]
  kq <- kappa_mcd(hc / (n * p), q)

  out$n_starts <- length(fits); out$n_conv <- sum(st == "converged"); out$n_cycle <- sum(st == "cycle")
  out$n_maxit <- sum(st == "maxit"); out$n_sing <- sum(st %in% c("singular", "error"))
  out$iter_med <- median(sapply(fits[ok], `[[`, "iter"))
  out$best_status <- fb$status; out$ff_status <- st[2]
  out$n_distinct <- length(unique(sapply(fits[ok & st == "converged"], function(f) paste(which(f$a), collapse = "."))))
  out$KLs_MLE <- kl_shape(ff$Sr, ff$Sc, S$Sr, S$Sc); out$KL_MLE <- kl_mn(ff$Sr, ff$Sc, S$Sr, S$Sc)
  out$KLs_R1 <- kl_shape(fb$Sr, fb$Sc, S$Sr, S$Sc); out$KL_R1 <- kl_mn(fb$Sr * kq, fb$Sc, S$Sr, S$Sc)
  out$KLs_R1ff <- if (ok[2]) kl_shape(fits[[2]]$Sr, fits[[2]]$Sc, S$Sr, S$Sc) else NA
  # cycling / non-converged best start: how far apart are the states it keeps visiting?
  out$amp_best <- if (fb$status == "converged") 0 else max(tail(fb$step, 50))

  mm <- tryCatch(robustmatrix::mmcd(X, alpha = beta, adapt_alpha = FALSE, reweight = FALSE),
                 error = function(e) NULL)
  out$KLs_MMCD <- if (is.null(mm)) NA else kl_shape(mm$cov_row, mm$cov_col, S$Sr, S$Sc)
  out$KL_MMCD <- if (is.null(mm)) NA else kl_mn(mm$cov_row, mm$cov_col, S$Sr, S$Sc)

  ## Q4 diagnostics on the ff-start run
  f2 <- fits[[2]]
  out$rho <- out$rho_ff <- out$Lr_inc <- out$Lc_inc <- out$J_inc <- out$gap_inc <- NA
  out$late_ratio_r <- out$late_ratio_c <- out$late_Jinc <- out$late_n <- NA
  if (!is.null(f2) && q > 1) {
    if (f2$status == "converged") {
      out$rho <- max(Mod(r1_jacobian(X, f2$Sc, f2$a, f2$b)))
      out$rho_ff <- max(Mod(r1_jacobian(X, ff$Sc, matrix(TRUE, n, q), matrix(TRUE, n, p))))
    }
    cr <- f2$crit; K <- f2$iter; tolr <- 1e-9 * max(abs(cr), 1)
    dr <- diff(cr[, "Lr"]); dc <- diff(cr[, "Lc"]); dJ <- diff(rowSums(cr))
    out$Lr_inc <- sum(dr > tolr); out$Lc_inc <- sum(dc > tolr); out$J_inc <- sum(dJ > tolr)
    out$gap_inc <- sum(diff(f2$gap) > tolr)
    # half-steps: odd = row step (own for Lr, cross for Lc), even = column step
    late <- which(f2$step >= 1e-9 & f2$step <= 1e-4)
    late <- late[late >= 2]
    if (length(late)) {
      own_r <- -dr[2 * late - 1]; cross_r <- dr[2 * late]
      own_c <- -dc[2 * late];     cross_c <- dc[2 * late - 1]
      out$late_ratio_r <- median(abs(cross_r) / pmax(own_r, 1e-300))
      out$late_ratio_c <- median(abs(cross_c) / pmax(own_c, 1e-300))
      out$late_Jinc <- mean(c(dJ[2 * late - 1], dJ[2 * late]) > 0)
      out$late_n <- length(late)
    }
  }
  out
}

t0 <- proc.time()[3]
cl <- makeCluster(15)
invisible(clusterEvalQ(cl, { source("R/load.R"); NULL }))
clusterExport(cl, "one")
res <- parLapplyLB(cl, split(grid, seq_len(nrow(grid))), function(g)
  tryCatch(one(g), error = function(e) conditionMessage(e)))
stopCluster(cl)
err <- !sapply(res, is.data.frame)
if (any(err)) print(table(unlist(res[err])))
res <- do.call(rbind, res[!err])
elapsed <- proc.time()[3] - t0
saveRDS(list(res = res, elapsed = elapsed), "output/step7.rds")
cat(sprintf("Step 7 done in %.0f s, %d data sets (%d failed)\n", elapsed, nrow(res), nrow(grid) - nrow(res)))
