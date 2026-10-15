# R1 study: slice-trimmed flip-flop for matrix-variate data. Base R only.
# Conventions: X is a p x q x n array; vec(X_i) ~ N(vec M, Sc %x% Sr), Sr = p x p (rows), Sc = q x q (columns).
# Functions never read globals: everything they need is an argument.
# Naming: local variables are snake_case (x = X, m = M, s_r = Sr, s_c = Sc, ...); list fields keep the paper's names.

library(robustmatrix)


## ---- 1. Data -------------------------------------------------------------------------------

# Clean matrix-normal sample with AR(1) row and column covariances.
#   X_i = M + Lr Z_i Lc',  Z_i iid N(0,1) entries,  Lr Lr' = Sr,  Lc Lc' = Sc
# Returns the data together with the true parameters (needed later for evaluation).
gen_matnorm <- function(n, p, q, rho_r, rho_c, m = matrix(0, p, q), seed = NULL) {
  stopifnot(n >= 1, p >= 1, q >= 1, abs(rho_r) < 1, abs(rho_c) < 1,
            is.matrix(m), all(dim(m) == c(p, q)))
  if (!is.null(seed)) set.seed(seed)
  ar1 <- function(k, rho) rho^abs(outer(seq_len(k), seq_len(k), "-"))
  s_r <- ar1(p, rho_r)
  s_c <- ar1(q, rho_c)
  l_r <- t(chol(s_r))
  l_c <- t(chol(s_c))
  x <- vapply(seq_len(n),                                   # p x q x n, also for p = 1 or q = 1
              function(i) m + l_r %*% matrix(rnorm(p * q), p, q) %*% t(l_c),
              matrix(0, p, q))
  list(X = x, M = m, Sr = s_r, Sc = s_c, n = n, p = p, q = q)
}

# Replace a fraction eps of the units of ONE type by draws with sd_factor times the standard deviation.
# The outliers keep the model's correlation structure, only their scale is inflated:
#   type = "row": eps of all n*p rows      row k of X_i  -> N_q(m_k, sd_factor^2 Sr[k,k] Sc)
#        = "col": eps of all n*q columns   col j of X_i  -> N_p(m_j, sd_factor^2 Sc[j,j] Sr)
#        = "mat": eps of all n matrices    X_i           -> MN(M, sd_factor^2 Sr, Sc)
# (A replaced row/column is drawn independently of the rest of its matrix.)
# Returns the clean list plus:
#   out  : the reference for this type - p x n (row), q x n (col) or length n (mat); TRUE = replaced
#   mask : p x q x n, TRUE = replaced cell (same information at cell level)
contaminate <- function(clean, type = c("row", "col", "mat"), eps = 0.1, sd_factor = 6, seed = NULL) {
  type <- match.arg(type)
  stopifnot(eps >= 0, eps <= 1, sd_factor > 0)
  if (!is.null(seed)) set.seed(seed)
  p <- clean$p
  q <- clean$q
  n <- clean$n
  m <- clean$M
  l_r <- t(chol(clean$Sr))
  l_c <- t(chol(clean$Sc))
  x <- clean$X
  mask <- array(FALSE, c(p, q, n))

  if (type == "row") {
    out <- matrix(FALSE, p, n)
    out[sample.int(p * n, round(eps * p * n))] <- TRUE
    idx <- which(out, arr.ind = TRUE)                       # columns: row k, observation i
    for (r in seq_len(nrow(idx))) {
      k <- idx[r, 1]
      i <- idx[r, 2]
      x[k, , i] <- m[k, ] + sd_factor * sqrt(clean$Sr[k, k]) * drop(l_c %*% rnorm(q))
      mask[k, , i] <- TRUE
    }
  } else if (type == "col") {
    out <- matrix(FALSE, q, n)
    out[sample.int(q * n, round(eps * q * n))] <- TRUE
    idx <- which(out, arr.ind = TRUE)                       # columns: column j, observation i
    for (r in seq_len(nrow(idx))) {
      j <- idx[r, 1]
      i <- idx[r, 2]
      x[, j, i] <- m[, j] + sd_factor * sqrt(clean$Sc[j, j]) * drop(l_r %*% rnorm(p))
      mask[, j, i] <- TRUE
    }
  } else {
    out <- logical(n)
    out[sample.int(n, round(eps * n))] <- TRUE
    for (i in which(out)) {
      x[, , i] <- m + sd_factor * l_r %*% matrix(rnorm(p * q), p, q) %*% t(l_c)
      mask[, , i] <- TRUE
    }
  }

  dat <- clean
  dat$X_clean <- clean$X
  dat$X       <- x
  dat$out     <- out
  dat$mask    <- mask
  dat$contam  <- list(type = type, eps = eps, sd_factor = sd_factor)
  dat
}

## ---- 1b. Look at the data: squared MDs per row, column and matrix --------------------------

# Squared Mahalanobis distances of the original rows, the original columns and the whole matrix,
# for given parameters (here: the truth; later: estimates). R_i = X_i - M, Or = Sr^-1, Oc = Sc^-1.
#   row k   : r_k' Oc r_k / Sr[k,k]   ~ chi2_q   (marginal distribution of row k is N_q(m_k, Sr[k,k] Sc))
#   column j: c_j' Or c_j / Sc[j,j]   ~ chi2_p   (marginal distribution of col j is N_p(m_j, Sc[j,j] Sr))
#   matrix  : tr(Oc R' Or R)          ~ chi2_pq
md_slices <- function(x, m, s_r, s_c) {
  o_r <- solve(s_r)
  o_c <- solve(s_c)
  n <- dim(x)[3]
  out <- list(row = matrix(0, nrow(s_r), n), col = matrix(0, nrow(s_c), n), mat = numeric(n))
  for (i in seq_len(n)) {
    r_mat <- x[, , i] - m
    out$row[, i] <- rowSums((r_mat %*% o_c) * r_mat) / diag(s_r)
    out$col[, i] <- colSums(r_mat * (o_r %*% r_mat)) / diag(s_c)
    out$mat[i]   <- sum(r_mat * (o_r %*% r_mat %*% o_c))
  }
  out
}

# Four panels, x = observation i, y = squared MD (log scale), computed with the TRUE parameters.
# Each data set is shown with the distance that matches its outlier type:
#   (1) clean: matrix MDs   (2) row outliers: row MDs   (3) column outliers: column MDs   (4) matrix outliers: matrix MDs
# Red = replaced row / column / matrix (the reference dat$out). Dashed line = chi2 0.975 quantile.
plot_md <- function(clean, row, col, mat) {
  stopifnot(row$contam$type == "row", col$contam$type == "col", mat$contam$type == "mat")
  p <- clean$p
  q <- clean$q
  n <- clean$n
  md <- function(d) md_slices(d$X, clean$M, clean$Sr, clean$Sc)

  panel <- function(d_mat, out, df, main) {                 # d_mat: k x n distances, out: k x n reference
    d_mat <- matrix(d_mat, ncol = n)
    out <- matrix(out, ncol = n)
    x <- col(d_mat)
    plot(x, d_mat, log = "y", type = "n", xlab = "observation i", ylab = "squared MD", main = main)
    points(x[!out], d_mat[!out], pch = 16, cex = 0.6, col = "grey55")
    points(x[out],  d_mat[out],  pch = 16, cex = 0.8, col = "firebrick")
    abline(h = qchisq(0.975, df), lty = 2)
  }
  lab <- function(d, unit) {
    sprintf("%s outliers (%.0f%% of %s): %s MDs, df = %d",
            d$contam$type, 100 * d$contam$eps, unit, d$contam$type,
            switch(d$contam$type, row = q, col = p, mat = p * q))
  }

  op <- par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1), oma = c(0, 0, 2, 0))
  on.exit(par(op))
  panel(md(clean)$mat, rep(FALSE, n), p * q, sprintf("clean data: matrix MDs, df = %d", p * q))
  legend("topright", c("clean", "replaced (outlier)", "chi2 0.975"), pch = c(16, 16, NA),
         lty = c(NA, NA, 2), col = c("grey55", "firebrick", "black"), bty = "n", cex = 0.8)
  panel(md(row)$row, row$out, q,     lab(row, "rows"))
  panel(md(col)$col, col$out, p,     lab(col, "columns"))
  panel(md(mat)$mat, mat$out, p * q, lab(mat, "matrices"))
  mtext(sprintf("n=%d p=%d q=%d | outliers have %g x the sd | distances with true parameters",
                n, p, q, row$contam$sd_factor), outer = TRUE, cex = 0.85)
  invisible(NULL)
}

## ---- Run: data, outliers, plot ------------------------------------------------------------
clean   <- gen_matnorm(n = 100, p = 3, q = 4, rho_r = 0.6, rho_c = 0.5, seed = 1)
dat_row <- contaminate(clean, type = "row", eps = 0.1, sd_factor = 6, seed = 2)
dat_col <- contaminate(clean, type = "col", eps = 0.1, sd_factor = 6, seed = 3)
dat_mat <- contaminate(clean, type = "mat", eps = 0.1, sd_factor = 6, seed = 4)
dir.create("output", showWarnings = FALSE)
png(file.path("output", "md_types.png"), width = 1300, height = 1000, res = 110)
plot_md(clean, dat_row, dat_col, dat_mat)
dev.off()


####################################################################################################
#####################################################################################################



# Both return list(M, Sr, Sc, md2): md2 = squared matrix MD of every observation, from our md_slices().
as_fit <- function(f, x) {
  list(M = f$mu, Sr = f$cov_row, Sc = f$cov_col,
       md2 = md_slices(x, f$mu, f$cov_row, f$cov_col)$mat)
}
fit_mle  <- function(x) as_fit(robustmatrix::mmle(x, silent = TRUE), x)
fit_mmcd <- function(x, alpha = 0.75) as_fit(robustmatrix::mmcd(x, alpha = alpha), x)

# Divergence of (Sr, Sc) from the truth, eq. (31) of Mayrhofer et al. (2025); 0 = perfect.
kl_div <- function(s_r, s_c, s_r0, s_c0) {
  m_r <- solve(s_r0, s_r)
  m_c <- solve(s_c0, s_c)
  sum(diag(m_r)) * sum(diag(m_c)) - ncol(s_c) * determinant(m_r)$modulus -
    ncol(s_r) * determinant(m_c)$modulus - ncol(s_r) * ncol(s_c)
}

## ---- Run: MLE and MMCD on the four data sets -----------------------------------------------
sets <- list(clean = clean, row = dat_row, col = dat_col, mat = dat_mat)
res <- t(sapply(sets, function(d) {
  f1 <- fit_mle(d$X)
  f2 <- fit_mmcd(d$X)
  c(KL_MLE  = kl_div(f1$Sr, f1$Sc, clean$Sr, clean$Sc),
    KL_MMCD = kl_div(f2$Sr, f2$Sc, clean$Sr, clean$Sc),
    flag_MLE  = sum(f1$md2 > qchisq(0.975, clean$p * clean$q)),
    flag_MMCD = sum(f2$md2 > qchisq(0.975, clean$p * clean$q)))
}))
print(round(res, 3))



#########################################################################################################
##############################################################################################################
## ---- 2. Trimmed MLE (slice-trimmed flip-flop) ---------------------------------------------

## Parameters
beta  <- 0.75     # fraction of slices kept in each step
maxit <- 50       # max number of flip-flop sweeps

## Helpers
# Symmetric inverse square root S^(-1/2).
isqrt <- function(s) {
  e <- eigen(s, symmetric = TRUE)
  e$vectors %*% (t(e$vectors) / sqrt(e$values))
}

# Indices of the h smallest distances.
keep_h <- function(d, h) order(d)[seq_len(h)]

## Algorithm
# X: p x q x n. M: centre (default: cellwise median, held fixed).
# Converged when the kept sets of both steps no longer change.
# Trimmed -2 log-likelihood per observation (up to constants), recorded after every half step:
#   nll = q log|Sr| + p log|Sc| + q * mean(d_ij kept)   after the row step
#       = q log|Sr| + p log|Sc| + p * mean(e_ik kept)   after the column step
#   (without trimming both equal the matrix-normal -2 log-lik / n; invariant to Sr*c, Sc/c)
# Consistency scaling and reweighting, named and built as in robustmatrix::mmcd() (Mayrhofer et al. 2025):
#   scale_consistency = "quant"  : scale Sc %x% Sr by the consistency factor under the matrix normal distribution
#                                  raw:        c = beta / F_{chi2_{q+2}}(chi2_{beta; q})
#                                  (MMCD uses alpha / F_{chi2_{pq+2}}(chi2_{alpha; pq}) because it trims whole
#                                   matrices; here the raw size is set by the LAST half step, which trims
#                                   q-dimensional units with fraction beta, hence df = q)
#                                  reweighted: c = alpha_rw / F_{chi2_{pq+2}}(chi2_{alpha_rw; pq}), alpha_rw = h_rw / n
#                                  (as in MMCD; the reweighting step works on whole matrices, so df = pq)
#                     = "mmd_med": c = median(MMD^2) / chi2_{pq; 0.5}  (from the squared matrix MDs)
#                     = "none"   : no scaling (not available in mmcd; switch-off flag)
#   reweight = TRUE   : weights w_i = 1 if MMD^2(X_i) < chi2_{pq; outlier_quant}, else 0, computed with the scaled
#                       raw estimates; then the matrix-normal MLE on the matrices with w_i = 1 (M, Sr, Sc are all
#                       re-estimated, as in MMCD), followed by the reweighted consistency factor.
#                       (MMCD also sets w_i = 1 for i in the final h-subset; there is no matrix h-subset here.)
#   outlier_quant     : chi-square quantile used in the reweighting step.
# The consistency factor multiplies Sc; only the product Sc %x% Sr matters.
# Returns list(M, Sr, Sc, md2, md2_raw, consistency_factors, h_subset_reweighted, iter, converged, path, switches, nll):
#   md2_raw: squared matrix MDs with the raw (scaled) estimates;  consistency_factors: c(raw, reweighted)
#   h_subset_reweighted: matrices with weight one in the reweighting step (NULL if reweight = FALSE)
#   path: (Sr, Sc) after every sweep (raw, unscaled);  switches: slices that entered the kept sets in that sweep
fit_tmle <- function(x, beta = 0.75, maxit = 50, m = apply(x, 1:2, median),
                     reweight = TRUE, scale_consistency = c("quant", "mmd_med", "none"),
                     outlier_quant = 0.975) {
  scale_consistency <- match.arg(scale_consistency)
  p <- dim(x)[1]
  q <- dim(x)[2]
  n <- dim(x)[3]
  hr <- floor(beta * n * q)
  hc <- floor(beta * n * p)
  x_c <- sweep(x, 1:2, m)
  s_r <- diag(p)
  s_c <- diag(q)
  a <- b <- integer(0)
  converged <- FALSE
  path <- list()
  switches <- integer(0)
  nll <- numeric(0)
  ld <- function(s) as.numeric(determinant(s)$modulus)

  for (it in seq_len(maxit)) {
    a_old <- a
    b_old <- b

    # 1. Row step: z_ij = column j of X_i Sc^(-1/2)  ->  Z is p x (q n)
    w_c <- isqrt(s_c)
    z  <- matrix(apply(x_c, 3, function(x_i) x_i %*% w_c), p)
    a  <- sort(keep_h(colSums(z * solve(s_r, z)), hr))
    s_r <- tcrossprod(z[, a, drop = FALSE]) / hr
    z_a <- z[, a, drop = FALSE]
    nll <- c(nll, q * ld(s_r) + p * ld(s_c) + q * mean(colSums(z_a * solve(s_r, z_a))))

    # 2. Column step: u_ik = row k of Sr^(-1/2) X_i  ->  U is q x (p n)
    w_r <- isqrt(s_r)
    u  <- matrix(apply(x_c, 3, function(x_i) t(w_r %*% x_i)), q)
    b  <- sort(keep_h(colSums(u * solve(s_c, u)), hc))
    s_c <- tcrossprod(u[, b, drop = FALSE]) / hc
    u_b <- u[, b, drop = FALSE]
    nll <- c(nll, q * ld(s_r) + p * ld(s_c) + p * mean(colSums(u_b * solve(s_c, u_b))))

    path[[it]] <- list(Sr = s_r, Sc = s_c)                       # trace for plotting
    switches[it] <- length(setdiff(a, a_old)) + length(setdiff(b, b_old))
    if (identical(a, a_old) && identical(b, b_old)) {
      converged <- TRUE
      break
    }
  }

  # Consistency factor: c_quant(a, df) = a / F_{chi2_{df+2}}(chi2_{a; df})
  pq <- p * q
  c_quant <- function(a, df) a / pchisq(qchisq(a, df), df + 2)
  c_fac   <- function(a, df, md2) {
    switch(scale_consistency,
           quant   = c_quant(a, df),
           mmd_med = median(md2) / qchisq(0.5, pq),
           none    = 1)
  }

  # Raw estimators, scaled
  c_raw   <- c_fac(beta, q, md_slices(x, m, s_r, s_c)$mat)
  s_c     <- c_raw * s_c
  md2_raw <- md_slices(x, m, s_r, s_c)$mat
  c_rw    <- NA_real_
  h_rw    <- NULL

  # Reweighted estimators
  if (reweight) {
    h_rw <- which(md2_raw < qchisq(outlier_quant, pq))
    g    <- robustmatrix::mmle(x[, , h_rw, drop = FALSE], silent = TRUE)
    m    <- g$mu
    s_r  <- g$cov_row
    s_c  <- g$cov_col
    c_rw <- c_fac(length(h_rw) / n, pq, md_slices(x, m, s_r, s_c)$mat)
    s_c  <- c_rw * s_c
  }

  list(M = m, Sr = s_r, Sc = s_c, md2 = md_slices(x, m, s_r, s_c)$mat, md2_raw = md2_raw,
       consistency_factors = c(raw = c_raw, reweighted = c_rw), h_subset_reweighted = h_rw,
       iter = it, converged = converged, path = path, switches = switches, nll = nll)
}

## ---- 2b. Plot one trimmed-MLE fit -----------------------------------------------------------

# Four panels for one data set d (clean or contaminated) and its fit f = fit_tmle(d$X):
#   (1) KL to the truth after every full sweep (numbers = slices that switched in)
#   (2) trimmed -2 log-likelihood after the row step and after the column step of every sweep
#   (3) squared MDs with the estimates after convergence, (4) the same with the true parameters.
#       The distance matches the outlier type: row MDs (row data), column MDs (col data),
#       matrix MDs (clean and mat data). Red = true outliers, dashed = chi2 0.975 quantile.
plot_tmle <- function(d, f, s_r0, s_c0) {
  p <- nrow(s_r0)
  q <- nrow(s_c0)
  n <- dim(d$X)[3]
  type <- if (is.null(d$contam)) "clean" else d$contam$type
  slot <- if (type == "clean") "mat" else type
  df   <- switch(slot, row = q, col = p, mat = p * q)
  out  <- if (type == "clean") rep(FALSE, n) else d$out

  d_est  <- matrix(md_slices(d$X, f$M, f$Sr, f$Sc)[[slot]], ncol = n)
  d_true <- matrix(md_slices(d$X, d$M, s_r0, s_c0)[[slot]], ncol = n)
  ylim   <- range(d_est, d_true)

  md_panel <- function(d_mat, main) {
    x <- col(d_mat)
    o <- matrix(out, ncol = n)
    plot(x, d_mat, log = "y", type = "n", ylim = ylim, xlab = "observation i",
         ylab = sprintf("squared %s MD", slot), main = main)
    points(x[!o], d_mat[!o], pch = 16, cex = 0.6, col = "grey55")
    points(x[o],  d_mat[o],  pch = 16, cex = 0.8, col = "firebrick")
    abline(h = qchisq(0.975, df), lty = 2)
  }

  op <- par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1), oma = c(0, 0, 2, 0))
  on.exit(par(op))

  # (1) KL to the truth per sweep
  kl <- sapply(f$path, function(s) kl_div(s$Sr, s$Sc, s_r0, s_c0))
  plot(seq_along(kl), kl, type = "b", pch = 16, log = "y", xlab = "full sweep", ylab = "KL to truth",
       main = sprintf("KL to truth: %d sweeps%s", f$iter, if (f$converged) "" else " (not converged)"))
  text(seq_along(kl), kl, f$switches, pos = 3, cex = 0.7, col = "grey40", xpd = NA)

  # (2) trimmed -2 log-likelihood per half step
  r <- f$nll[c(TRUE, FALSE)]
  cc <- f$nll[c(FALSE, TRUE)]
  matplot(seq_along(r), cbind(r, cc), type = "b", pch = c(1, 4), lty = 1, col = c("black", "steelblue"),
          xlab = "full sweep", ylab = "-2 trimmed log-lik / n", main = "trimmed -2 log-likelihood")
  legend("topright", c("after row step", "after column step"), pch = c(1, 4), lty = 1,
         col = c("black", "steelblue"), bty = "n", cex = 0.8)

  # (3)-(4) MDs, estimated vs true parameters
  md_panel(d_est,  sprintf("%s MDs, estimated Sr, Sc (df = %d)", slot, df))
  md_panel(d_true, sprintf("%s MDs, true Sr, Sc (df = %d)", slot, df))

  mtext(sprintf("trimmed MLE, %s data | n=%d p=%d q=%d | red = true outliers", type, n, p, q),
        outer = TRUE, cex = 0.9)
  invisible(NULL)
}

## ---- Run: MLE, MMCD and trimmed MLE on the four data sets ----------------------------------
res <- t(sapply(sets, function(d) {
  f1 <- fit_mle(d$X)
  f2 <- fit_mmcd(d$X)
  f3 <- fit_tmle(d$X, beta, maxit)
  cut <- qchisq(0.975, clean$p * clean$q)
  c(KL_MLE    = kl_div(f1$Sr, f1$Sc, clean$Sr, clean$Sc),
    KL_MMCD   = kl_div(f2$Sr, f2$Sc, clean$Sr, clean$Sc),
    KL_TMLE   = kl_div(f3$Sr, f3$Sc, clean$Sr, clean$Sc),
    flag_MLE  = sum(f1$md2 > cut),
    flag_MMCD = sum(f2$md2 > cut),
    flag_TMLE = sum(f3$md2 > cut),
    iter_TMLE = f3$iter,
    conv_TMLE = f3$converged)
}))
print(round(res, 3))

## ---- Run: one figure per data set -----------------------------------------------------------
for (nm in names(sets)) {
  f <- fit_tmle(sets[[nm]]$X, beta, maxit)
  png(file.path("output", sprintf("tmle_%s.png", nm)), width = 1300, height = 1000, res = 110)
  plot_tmle(sets[[nm]], f, clean$Sr, clean$Sc)
  dev.off()
}
