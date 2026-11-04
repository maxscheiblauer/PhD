# R1 study, cellwise outliers only: slice-trimmed flip-flop for matrix-variate data. Base R only.
# Almost a copy of r1_playground.R; the row / column / matrix outliers are replaced by cellwise outliers.
# Conventions: X is a p x q x n array; vec(X_i) ~ N(vec M, Sc %x% Sr), Sr = p x p (rows), Sc = q x q (columns).
# Functions never read globals: everything they need is an argument.
# Naming: local variables are snake_case (x = X, m = M, s_r = Sr, s_c = Sc, ...); list fields keep the paper's names.

library(robustmatrix)


## ---- 1. Data -------------------------------------------------------------------------------

# Clean matrix-normal sample with AR(1) row and column covariances.
#   X_i = M + Lr Z_i Lc',  Z_i iid N(0,1) entries,  Lr Lr' = Sr,  Lc Lc' = Sc
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

# Replace a fraction eps of ALL n*p*q cells independently (cellwise contamination).
# The replaced cell is drawn independently of the rest of its matrix; sd_s = sqrt(Sr[k,k] Sc[j,j]) is its model sd:
#   type = "scale": x_kji -> m_kj + sd_factor * sd_s * N(0,1)            (inflated variance)
#        = "shift": x_kji -> m_kj + sd_factor * sd_s + sd_s * N(0,1)     (mean shifted by sd_factor sds)
# Returns the clean list plus:
#   mask : p x q x n, TRUE = replaced cell
#   out  : length n,  TRUE = matrix contains at least one replaced cell
#   rows : p x n,     TRUE = row of X_i contains at least one replaced cell
#   cols : q x n,     TRUE = column of X_i contains at least one replaced cell
contaminate_cell <- function(clean, eps = 0.05, sd_factor = 6, type = c("scale", "shift"), seed = NULL) {
  type <- match.arg(type)
  stopifnot(eps >= 0, eps <= 1, sd_factor > 0)
  if (!is.null(seed)) set.seed(seed)
  p <- clean$p
  q <- clean$q
  n <- clean$n
  sd_s <- sqrt(outer(diag(clean$Sr), diag(clean$Sc)))       # p x q
  mask <- array(FALSE, c(p, q, n))
  mask[sample.int(p * q * n, round(eps * p * q * n))] <- TRUE
  x <- clean$X
  for (i in which(apply(mask, 3, any))) {
    z  <- matrix(rnorm(p * q), p, q)
    xi <- x[, , i]
    mk <- mask[, , i]
    xi[mk] <- switch(type,
                     scale = (clean$M + sd_factor * sd_s * z)[mk],
                     shift = (clean$M + sd_factor * sd_s + sd_s * z)[mk])
    x[, , i] <- xi
  }

  dat <- clean
  dat$X_clean <- clean$X
  dat$X       <- x
  dat$mask    <- mask
  dat$out     <- apply(mask, 3, any)
  dat$rows    <- apply(mask, c(1, 3), any)
  dat$cols    <- apply(mask, c(2, 3), any)
  dat$contam  <- list(type = "cell", cell_type = type, eps = eps, sd_factor = sd_factor)
  dat
}

## ---- 1b. Look at the data: squared MDs per row, column and matrix --------------------------

# Squared Mahalanobis distances of the original rows, the original columns and the whole matrix,
# for given parameters (here: the truth; later: estimates). R_i = X_i - M, Or = Sr^-1, Oc = Sc^-1.
#   row k   : r_k' Oc r_k / Sr[k,k]   ~ chi2_q
#   column j: c_j' Or c_j / Sc[j,j]   ~ chi2_p
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

# Three panels for one cellwise-contaminated data set, distances with the TRUE parameters:
#   row MDs, column MDs, matrix MDs. Red = row / column / matrix that contains at least one replaced cell.
plot_md_cell <- function(d) {
  p <- d$p
  q <- d$q
  n <- d$n
  md <- md_slices(d$X, d$M, d$Sr, d$Sc)
  panel <- function(d_mat, out, df, main) {
    d_mat <- matrix(d_mat, ncol = n)
    out <- matrix(out, ncol = n)
    x <- col(d_mat)
    plot(x, d_mat, log = "y", type = "n", xlab = "observation i", ylab = "squared MD", main = main)
    points(x[!out], d_mat[!out], pch = 16, cex = 0.6, col = "grey55")
    points(x[out],  d_mat[out],  pch = 16, cex = 0.8, col = "firebrick")
    abline(h = qchisq(0.975, df), lty = 2)
  }
  op <- par(mfrow = c(1, 3), mar = c(4, 4.5, 3, 1), oma = c(0, 0, 2, 0))
  on.exit(par(op))
  panel(md$row, d$rows, q,     sprintf("row MDs, df = %d", q))
  panel(md$col, d$cols, p,     sprintf("column MDs, df = %d", p))
  panel(md$mat, d$out,  p * q, sprintf("matrix MDs, df = %d", p * q))
  legend("topright", c("no replaced cell", "contains replaced cell", "chi2 0.975"), pch = c(16, 16, NA),
         lty = c(NA, NA, 2), col = c("grey55", "firebrick", "black"), bty = "n", cex = 0.8)
  mtext(sprintf("cellwise outliers (%s, %.0f%% of cells, %g sd) | n=%d p=%d q=%d | true parameters",
                d$contam$cell_type, 100 * d$contam$eps, d$contam$sd_factor, n, p, q),
        outer = TRUE, cex = 0.85)
  invisible(NULL)
}

## ---- Run: data, outliers, plot ------------------------------------------------------------
clean    <- gen_matnorm(n = 100, p = 3, q = 4, rho_r = 0.6, rho_c = 0.5, seed = 1)
dat_cell <- contaminate_cell(clean, eps = 0.05, sd_factor = 6, type = "scale", seed = 2)
dat_shft <- contaminate_cell(clean, eps = 0.05, sd_factor = 6, type = "shift", seed = 3)
dir.create("output", showWarnings = FALSE)
png(file.path("output", "md_cellwise.png"), width = 1500, height = 560, res = 110)
plot_md_cell(dat_cell)
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


#########################################################################################################
##############################################################################################################
## ---- 2. Trimmed MLE (slice-trimmed flip-flop) ---------------------------------------------
# Unchanged copy of the algorithm in r1_playground.R.

## Parameters
beta  <- 0.5      # fraction of slices kept in each step
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
# Consistency scaling and reweighting, named and built as in robustmatrix::mmcd() (Mayrhofer et al. 2025):
#   scale_consistency = "quant"  : raw: c = beta / F_{chi2_{q+2}}(chi2_{beta; q})
#                                  reweighted: c = alpha_rw / F_{chi2_{pq+2}}(chi2_{alpha_rw; pq}), alpha_rw = h_rw / n
#                     = "mmd_med": c = median(MMD^2) / chi2_{pq; 0.5}
#                     = "none"   : no scaling
#   reweight = TRUE   : w_i = 1 if MMD^2(X_i) < chi2_{pq; outlier_quant}; matrix-normal MLE on those matrices.
# The consistency factor multiplies Sc; only the product Sc %x% Sr matters.
# Returns list(M, Sr, Sc, md2, md2_raw, consistency_factors, h_subset_reweighted, iter, converged, path, switches, nll)
fit_tmle <- function(x, beta = 0.5, maxit = 50, m = apply(x, 1:2, median),
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

## ---- 2a. Row / column reweighting (ECM) ---------------------------------------------------------
# Alternative to the matrix reweighting in fit_tmle(): instead of dropping whole matrices, drop the
# flagged rows OR the flagged columns of each matrix (never both) and fit the matrix normal to the rest.

# Flags from the marginal row / column MDs of md_slices() (row k ~ chi2_q, column j ~ chi2_p).
# A row MD uses only that row, so one bad cell does not spill into the other rows (same for columns).
# A cell outlier flags its row and its column; each matrix drops the cheaper side (#rows * q vs #cols * p).
# Returns drop_r (p x n) and drop_c (q x n); in every matrix at most one of the two is non-empty.
flag_rowcol <- function(x, m, s_r, s_c, quant = 0.99) {
  p <- dim(x)[1]
  q <- dim(x)[2]
  md <- md_slices(x, m, s_r, s_c)
  fr <- md$row > qchisq(quant, q)
  fc <- md$col > qchisq(quant, p)
  use_r <- colSums(fr) * q <= colSums(fc) * p
  fr[, !use_r] <- FALSE
  fc[, use_r]  <- FALSE
  list(drop_r = fr, drop_c = fc)
}

# ECM for the matrix normal when whole rows OR whole columns of X_i are missing.
# Rows m missing (o = observed), Or = Sr^-1:
#   X_m | X_o ~ MN(M_m - C_r Or[m,o] (X_o - M_o), C_r, Sc),  C_r = (Or[m,m])^-1   (partitioned inverse)
#   E[R Oc R'] = Rh Oc Rh' + tr(Oc Sc_old) pad(C_r),   E[R' Or R] = Rh' Or Rh + tr(Or[m,m] C_r) Sc_old
# Columns missing: the same with the roles swapped. Per matrix only an |m| x |m| inverse, never pq x pq,
# so one iteration costs O(n (p^2 q + p q^2) + p^3 + q^3), the order of one flip-flop sweep.
# One iteration = E-step, M, Sr given the old Sc, Sc given the new Sr. Returns list(M, Sr, Sc, iter).
ecm_rowcol <- function(x, drop_r, drop_c, m, s_r, s_c, maxit = 200, tol = 1e-8) {
  p <- dim(x)[1]
  q <- dim(x)[2]
  n <- dim(x)[3]
  ir <- which(colSums(drop_r) > 0)
  ic <- which(colSums(drop_c) > 0)
  for (it in seq_len(maxit)) {
    o_r <- solve(s_r)
    o_c <- solve(s_c)
    xh <- x
    c_r <- c_c <- vector("list", n)

    # E-step: conditional means of the dropped rows / columns, conditional row / column covariances
    for (i in ir) {
      mi <- drop_r[, i]
      c_r[[i]] <- solve(o_r[mi, mi, drop = FALSE])
      xh[mi, , i] <- m[mi, , drop = FALSE] -
        c_r[[i]] %*% o_r[mi, !mi, drop = FALSE] %*% (x[!mi, , i] - m[!mi, , drop = FALSE])
    }
    for (i in ic) {
      mi <- drop_c[, i]
      c_c[[i]] <- solve(o_c[mi, mi, drop = FALSE])
      xh[, mi, i] <- m[, mi, drop = FALSE] -
        (x[, !mi, i] - m[, !mi, drop = FALSE]) %*% o_c[!mi, mi, drop = FALSE] %*% c_c[[i]]
    }
    m_new <- apply(xh, 1:2, mean)
    r <- sweep(xh, 1:2, m_new)

    # Sr given Sc_old (tr(Oc Sc_old) = q)
    w_c <- isqrt(s_c)
    z <- matrix(apply(r, 3, function(r_i) r_i %*% w_c), p)
    a_r <- tcrossprod(z)
    for (i in ir) {
      mi <- drop_r[, i]
      a_r[mi, mi] <- a_r[mi, mi] + q * c_r[[i]]
    }
    for (i in ic) a_r <- a_r + sum(o_c[drop_c[, i], drop_c[, i]] * c_c[[i]]) * s_r
    s_r_new <- a_r / (n * q)

    # Sc given Sr_new
    o_rn <- solve(s_r_new)
    w_r <- isqrt(s_r_new)
    u <- matrix(apply(r, 3, function(r_i) t(w_r %*% r_i)), q)
    a_c <- tcrossprod(u)
    for (i in ir) a_c <- a_c + sum(o_rn[drop_r[, i], drop_r[, i]] * c_r[[i]]) * s_c
    for (i in ic) {
      mi <- drop_c[, i]
      a_c[mi, mi] <- a_c[mi, mi] + sum(o_rn * s_r) * c_c[[i]]
    }
    s_c_new <- a_c / (n * p)

    done <- max(abs(s_r_new - s_r)) + max(abs(s_c_new - s_c)) + max(abs(m_new - m)) < tol
    m <- m_new
    s_r <- s_r_new
    s_c <- s_c_new
    if (done) break
  }
  list(M = m, Sr = s_r, Sc = s_c, iter = it)
}

# Raw R1 (scaled, no matrix reweighting) -> flag rows / columns -> ECM; n_rw rounds of re-flagging with the new fit.
# Returns the raw fit (path, nll, ... for plot_tmle_cell) with M, Sr, Sc, md2 replaced, plus
#   frac_dropped: share of all cells dropped in the last round;  iter_ecm: ECM iterations of the last round.
fit_tmle_rowcol <- function(x, beta = 0.5, maxit = 50, quant = 0.99, n_rw = 2) {
  f <- fit_tmle(x, beta, maxit, reweight = FALSE)
  g <- f
  for (k in seq_len(n_rw)) {
    fl <- flag_rowcol(x, g$M, g$Sr, g$Sc, quant)
    g  <- ecm_rowcol(x, fl$drop_r, fl$drop_c, g$M, g$Sr, g$Sc)
  }
  f$M  <- g$M
  f$Sr <- g$Sr
  f$Sc <- g$Sc
  f$md2 <- md_slices(x, g$M, g$Sr, g$Sc)$mat
  f$frac_dropped <- (sum(fl$drop_r) * dim(x)[2] + sum(fl$drop_c) * dim(x)[1]) / length(x)
  f$iter_ecm <- g$iter
  f
}

## ---- 2b. Plot one trimmed-MLE fit on cellwise data -------------------------------------------

# Four panels for one cellwise data set d and its fit f = fit_tmle(d$X):
#   (1) KL to the truth after every full sweep (numbers = slices that switched in)
#   (2) trimmed -2 log-likelihood after the row step and after the column step of every sweep
#   (3) squared matrix MDs with the estimates, (4) the same with the true parameters.
#       Red = matrices containing at least one replaced cell, dashed = chi2 0.975 quantile.
plot_tmle_cell <- function(d, f, s_r0, s_c0, label = "trimmed MLE") {
  p <- nrow(s_r0)
  q <- nrow(s_c0)
  n <- dim(d$X)[3]
  df <- p * q
  d_est  <- md_slices(d$X, f$M, f$Sr, f$Sc)$mat
  d_true <- md_slices(d$X, d$M, s_r0, s_c0)$mat
  ylim   <- range(d_est, d_true)

  md_panel <- function(dist, main) {
    plot(seq_len(n), dist, log = "y", type = "n", ylim = ylim, xlab = "observation i",
         ylab = "squared matrix MD", main = main)
    points(which(!d$out), dist[!d$out], pch = 16, cex = 0.6, col = "grey55")
    points(which(d$out),  dist[d$out],  pch = 16, cex = 0.8, col = "firebrick")
    abline(h = qchisq(0.975, df), lty = 2)
  }

  op <- par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1), oma = c(0, 0, 2, 0))
  on.exit(par(op))

  kl <- sapply(f$path, function(s) kl_div(s$Sr, s$Sc, s_r0, s_c0))
  plot(seq_along(kl), kl, type = "b", pch = 16, log = "y", xlab = "full sweep", ylab = "KL to truth",
       main = sprintf("KL to truth: %d sweeps%s", f$iter, if (f$converged) "" else " (not converged)"))
  text(seq_along(kl), kl, f$switches, pos = 3, cex = 0.7, col = "grey40", xpd = NA)

  r <- f$nll[c(TRUE, FALSE)]
  cc <- f$nll[c(FALSE, TRUE)]
  matplot(seq_along(r), cbind(r, cc), type = "b", pch = c(1, 4), lty = 1, col = c("black", "steelblue"),
          xlab = "full sweep", ylab = "-2 trimmed log-lik / n", main = "trimmed -2 log-likelihood")
  legend("topright", c("after row step", "after column step"), pch = c(1, 4), lty = 1,
         col = c("black", "steelblue"), bty = "n", cex = 0.8)

  md_panel(d_est,  sprintf("matrix MDs, estimated Sr, Sc (df = %d)", df))
  md_panel(d_true, sprintf("matrix MDs, true Sr, Sc (df = %d)", df))

  mtext(sprintf("%s, cellwise %s outliers (%.0f%% of cells) | n=%d p=%d q=%d | red = matrix with replaced cell",
                label, d$contam$cell_type, 100 * d$contam$eps, n, p, q), outer = TRUE, cex = 0.9)
  invisible(NULL)
}

## ---- Run: MLE, MMCD and trimmed MLE on the clean and the two cellwise data sets -------------
# TMLE = R1 + matrix reweighting, RAW = R1 without reweighting, RC = R1 + row / column reweighting (ECM).
sets <- list(clean = clean, cell_scale = dat_cell, cell_shift = dat_shft)
cut  <- qchisq(0.975, clean$p * clean$q)
res <- t(sapply(sets, function(d) {
  f1 <- fit_mle(d$X)
  f2 <- fit_mmcd(d$X)
  f3 <- fit_tmle(d$X, beta, maxit)
  f4 <- fit_tmle(d$X, beta, maxit, reweight = FALSE)
  f5 <- fit_tmle_rowcol(d$X, beta, maxit)
  c(KL_MLE    = kl_div(f1$Sr, f1$Sc, clean$Sr, clean$Sc),
    KL_MMCD   = kl_div(f2$Sr, f2$Sc, clean$Sr, clean$Sc),
    KL_TMLE   = kl_div(f3$Sr, f3$Sc, clean$Sr, clean$Sc),
    KL_RAW    = kl_div(f4$Sr, f4$Sc, clean$Sr, clean$Sc),
    KL_RC     = kl_div(f5$Sr, f5$Sc, clean$Sr, clean$Sc),
    flag_MLE  = sum(f1$md2 > cut),
    flag_MMCD = sum(f2$md2 > cut),
    flag_TMLE = sum(f3$md2 > cut),
    flag_RC   = sum(f5$md2 > cut),
    drop_RC   = f5$frac_dropped,
    iter_TMLE = f3$iter,
    conv_TMLE = f3$converged)
}))
print(round(res, 3))

## ---- Run: one figure per data set -----------------------------------------------------------
for (nm in names(sets)) {
  if (nm == "clean") next
  f <- fit_tmle(sets[[nm]]$X, beta, maxit)
  png(file.path("output", sprintf("tmle_%s.png", nm)), width = 1300, height = 1000, res = 110)
  plot_tmle_cell(sets[[nm]], f, clean$Sr, clean$Sc)
  dev.off()
  f <- fit_tmle_rowcol(sets[[nm]]$X, beta, maxit)
  png(file.path("output", sprintf("tmle_rowcol_%s.png", nm)), width = 1300, height = 1000, res = 110)
  plot_tmle_cell(sets[[nm]], f, clean$Sr, clean$Sc, label = "trimmed MLE + row/col reweighting (ECM)")
  dev.off()
}


#########################################################################################################
## ---- 3. Scan over the share of contaminated cells (eps) -----------------------------------------
# Median over reps of the KL to the truth; all methods are fitted on the SAME contaminated sample in every rep.
# frac_mat = share of matrices that contain a replaced cell. RAW / RC = R1 with beta = 0.5 without reweighting /
# with row / column reweighting (ECM).
scan_cell <- function(eps, type, betas = c(0.5, 0.75, 0.9), reps = 20, n = 100, p = 3, q = 4, sd_factor = 6) {
  r <- replicate(reps, {
    cl <- gen_matnorm(n, p, q, 0.6, 0.5)
    d  <- contaminate_cell(cl, eps, sd_factor, type)
    kl <- function(f) kl_div(f$Sr, f$Sc, cl$Sr, cl$Sc)
    c(frac_mat = mean(d$out),
      setNames(sapply(betas, function(b) kl(fit_tmle(d$X, b))), paste0("TMLE_", betas)),
      RAW_0.5 = kl(fit_tmle(d$X, 0.5, reweight = FALSE)), RC_0.5 = kl(fit_tmle_rowcol(d$X, 0.5)),
      MMCD = kl(fit_mmcd(d$X)), MLE = kl(fit_mle(d$X)))
  })
  apply(r, 1, median)
}

set.seed(11)
eps_grid <- c(0.01, 0.02, 0.05, 0.10)
for (type in c("scale", "shift")) {
  cat(sprintf("
=== cellwise %s outliers, 6 sd: median KL to truth over 20 reps ===
", type))
  tab <- cbind(eps = eps_grid, t(sapply(eps_grid, scan_cell, type = type)))
  print(round(tab, 3), row.names = FALSE)
}
