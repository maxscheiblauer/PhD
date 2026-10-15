# RJ, simplest version: one trimmed likelihood J*, minimized over (kept set, Sr, A). No reweighting.
# Same data and study as R1_simple.R; only the algorithm (section 2) is different.
# Conventions: X is p x q x n, vec(X_i) ~ N(vec M, Sc %x% Sr), Sr = p x p, Sc = q x q.
# Whitening matrix A (q x q) with A A' = Sc^-1; it is a free parameter, Sc = (A A')^-1.

library(robustmatrix)


## ---- 1. Data (copied from R1_simple.R) -------------------------------------------------------

gen_matnorm <- function(n, p, q, rho_r, rho_c, m = matrix(0, p, q)) {
  ar1 <- function(k, rho) rho^abs(outer(seq_len(k), seq_len(k), "-"))
  s_r <- ar1(p, rho_r)
  s_c <- ar1(q, rho_c)
  l_r <- t(chol(s_r))
  l_c <- t(chol(s_c))
  x <- vapply(seq_len(n), function(i) m + l_r %*% matrix(rnorm(p * q), p, q) %*% t(l_c), matrix(0, p, q))
  list(X = x, M = m, Sr = s_r, Sc = s_c, n = n, p = p, q = q)
}

# Replace a fraction eps of all n*p*q cells; sd_s = model sd of the cell.
#   scale: m + sd_factor * sd_s * N(0,1),  shift: m + sd_factor * sd_s + sd_s * N(0,1)
contaminate_cell <- function(clean, eps, sd_factor = 6, type = c("scale", "shift")) {
  type <- match.arg(type)
  p <- clean$p
  q <- clean$q
  n <- clean$n
  sd_s <- sqrt(outer(diag(clean$Sr), diag(clean$Sc)))
  mask <- array(FALSE, c(p, q, n))
  mask[sample.int(p * q * n, round(eps * p * q * n))] <- TRUE
  x <- clean$X
  for (i in which(apply(mask, 3, any))) {
    xi <- x[, , i]
    mk <- mask[, , i]
    z  <- matrix(rnorm(p * q), p, q)
    xi[mk] <- switch(type,
                     scale = (clean$M + sd_factor * sd_s * z)[mk],
                     shift = (clean$M + sd_factor * sd_s + sd_s * z)[mk])
    x[, , i] <- xi
  }
  clean$X <- x
  clean
}

# Divergence of (Sr, Sc) from the truth, eq. (31) of Mayrhofer et al. (2025); 0 = perfect.
kl_div <- function(s_r, s_c, s_r0, s_c0) {
  m_r <- solve(s_r0, s_r)
  m_c <- solve(s_c0, s_c)
  sum(diag(m_r)) * sum(diag(m_c)) - ncol(s_c) * determinant(m_r)$modulus -
    ncol(s_r) * determinant(m_c)$modulus - ncol(s_r) * ncol(s_c)
}


###################################################################################################
###################################################################################################
###################################################################################################
## ---- 2. RJ algorithm -------------------------------------------------------------------------
#
# Whitened columns: z_ij = X_i a_j (a_j = column j of A), distance d_ij = z_ij' Sr^-1 z_ij.
# Kept columns: keep[i, j] = TRUE. Quota: in every column position j exactly h = beta n of the n columns
# are kept (without it J* -> -Inf by dropping a whole position j and letting a_j grow).
#
# One objective (H = h q kept columns, each carrying log|Sr| + (p/q) log|Sc|, and log|Sc| = -2 log|det A|):
#   J* = H log|Sr| - 2 p h log|det A| + sum of the kept d_ij
# With beta = 1 this is -2 x the matrix-normal log-likelihood (up to constants).
#
# One sweep = three steps, each the exact minimizer of J* given the rest, so J* never goes up:
#   1. C-step : per position j keep the h smallest d_ij
#   2. Sr step: Sr = (sum of kept z_ij z_ij') / H
#   3. A step : column by column, a_j = sqrt(p h / (w' V^-1 w)) V^-1 w,
#               V = sum_{i kept in j} X_i' Sr^-1 X_i,  w = column j of A^-T
# Then one rescale: det(Sr) = 1, A divided by the same amount (product Sc %x% Sr and J* unchanged).
# Stop when the kept set does not change and J* drops by less than tol.
# End: Sc = (A A')^-1, multiplied by the consistency factor kappa(beta, p) (kept units are p-vectors).

keep_h <- function(d, h) order(d)[seq_len(h)]

# n x q matrix of the distances d_ij.
rj_dist <- function(x_c, a_mat, s_r) {
  o_r <- solve(s_r)
  t(apply(x_c, 3, function(x_i) {
    z <- x_i %*% a_mat                      # p x q, column j = z_ij
    colSums(z * (o_r %*% z))
  }))
}

# The objective J* for the current parameters and kept set.
rj_crit <- function(x_c, a_mat, s_r, keep) {
  p <- dim(x_c)[1]
  q <- dim(x_c)[2]
  h <- sum(keep[, 1])
  d <- rj_dist(x_c, a_mat, s_r)
  h * q * as.numeric(determinant(s_r)$modulus) - 2 * p * h * as.numeric(determinant(a_mat)$modulus) +
    sum(d[keep])
}

# Returns list(M, Sr, Sc, A, keep, crit (J* after every sweep), iter, converged).
fit_rj <- function(x, beta = 0.5, maxit = 500, tol = 1e-9) {
  p <- dim(x)[1]
  q <- dim(x)[2]
  n <- dim(x)[3]
  h <- floor(beta * n)
  m <- apply(x, 1:2, median)
  x_c <- sweep(x, 1:2, m)
  a_mat <- diag(q)
  s_r <- diag(p)
  keep <- matrix(FALSE, n, q)
  crit <- numeric(0)
  converged <- FALSE

  for (it in seq_len(maxit)) {
    keep_old <- keep

    # 1. C-step: per position j, keep the h smallest distances
    d <- rj_dist(x_c, a_mat, s_r)
    keep <- matrix(FALSE, n, q)
    for (j in seq_len(q)) keep[keep_h(d[, j], h), j] <- TRUE

    # 2. Sr step: covariance of the kept whitened columns
    s_r <- matrix(0, p, p)
    for (j in seq_len(q)) {
      for (i in which(keep[, j])) s_r <- s_r + tcrossprod(x_c[, , i] %*% a_mat[, j])
    }
    s_r <- s_r / (h * q)

    # 3. A step: one column at a time
    o_r <- solve(s_r)
    for (j in seq_len(q)) {
      v  <- Reduce(`+`, lapply(which(keep[, j]), function(i) crossprod(x_c[, , i], o_r %*% x_c[, , i])))
      w  <- t(solve(a_mat))[, j]
      vw <- solve(v, w)
      a_mat[, j] <- sqrt(p * h / sum(w * vw)) * vw
    }

    # Rescale: det(Sr) = 1
    s <- det(s_r)^(1 / p)
    s_r <- s_r / s
    a_mat <- a_mat / sqrt(s)

    crit[it] <- rj_crit(x_c, a_mat, s_r, keep)
    if (it > 1 && identical(keep, keep_old) && crit[it - 1] - crit[it] < tol * abs(crit[it])) {
      converged <- TRUE
      break
    }
  }

  kappa <- beta / pchisq(qchisq(beta, p), p + 2)
  list(M = m, Sr = s_r, Sc = kappa * solve(tcrossprod(a_mat)), A = a_mat, keep = keep,
       crit = crit, iter = it, converged = converged)
}
###################################################################################################
###################################################################################################
###################################################################################################


## ---- 3. Study: KL to the truth over the share of contaminated cells ---------------------------
# Every method is fitted on the SAME sample in each rep. MLE / MMCD from robustmatrix (MMCD alpha = 0.75,
# with its own reweighting). conv = share of reps that converged, iter = mean number of sweeps,
# up = mean number of sweeps in which J* went up (must be 0). KL columns: median, others: mean.

fit_mle <- function(x) {
  f <- robustmatrix::mmle(x, silent = TRUE)
  list(Sr = f$cov_row, Sc = f$cov_col)
}
fit_mmcd <- function(x) {
  f <- robustmatrix::mmcd(x, alpha = 0.75)
  list(Sr = f$cov_row, Sc = f$cov_col)
}

one_rep <- function(eps, type, n, p, q, betas) {
  cl <- gen_matnorm(n, p, q, 0.6, 0.5)
  d  <- if (eps > 0) contaminate_cell(cl, eps, 6, type) else cl
  kl <- function(f) kl_div(f$Sr, f$Sc, cl$Sr, cl$Sc)
  rj <- lapply(betas, function(b) fit_rj(d$X, b))
  c(MLE  = kl(fit_mle(d$X)),
    MMCD = kl(fit_mmcd(d$X)),
    setNames(sapply(rj, kl), paste0("RJ_", betas)),
    setNames(sapply(rj, `[[`, "converged"), paste0("conv_", betas)),
    setNames(sapply(rj, `[[`, "iter"), paste0("iter_", betas)),
    setNames(sapply(rj, function(f) sum(diff(f$crit) > 1e-8 * abs(f$crit[-1]))), paste0("up_", betas)))
}

scan_rj <- function(eps_grid, type, reps = 50, n = 100, p = 3, q = 4, betas = c(0.5, 0.75)) {
  t(sapply(eps_grid, function(eps) {
    r <- replicate(reps, one_rep(eps, type, n, p, q, betas))
    kl_rows <- grepl("^(MLE|MMCD|RJ)", rownames(r))
    c(eps = eps, apply(r[kl_rows, ], 1, median), rowMeans(r[!kl_rows, ]))
  }))
}

options(width = 200)
set.seed(2026)
eps_grid <- c(0, 0.02, 0.05, 0.08, 0.10)
for (type in c("scale", "shift")) {
  tab <- scan_rj(eps_grid, type)
  cat(sprintf("\n=== cellwise %s outliers, 6 sd | n=100 p=3 q=4 | median KL to truth over 50 reps ===\n", type))
  print(round(tab, 3), row.names = FALSE)
}
