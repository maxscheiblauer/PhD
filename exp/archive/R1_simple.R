# R1, simplest version: slice-trimmed flip-flop = flip step, flop step, rescale. No reweighting.
# Multi-start: n_start starts, keep the one with the smallest trimmed criterion (see fit_r1()).
# Compared with the plain matrix-normal MLE and MMCD under cellwise contamination. Base R + robustmatrix.
# Conventions as in R1_cellwise.R: X is p x q x n, vec(X_i) ~ N(vec M, Sc %x% Sr), Sr = p x p, Sc = q x q.
# Whitening always uses the symmetric inverse square root (isqrt), never another root.

library(robustmatrix)


## ---- 1. Data (copied from R1_cellwise.R) -----------------------------------------------------

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
# Depends on the product Sc %x% Sr only.
kl_div <- function(s_r, s_c, s_r0, s_c0) {
  m_r <- solve(s_r0, s_r)
  m_c <- solve(s_c0, s_c)
  sum(diag(m_r)) * sum(diag(m_c)) - ncol(s_c) * determinant(m_r)$modulus -
    ncol(s_r) * determinant(m_c)$modulus - ncol(s_r) * ncol(s_c)
}


## ---- 2. R1: flip, flop, rescale ---------------------------------------------------------------

# Symmetric inverse square root S^(-1/2).
isqrt <- function(s) {
  e <- eigen(s, symmetric = TRUE)
  e$vectors %*% (t(e$vectors) / sqrt(e$values))
}

keep_h <- function(d, h) order(d)[seq_len(h)]

# One run from a given start s_c. M fixed (x_c = X - M already centred).
#   flip: z_ij = column j of X_i Sc^(-1/2) (p-vectors), keep the h_r = beta n q smallest d_ij, Sr = their cov
#   flop: u_ik = row k of Sr^(-1/2) X_i  (q-vectors), keep the h_c = beta n p smallest e_ik, Sc = their cov
# Converged when both kept sets stop changing. After a flop step mean(e kept) = q, so the trimmed
# -2 log-lik is q log|Sr| + p log|Sc| + pq = log|Sc %x% Sr| + pq: `crit` = log|Sc %x% Sr| (MCD-like).
r1_run <- function(x_c, s_c, beta, maxit) {
  p <- dim(x_c)[1]
  q <- dim(x_c)[2]
  n <- dim(x_c)[3]
  hr <- floor(beta * n * q)
  hc <- floor(beta * n * p)
  s_r <- diag(p)
  a <- b <- integer(0)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    a_old <- a
    b_old <- b

    # flip
    w_c <- isqrt(s_c)
    z   <- matrix(apply(x_c, 3, function(x_i) x_i %*% w_c), p)
    a   <- sort(keep_h(colSums(z * solve(s_r, z)), hr))
    s_r <- tcrossprod(z[, a, drop = FALSE]) / hr

    # flop
    w_r <- isqrt(s_r)
    u   <- matrix(apply(x_c, 3, function(x_i) t(w_r %*% x_i)), q)
    b   <- sort(keep_h(colSums(u * solve(s_c, u)), hc))
    s_c <- tcrossprod(u[, b, drop = FALSE]) / hc

    if (identical(a, a_old) && identical(b, b_old)) {
      converged <- TRUE
      break
    }
  }
  crit <- q * as.numeric(determinant(s_r)$modulus) + p * as.numeric(determinant(s_c)$modulus)
  list(Sr = s_r, Sc = s_c, crit = crit, iter = it, converged = converged)
}

# Multi-start R1. M = cellwise median (fixed). Starts for Sc:
#   start 1: identity;  starts 2..n_start: crude Sc from a random subset of n_sub matrices,
#   sum (X_i - M)'(X_i - M) / (p n_sub)  (i.e. assuming Sr = I).
# Keep the run with the smallest crit, then rescale the product by c = beta / F_{chi2_{q+2}}(chi2_{beta; q})
# (the last step is a flop step on chi2_q distances). Returns list(M, Sr, Sc, crit, converged, n_conv, n_sol,
# single): n_sol = number of distinct end points (crit rounded to 1e-6), single = rescaled fit of start 1 alone.
fit_r1 <- function(x, beta = 0.5, n_start = 10, n_sub = 10, maxit = 100) {
  p <- dim(x)[1]
  q <- dim(x)[2]
  n <- dim(x)[3]
  m <- apply(x, 1:2, median)
  x_c <- sweep(x, 1:2, m)
  starts <- c(list(diag(q)), lapply(seq_len(n_start - 1), function(k) {
    idx <- sample.int(n, n_sub)
    Reduce(`+`, lapply(idx, function(i) crossprod(x_c[, , i]))) / (p * n_sub)
  }))
  runs <- lapply(starts, function(s_c) r1_run(x_c, s_c, beta, maxit))
  crit <- sapply(runs, `[[`, "crit")
  best <- runs[[which.min(crit)]]
  c_fac <- beta / pchisq(qchisq(beta, q), q + 2)
  list(M = m, Sr = best$Sr, Sc = c_fac * best$Sc, crit = best$crit, converged = best$converged,
       n_conv = sum(sapply(runs, `[[`, "converged")), n_sol = length(unique(round(crit, 6))),
       single = list(Sr = runs[[1]]$Sr, Sc = c_fac * runs[[1]]$Sc))
}


## ---- 3. Study: KL to the truth over the share of contaminated cells ---------------------------
# Every method is fitted on the SAME sample in each rep. MLE / MMCD from robustmatrix (MMCD alpha = 0.75,
# with its own reweighting). conv = share of reps whose chosen R1 start converged (kept sets fixed),
# starts_conv = mean number of the n_start starts that converged, n_sol = mean number of distinct end points,
# R1single = R1 from the identity start only (to see what the restarts add). KL columns: median, others: mean.

fit_mle <- function(x) {
  f <- robustmatrix::mmle(x, silent = TRUE)
  list(Sr = f$cov_row, Sc = f$cov_col)
}
fit_mmcd <- function(x) {
  f <- robustmatrix::mmcd(x, alpha = 0.75)
  list(Sr = f$cov_row, Sc = f$cov_col)
}

one_rep <- function(eps, type, n, p, q, betas, n_start) {
  cl <- gen_matnorm(n, p, q, 0.6, 0.5)
  d  <- if (eps > 0) contaminate_cell(cl, eps, 6, type) else cl
  kl <- function(f) kl_div(f$Sr, f$Sc, cl$Sr, cl$Sc)
  r1 <- lapply(betas, function(b) fit_r1(d$X, b, n_start))
  c(MLE  = kl(fit_mle(d$X)),
    MMCD = kl(fit_mmcd(d$X)),
    setNames(sapply(r1, kl), paste0("R1_", betas)),
    setNames(sapply(r1, function(f) kl(f$single)), paste0("R1single_", betas)),
    setNames(sapply(r1, `[[`, "converged"), paste0("conv_", betas)),
    setNames(sapply(r1, `[[`, "n_conv"), paste0("starts_conv_", betas)),
    setNames(sapply(r1, `[[`, "n_sol"), paste0("n_sol_", betas)))
}

scan_r1 <- function(eps_grid, type, reps = 50, n = 100, p = 3, q = 4, betas = c(0.5, 0.75), n_start = 10) {
  t(sapply(eps_grid, function(eps) {
    r <- replicate(reps, one_rep(eps, type, n, p, q, betas, n_start))
    kl_rows <- grepl("^(MLE|MMCD|R1)", rownames(r))
    c(eps = eps, apply(r[kl_rows, ], 1, median), rowMeans(r[!kl_rows, ]))
  }))
}

options(width = 200)
set.seed(2026)
eps_grid <- c(0, 0.02, 0.05, 0.08, 0.10)
for (type in c("scale", "shift")) {
  tab <- scan_r1(eps_grid, type)
  cat(sprintf("\n=== cellwise %s outliers, 6 sd | n=100 p=3 q=4 | median KL to truth over 50 reps ===\n", type))
  print(round(tab, 3), row.names = FALSE)
}
