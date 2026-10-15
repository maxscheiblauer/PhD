# Core matrix-normal utilities: sampler, plain flip-flop MLE, roots, kappa.
# Data convention throughout: X is a p x q x n array (same as robustmatrix).

# Sampler: vec(X_i) ~ N(vec(M), Sc %x% Sr), X_i = M + Lr E Rc with Lr Lr' = Sr, Rc'Rc = Sc.
rmatnorm <- function(n, M = NULL, Sr, Sc) {
  p <- nrow(Sr); q <- nrow(Sc)
  if (is.null(M)) M <- matrix(0, p, q)
  Lr <- t(chol(Sr)); Rc <- chol(Sc)
  X <- array(rnorm(p * q * n), c(p, q, n))
  for (i in seq_len(n)) X[, , i] <- M + Lr %*% X[, , i] %*% Rc
  X
}

# Scale identification of (Sr, Sc); the Kronecker product is unchanged.
scale_ident <- function(Sr, Sc, scale_id = c("det", "none", "tr", "c11")) {
  scale_id <- match.arg(scale_id)
  q <- nrow(Sc)
  cc <- switch(scale_id,
    none = 1,
    det = exp(as.numeric(determinant(Sc)$modulus) / q),
    tr = sum(diag(Sc)) / q,
    c11 = Sc[1, 1])
  list(Sr = Sr * cc, Sc = Sc / cc, c = cc)
}

# Plain flip-flop MLE. center = TRUE estimates M by the sample mean (as robustmatrix::mmle).
flipflop <- function(X, center = FALSE, init = NULL, maxit = 500, tol = 1e-10,
                     scale_id = "c11") {
  p <- dim(X)[1]; q <- dim(X)[2]; n <- dim(X)[3]
  M <- matrix(0, p, q)
  if (center) {
    M <- apply(X, 1:2, mean)
    X <- sweep(X, 1:2, M)
  }
  Sr <- if (is.null(init)) diag(p) else init$Sr
  Sc <- if (is.null(init)) diag(q) else init$Sc
  S_old <- Sc %x% Sr
  conv <- FALSE
  for (it in seq_len(maxit)) {
    # sum_i X_i Sc^{-1} X_i' = Z Z' with Z the whitened columns (zmat/umat in R/tff.R)
    Sr <- tcrossprod(zmat(X, root(Sc, "chol"))) / (n * q)
    Sc <- tcrossprod(umat(X, t(root(Sr, "chol")))) / (n * p)
    s <- scale_ident(Sr, Sc, scale_id); Sr <- s$Sr; Sc <- s$Sc
    S_new <- Sc %x% Sr
    rel <- norm(S_new - S_old, "F") / norm(S_old, "F")
    S_old <- S_new
    if (rel < tol) { conv <- TRUE; break }
  }
  list(Sr = Sr, Sc = Sc, M = M, iter = it, converged = conv)
}

# Random orthogonal matrix (Haar).
rorth <- function(m) {
  qrd <- qr(matrix(rnorm(m * m), m))
  Q <- qr.Q(qrd)
  Q %*% diag(sign(diag(qr.R(qrd))), m)
}

# Returns A with A A' = S^{-1}.
#   sym:     S^{-1/2}
#   chol:    lower-triangular L with L L' = S^{-1}
#   randrot: S^{-1/2} Q for a fixed orthogonal Q (drawn if not given)
root <- function(S, type = c("sym", "chol", "randrot"), Q = NULL) {
  type <- match.arg(type)
  if (type == "chol") return(t(chol(solve(S))))
  ev <- eigen(S, symmetric = TRUE)
  A <- ev$vectors %*% (t(ev$vectors) / sqrt(ev$values))
  if (type == "randrot") {
    if (is.null(Q)) Q <- rorth(nrow(S))
    A <- A %*% Q
  }
  A
}

# MCD-type consistency factor kappa(beta, m) = beta / F_{chi2_{m+2}}(chi2_{m, beta}).
# (Named kappa_mcd because base::kappa already exists.)
kappa_mcd <- function(beta, m) {
  if (beta >= 1) return(1)
  beta / pchisq(qchisq(beta, m), m + 2)
}
