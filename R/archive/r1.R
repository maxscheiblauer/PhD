# Lean R1 iteration for the convergence study (exp/step7-9).
# Same algorithm as tff(rule = "R1", kappa = FALSE, scale_id = "det", root_type = "sym"),
# but with O(1) cycle detection (hash of the kept sets) and per-half-step diagnostics.
#
# Criteria. Each factor's criterion gets its share of the Jacobian (as in RJ), so both are
# -2 x a trimmed Gaussian log-likelihood of the product with the same log-det part
#   B = hr log|Sr| + hc log|Sc|   (hr p = hc q, so B = beta n log|Sc %x% Sr|),
# and differ only in which distances are trimmed. Both are then profiled over the overall scale
# of the product (so the scale ping-pong between the half-steps does not count):
#   Lr~ = B + H log(s_d / H) + H,  s_d = sum of the hr smallest d_ij,  H = hr p
#   Lc~ = B + H log(s_e / H) + H,  s_e = sum of the hc smallest e_ik,  H = hc q
# Both are functions of the product Sc %x% Sr only. The row step never increases Lr~ for fixed Sc,
# the column step never increases Lc~ for fixed Sr (the added log|.| term of the other factor is
# constant in that step).

trim_sum <- function(dist, h) sum(sort.int(dist, partial = h)[seq_len(h)])

prof_crit <- function(Sr, Sc, hr, hc, s) {
  H <- hr * nrow(Sr)
  hr * logdet(Sr) + hc * logdet(Sc) + H * log(s / H) + H
}

# Both profiled criteria at (Sr, Sc).
r1_crit <- function(X, Sr, Sc, hr, hc) {
  D <- distances(X, Sr, Sc)
  c(Lr = prof_crit(Sr, Sc, hr, hc, trim_sum(D$d, hr)),
    Lc = prof_crit(Sr, Sc, hr, hc, trim_sum(D$e, hc)))
}

# One-step residual ("Nash gap" lower bound): what a C-step of each factor, from the same
# state and with the other factor held fixed, would still gain. >= 0, and 0 iff (Sr, Sc) is a
# fixed point of both C-steps (ties aside).
r1_gap <- function(X, Sr, Sc, hr, hc) {
  n <- dim(X)[3]; p <- nrow(Sr); q <- nrow(Sc)
  Z <- zmat(X, root(Sc)); d <- qdist(Z, Sr)
  Srn <- tcrossprod(Z[, csel(matrix(d, n), hr), drop = FALSE]) / hr
  U <- umat(X, t(root(Sr))); e <- qdist(U, Sc)
  Scn <- tcrossprod(U[, csel(matrix(e, n), hc), drop = FALSE]) / hc
  c(gap_r = prof_crit(Sr, Sc, hr, hc, trim_sum(d, hr)) - prof_crit(Srn, Sc, hr, hc, trim_sum(qdist(Z, Srn), hr)),
    gap_c = prof_crit(Sr, Sc, hr, hc, trim_sum(e, hc)) - prof_crit(Sr, Scn, hr, hc, trim_sum(qdist(U, Scn), hc)))
}

# R1 from (Sr, Sc). diag = TRUE records Lr~, Lc~ after every half-step (row 1 = start,
# then row step, column step, row step, ...) and the one-step residual r1_gap after every
# iteration. Always recorded: the shape-KL between successive products (step).
r1_fit <- function(X, beta = 0.75, Sr, Sc, maxit = 500, tol = 1e-10, diag = TRUE,
                   hr = round(beta * dim(X)[3] * dim(X)[2]), hc = round(beta * dim(X)[3] * dim(X)[1])) {
  p <- dim(X)[1]; q <- dim(X)[2]; n <- dim(X)[3]
  s <- scale_ident(Sr, Sc, "det"); Sr <- s$Sr; Sc <- s$Sc
  seen <- new.env(hash = TRUE); hist <- vector("list", maxit)
  crit <- if (diag) matrix(NA_real_, 2 * maxit + 1, 2, dimnames = list(NULL, c("Lr", "Lc")))
  if (diag) crit[1, ] <- r1_crit(X, Sr, Sc, hr, hc)
  step <- nchg <- gap <- rep(NA_real_, maxit)
  a_old <- b_old <- NULL; status <- "maxit"; period <- NA
  for (k in seq_len(maxit)) {
    Sr0 <- Sr; Sc0 <- Sc
    Z <- zmat(X, root(Sc)); a <- csel(matrix(qdist(Z, Sr), n), hr)
    Sr <- tcrossprod(Z[, a, drop = FALSE]) / hr
    if (diag) crit[2 * k, ] <- r1_crit(X, Sr, Sc, hr, hc)
    U <- umat(X, t(root(Sr))); b <- csel(matrix(qdist(U, Sc), n), hc)
    Sc <- tcrossprod(U[, b, drop = FALSE]) / hc
    s <- scale_ident(Sr, Sc, "det"); Sr <- s$Sr; Sc <- s$Sc
    if (!all(is.finite(Sr)) || !all(is.finite(Sc)) ||
        min(eigen(Sr, TRUE, TRUE)$values, eigen(Sc, TRUE, TRUE)$values) < 1e-12) {
      status <- "singular"; break
    }
    if (diag) { crit[2 * k + 1, ] <- r1_crit(X, Sr, Sc, hr, hc); gap[k] <- sum(r1_gap(X, Sr, Sc, hr, hc)) }
    step[k] <- kl_shape(Sr, Sc, Sr0, Sc0)
    nchg[k] <- if (is.null(a_old)) NA else sum(a != a_old) + sum(b != b_old)
    a_old <- a; b_old <- b
    hist[[k]] <- list(Sr = Sr, Sc = Sc)
    if (k > 1 && nchg[k] == 0 && step[k] < tol) { status <- "converged"; break }
    key <- paste(c(which(a), 0L, which(b)), collapse = ",")
    k0 <- seen[[key]]
    if (!is.null(k0) && k0 < k - 1 && nchg[k] > 0 &&
        kl_shape(Sr, Sc, hist[[k0]]$Sr, hist[[k0]]$Sc) < 1e-12) {
      status <- "cycle"; period <- k - k0; break
    }
    assign(key, k, envir = seen)
  }
  list(Sr = Sr, Sc = Sc, a = a, b = b, status = status, iter = k, period = period,
       n_sets = length(ls(seen)), step = step[seq_len(k)], nchange = nchg[seq_len(k)],
       crit = if (diag) crit[seq_len(2 * k + 1), , drop = FALSE], gap = if (diag) gap[seq_len(k)],
       hr = hr, hc = hc)
}

# Local contraction at a fixed point: Jacobian of Sc -> Sc+ (one full R1 iteration, det Sc = 1)
# with the kept sets a, b frozen, by central differences in the trace-free log-coordinates of Sc.
# Returns the eigenvalues; rho = max |eigenvalue| < 1 means the fixed point is locally attracting.
r1_jacobian <- function(X, Sc, a, b, hr = sum(a), hc = sum(b), eps = 1e-6) {
  q <- nrow(Sc)
  Tmap <- function(Sc) {
    Z <- zmat(X, root(Sc)); Sr <- tcrossprod(Z[, a, drop = FALSE]) / hr
    U <- umat(X, t(root(Sr))); Sc <- tcrossprod(U[, b, drop = FALSE]) / hc
    Sc / exp(logdet(Sc) / q)
  }
  ev <- eigen(Sc, TRUE); L <- ev$vectors %*% (log(ev$values) * t(ev$vectors))
  mexp <- function(S) { e <- eigen(S, TRUE); e$vectors %*% (exp(e$values) * t(e$vectors)) }
  mlog <- function(S) { e <- eigen(S, TRUE); e$vectors %*% (log(e$values) * t(e$vectors)) }
  # orthonormal basis of trace-free symmetric matrices
  B <- list()
  for (i in seq_len(q)) for (j in i:q) {
    E <- matrix(0, q, q); E[i, j] <- E[j, i] <- if (i == j) 1 else 1 / sqrt(2)
    B[[length(B) + 1]] <- E
  }
  B <- lapply(B, function(E) E - diag(sum(diag(E)) / q, q))
  M <- sapply(B, as.vector); Q <- qr.Q(qr(M))[, seq_len(length(B) - 1), drop = FALSE]
  J <- sapply(seq_len(ncol(Q)), function(k) {
    E <- matrix(Q[, k], q)
    (as.vector(mlog(Tmap(mexp(L + eps * E)))) - as.vector(mlog(Tmap(mexp(L - eps * E))))) / (2 * eps)
  })
  eigen(crossprod(Q, J), only.values = TRUE)$values
}
