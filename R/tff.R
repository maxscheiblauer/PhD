# Column-trimmed flip-flop (spec notes/trimmed_flip_flop_spec.md).
#
# Rules: R0, R1, R2 (= MMCD-type, observation trimming), R3 (spec section 5),
#        RJ  (joint trimmed likelihood, Step 4 alternative, see below),
#        none (no trimming; equals plain flip-flop).
#
# RJ: parameters (Sr, A) with Sc = (A A')^{-1}. Criterion
#   J*(Sr, A; a) = sum_ij a_ij [ log|Sr| + (p/q) log|Sc| + d_ij ],  d_ij = z_ij' Sr^{-1} z_ij,  z_ij = X_i A e_j.
# With a == 1 this is -2 x the matrix-normal log-likelihood (the Jacobian |A|^p of X_i -> Z_i is
# split equally over the q columns). a-step (C-step), Sr-step (closed form) and A-step (exact
# column-wise coordinate descent) each decrease J*, so J* is monotone by construction.
#   RJu: C-step over all n q distances. If every z_ij of one column index j gets trimmed,
#        J* is unbounded below (alpha_j -> infinity); the run is then stopped as "degenerate".
#   RJ:  balanced C-step, keeps the round(beta n) smallest d_ij within each column index j.

# Whitened columns: p x (n q) matrix, column i + (j-1) n holds z_ij.
zmat <- function(X, Ac) {
  p <- dim(X)[1]; q <- dim(X)[2]; n <- dim(X)[3]
  matrix(matrix(aperm(X, c(1, 3, 2)), p * n, q) %*% Ac, p, n * q)
}

# Whitened rows: q x (n p) matrix, column i + (k-1) n holds u_ik.
umat <- function(X, Ar) {
  p <- dim(X)[1]; q <- dim(X)[2]; n <- dim(X)[3]
  U <- array(Ar %*% matrix(X, p, q * n), c(p, q, n))
  matrix(aperm(U, c(2, 3, 1)), q, n * p)
}

# Squared distances of the columns of a matrix W w.r.t. S.
qdist <- function(W, S) colSums(W * solve(S, W))

csel <- function(d, h) {
  a <- matrix(FALSE, nrow(d), ncol(d))
  a[order(d)[seq_len(h)]] <- TRUE
  a
}

logdet <- function(S) as.numeric(determinant(S)$modulus)

# Distances d (n x q) and e (n x p) at (Sr, Sc) with the given roots.
distances <- function(X, Sr, Sc, root_type = "sym", rot = NULL, Ac = NULL) {
  n <- dim(X)[3]
  if (is.null(Ac)) Ac <- root(Sc, root_type, rot$c)
  Ar <- t(root(Sr, root_type, rot$r))
  list(d = matrix(qdist(zmat(X, Ac), Sr), n),
       e = matrix(qdist(umat(X, Ar), Sc), n))
}

init_elemental <- function(X, m = NULL) {
  p <- dim(X)[1]; q <- dim(X)[2]; n <- dim(X)[3]
  if (is.null(m)) m <- max(2, ceiling(max(p, q) / min(p, q)) + 1)
  for (k in 1:20) {
    f <- try(flipflop(X[, , sample(n, m), drop = FALSE], maxit = 50, tol = 1e-8), silent = TRUE)
    if (!inherits(f, "try-error") && all(is.finite(f$Sr)) &&
        min(eigen(f$Sc %x% f$Sr, TRUE, TRUE)$values) > 1e-10) return(f[c("Sr", "Sc")])
  }
  stop("no non-singular elemental start found")
}

make_init <- function(X, init) {
  if (is.list(init)) return(init)
  p <- dim(X)[1]; q <- dim(X)[2]
  switch(init,
    identity = list(Sr = diag(p), Sc = diag(q)),
    ff = flipflop(X, maxit = 500, tol = 1e-10)[c("Sr", "Sc")],
    elemental = init_elemental(X),
    stop("unknown init"))
}

# Exact coordinate descent on the columns of A for
#   f(A) = -(2 h p / q) log|det A| + sum_j alpha_j' V_j alpha_j.
rj_update_A <- function(A, V, h, p, q, sweeps = 20, tol = 1e-12) {
  cst <- 2 * h * p / q
  for (s in seq_len(sweeps)) {
    A_old <- A
    for (j in seq_len(q)) {
      cj <- solve(A)[j, ]                    # det A is linear in alpha_j along this vector
      Vc <- solve(V[, , j], cj)
      A[, j] <- sqrt(cst / (2 * sum(cj * Vc))) * Vc
    }
    if (max(abs(A - A_old)) < tol * max(abs(A))) break
  }
  A
}

tff <- function(X, rule = c("R1", "R0", "R2", "R3", "RJ", "RJu", "none"), beta = 0.75,
                hr = NULL, hc = NULL, denom = c("h", "n"), kappa = TRUE,
                scale_id = c("det", "none", "tr", "c11"),
                root_type = c("sym", "chol", "randrot"), rot = NULL,
                init = "ff", t_flag = NULL, Sc_fixed = NULL,
                maxit = 200, tol = 1e-8, trace = TRUE) {
  rule <- match.arg(rule); denom <- match.arg(denom)
  scale_id <- match.arg(scale_id); root_type <- match.arg(root_type)
  p <- dim(X)[1]; q <- dim(X)[2]; n <- dim(X)[3]
  if (rule == "none") beta <- 1
  if (is.null(hr)) hr <- round(beta * n * q)
  if (is.null(hc)) hc <- round(beta * n * p)
  h_obs <- round(beta * n)                      # R2
  if (is.null(t_flag)) t_flag <- floor(q * (1 - beta))   # R3: flag X_i if > t_flag columns trimmed
  if (root_type == "randrot" && is.null(rot)) rot <- list(c = rorth(q), r = rorth(p))
  if (!is.null(Sc_fixed)) scale_id <- "none"

  st <- make_init(X, init)
  Sr <- st$Sr; Sc <- if (is.null(Sc_fixed)) st$Sc else Sc_fixed
  s <- scale_ident(Sr, Sc, if (is.null(Sc_fixed)) scale_id else "none"); Sr <- s$Sr; Sc <- s$Sc
  is_rj <- rule %in% c("RJ", "RJu")
  if (is_rj) A <- root(Sc, root_type, rot$c)

  S_old <- Sc %x% Sr
  a_old <- b_old <- matrix(NA, n, 1)
  Jv <- delta <- nchg <- rep(NA_real_, maxit)
  keys <- character(maxit); hist <- vector("list", maxit)
  a_tr <- b_tr <- if (trace) vector("list", maxit) else NULL
  status <- "maxit"; period <- NA

  for (k in seq_len(maxit)) {
    ## ---- Sigma^r step
    Ac <- if (is_rj) A else root(Sc, root_type, rot$c)
    Z <- zmat(X, Ac)
    d <- matrix(qdist(Z, Sr), n)
    a <- switch(rule,
      none = matrix(TRUE, n, q),
      R2 = { sel <- matrix(FALSE, n, q); sel[order(rowSums(d))[seq_len(h_obs)], ] <- TRUE; sel },
      RJ = apply(d, 2, function(dj) rank(dj, ties.method = "first") <= h_obs),  # h_obs per column index
      csel(d, hr))
    ha <- sum(a)
    Nr <- if (denom == "h") ha else n * q
    kr <- if (kappa && rule %in% c("R0", "R1", "R3")) kappa_mcd(ha / (n * q), p) else 1
    Za <- Z[, as.vector(a), drop = FALSE]
    Sr <- kr / Nr * tcrossprod(Za)

    ## ---- Sigma^c step
    b <- matrix(TRUE, n, p)
    if (!is.null(Sc_fixed)) {
      # oracle: Sc stays at the truth
    } else if (is_rj) {
      Sri <- solve(Sr)
      V <- array(0, c(q, q, q))
      for (i in seq_len(n)) {
        Wi <- crossprod(X[, , i], Sri %*% X[, , i])
        for (j in which(a[i, ])) V[, , j] <- V[, , j] + Wi
      }
      # RJu: an entirely trimmed column index j makes V_j singular and J* unbounded below
      if (min(apply(V, 3, function(v) min(eigen(v, TRUE, TRUE)$values))) <= 1e-10) {
        status <- "degenerate"; break
      }
      A <- rj_update_A(A, V, ha, p, q)
      Sc <- solve(tcrossprod(A))
    } else {
      Ar <- t(root(Sr, root_type, rot$r))
      U <- umat(X, Ar)
      e <- matrix(qdist(U, Sc), n)
      b <- switch(rule,
        none = , R0 = matrix(TRUE, n, p),
        R1 = csel(e, hc),
        R2 = matrix(a[, 1], n, p),
        R3 = matrix(rowSums(!a) <= t_flag, n, p))
      hb <- sum(b)
      if (hb < q) { status <- "degenerate"; break }
      Nc <- if (denom == "h") hb else n * p
      kc <- if (!kappa) 1 else switch(rule,
        R1 = , R3 = kappa_mcd(hb / (n * p), q),
        R2 = kappa_mcd(h_obs / n, p * q),     # MMCD: consistency on the product, dim pq
        1)
      Sc <- kc / Nc * tcrossprod(U[, as.vector(b), drop = FALSE])
    }
    # RJ: kappa(beta, p) would cancel inside the alternation, so it is applied once at the end.

    ## ---- scale identification
    s <- scale_ident(Sr, Sc, scale_id); Sr <- s$Sr; Sc <- s$Sc
    if (is_rj) A <- A * sqrt(s$c)

    # per factor: with scale_id = "none" the factor scales may drift apart geometrically
    rc <- sapply(list(Sr, Sc), function(M) { ev <- eigen(M, TRUE, TRUE)$values; min(ev) / max(ev) })
    if (!all(is.finite(rc)) || min(rc) <= 1e-12) { status <- "singular"; break }

    ## ---- objective, change, bookkeeping
    if (is_rj) {
      dd <- matrix(qdist(zmat(X, A), Sr), n)
      Jv[k] <- ha * logdet(Sr) + ha * p / q * logdet(Sc) + sum(dd[a])
    } else {
      D <- distances(X, Sr, Sc, root_type, rot)
      Jv[k] <- ha * logdet(Sr) + sum(D$d[a]) + sum(b) * logdet(Sc) + sum(D$e[b])
    }
    S_new <- Sc %x% Sr
    delta[k] <- norm(S_new - S_old, "F") / norm(S_old, "F")
    nchg[k] <- if (k == 1) NA else sum(a != a_old) + sum(b != b_old)
    keys[k] <- paste0(as.integer(c(a, b)), collapse = "")
    hist[[k]] <- S_new
    if (trace) { a_tr[[k]] <- a; b_tr[[k]] <- b }
    S_old <- S_new; a_old <- a; b_old <- b

    if (k > 1 && delta[k] < tol && nchg[k] == 0) { status <- "converged"; break }
    # exact limit cycle: same weights and same estimate as L >= 2 iterations ago
    if (k > 2) {
      prev <- which(keys[seq_len(k - 2)] == keys[k])
      prev <- prev[prev < k - 1 & keys[k - 1] != keys[k]]
      for (k0 in rev(prev)) {
        if (norm(S_new - hist[[k0]], "F") / norm(S_new, "F") < 1e-6) {
          status <- "cycle"; period <- k - k0; break
        }
      }
      if (status == "cycle") break
    }
  }
  kk <- k
  if (is_rj && kappa) Sr <- Sr * kappa_mcd(sum(a) / (n * q), p)
  revisit <- any(sapply(seq_len(kk)[-1], function(m)
    keys[m] != keys[m - 1] && keys[m] %in% keys[seq_len(max(0, m - 2))]))
  list(Sr = Sr, Sc = Sc, a = a, b = b, iter = kk, status = status,
       converged = status == "converged", period = period, revisit = revisit,
       J = Jv[seq_len(kk)], delta = delta[seq_len(kk)], nchange = nchg[seq_len(kk)],
       a_trace = if (trace) a_tr[seq_len(kk)], b_trace = if (trace) b_tr[seq_len(kk)],
       rule = rule, beta = beta, root_type = root_type, rot = rot, A = if (is_rj) A)
}

# Several starts; returns all fits and the index of the one with the smallest final J.
tff_multi <- function(X, starts = c("identity", "ff", rep("elemental", 8)), ...) {
  fits <- lapply(starts, function(s) {
    f <- try(tff(X, init = s, trace = FALSE, ...), silent = TRUE)
    if (inherits(f, "try-error")) NULL else f
  })
  ok <- which(!sapply(fits, is.null))
  Jfin <- sapply(fits, function(f) if (is.null(f) || f$status %in% c("singular", "degenerate"))
    Inf else tail(f$J, 1))
  list(fits = fits, best = which.min(Jfin), Jfinal = Jfin)
}
