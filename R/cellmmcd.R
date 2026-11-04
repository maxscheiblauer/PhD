# cellMMCD: cellwise matrix MCD. Trimmed observed-cell likelihood with Kronecker covariance.
# Theory and derivations: notes/cellmmcd/JOURNAL.md (2026-11-04).
# Data X: array p x q x n. A cell is addressed by its vec index c = k + (l - 1) p.
# W: logical array p x q x n, TRUE = cell kept (observed), FALSE = flagged.

## ---- 1. Pairwise Kronecker start (no iteration, no imputation) -------------------
# Scale per cell position by the MAD; log s^2_kl = a_k + b_l fitted by median polish;
# correlations by the Gnanadesikan-Kettenring identity on standardized cells, using that
# corr(X_kl, X_k'l) = rho^r_kk' for every l (median over l), and likewise for columns.
kpair <- function(X, M = NULL, pd_min = 1e-3) {
  p <- dim(X)[1]; q <- dim(X)[2]
  if (is.null(M)) M <- apply(X, c(1, 2), median)
  Xc <- sweep(X, c(1, 2), M)
  s <- apply(Xc, c(1, 2), mad)
  if (any(s <= 0)) s[s <= 0] <- min(s[s > 0])
  L <- log(s^2); a <- rep(0, p); b <- rep(0, q)
  for (it in 1:3) {
    a <- a + apply(L - outer(a, b, "+"), 1, median)
    b <- b + apply(L - outer(a, b, "+"), 2, median)
  }
  U <- sweep(Xc, c(1, 2), s, "/")
  gk <- function(u, v) { sp <- mad(u + v)^2; sm <- mad(u - v)^2; (sp - sm) / (sp + sm) }
  Rr <- diag(p); Rc <- diag(q)
  if (p > 1) for (k in 1:(p - 1)) for (k2 in (k + 1):p)
    Rr[k, k2] <- Rr[k2, k] <- median(vapply(1:q, function(l) gk(U[k, l, ], U[k2, l, ]), 0))
  if (q > 1) for (l in 1:(q - 1)) for (l2 in (l + 1):q)
    Rc[l, l2] <- Rc[l2, l] <- median(vapply(1:p, function(k) gk(U[k, l, ], U[k, l2, ]), 0))
  pdfix <- function(R) {
    e <- eigen(R, symmetric = TRUE); v <- pmax(e$values, pd_min * max(e$values))
    R2 <- e$vectors %*% (v * t(e$vectors)); d <- sqrt(diag(R2)); R2 / outer(d, d)
  }
  dr <- exp(a / 2); dc <- exp(b / 2)
  list(Sr = dr * pdfix(Rr) * rep(dr, each = p), Sc = dc * pdfix(Rc) * rep(dc, each = q), M = M, s = s)
}

## ---- 2. Objective via the Schur identity ------------------------------------------
# Per observation: d_i = q log|Sr| + p log|Sc| + log det O_MM + tr(Or Y Oc Y') - u' O_MM^{-1} u
#                        + |O_i| log 2pi + sum lambda[flagged].   Xc is centered.
# Observation weights v (block 2): F = sum_i v_i d_i + mu * sum_i (1 - v_i). v = NULL: all kept.
cm_obj_i <- function(Xc, Sr, Sc, W, lambda) {
  p <- dim(Xc)[1]; q <- dim(Xc)[2]; n <- dim(Xc)[3]
  Or <- solve(Sr); Oc <- solve(Sc); base <- q * logdet(Sr) + p * logdet(Sc)
  d <- numeric(n)
  for (i in 1:n) {
    Y <- Xc[, , i]; w <- W[, , i]; Y[!w] <- 0
    Tm <- Or %*% Y %*% Oc
    val <- base + sum(Tm * Y) + sum(w) * log(2 * pi) + sum(lambda[!w])
    Mi <- which(!w)
    if (length(Mi)) {
      K <- (Mi - 1) %% p + 1; L <- (Mi - 1) %/% p + 1
      Omm <- Or[K, K, drop = FALSE] * Oc[L, L, drop = FALSE]
      u <- Tm[Mi]
      val <- val + logdet(Omm) - sum(u * solve(Omm, u))
    }
    d[i] <- val
  }
  d
}
cm_obj <- function(Xc, Sr, Sc, W, lambda, v = NULL, mu = 0) {
  d <- cm_obj_i(Xc, Sr, Sc, W, lambda)
  if (is.null(v)) sum(d) else sum(d[v]) + if (is.finite(mu)) mu * sum(!v) else 0
}

# Delta_c = -2 log f(x_c | other observed cells) for cell c of one observation, with flagged set Mi
# (c may or may not be in Mi). Tf = Or Xc_i Oc on the full matrix.
cm_delta1 <- function(Tf, Or, Oc, Mi, c, p) {
  S <- c(setdiff(Mi, c), c); jc <- length(S)
  if (jc == 1) {
    k <- (c - 1) %% p + 1; l <- (c - 1) %/% p + 1
    s2 <- 1 / (Or[k, k] * Oc[l, l])
    return(log(2 * pi * s2) + Tf[c]^2 * s2)
  }
  K <- (S - 1) %% p + 1; L <- (S - 1) %/% p + 1
  inv <- solve(Or[K, K] * Oc[L, L]); g <- inv %*% Tf[S]
  s2 <- inv[jc, jc]
  log(2 * pi * s2) + g[jc]^2 / s2
}

## ---- 3. Cell step: exact block coordinate descent over positions -------------------
# For each position (k,l): among the kept observations (v), flag those with the largest
# Delta - lambda > 0, at most n_v - h of them (forced flags in !Wfix count against the budget).
# h = ceiling(alpha * n_v) is passed as the quota for the kept observations.
cm_wstep <- function(Xc, Sr, Sc, W, Wfix, lambda, h, v = NULL) {
  p <- dim(Xc)[1]; q <- dim(Xc)[2]; n <- dim(Xc)[3]
  if (is.null(v)) v <- rep(TRUE, n); nv <- sum(v)
  Or <- solve(Sr); Oc <- solve(Sc)
  Tf <- array(0, dim(Xc)); for (i in which(v)) Tf[, , i] <- Or %*% Xc[, , i] %*% Oc
  for (c in 1:(p * q)) {
    k <- (c - 1) %% p + 1; l <- (c - 1) %/% p + 1
    free <- which(Wfix[k, l, ] & v)
    budget <- nv - h - (nv - length(free))
    W[k, l, v] <- Wfix[k, l, v]
    if (budget <= 0) next
    dec <- vapply(free, function(i) cm_delta1(Tf[, , i], Or, Oc, which(!W[, , i]), c, p), 0) - lambda[k, l]
    flag <- free[dec > 0]
    if (length(flag) > budget) flag <- flag[order(dec[dec > 0], decreasing = TRUE)[1:budget]]
    W[k, l, flag] <- FALSE
  }
  W
}

## ---- 3b. Observation step (block 2): exact minimisation over v given (Sigma, W) ------------
# Trim the observations with the largest d_i - mu > 0, at most n - h_obs of them.
# mu = -Inf gives the hard quota (classical trimmed likelihood, MMCD-type); a finite mu is the
# observation analogue of lambda (penalised trimming, cellMCD-type).
cm_vstep <- function(Xc, Sr, Sc, W, lambda, mu, h_obs) {
  n <- dim(Xc)[3]
  d <- cm_obj_i(Xc, Sr, Sc, W, lambda) - mu
  v <- rep(TRUE, n); budget <- n - h_obs
  if (budget <= 0) return(v)
  cand <- which(d > 0)
  if (length(cand) > budget) cand <- cand[order(d[cand], decreasing = TRUE)[1:budget]]
  v[cand] <- FALSE
  v
}

## ---- 4. Parameter step: fixed point of the observed-data score (one flip-flop sweep) -----
# E[X Oc X' | obs] = Xh Oc Xh' + A, A[k,k'] = sum_{a,b in M} C_ab Oc[l_a,l_b] 1[k_a=k] 1[k_b=k'],
# with Xh the conditional mean (x_S - O_SS^{-1} T_S on S) and C = O_MM^{-1}. GEM step: F does not increase.
# Only the kept observations (v) enter.
cm_pstep <- function(Xc, Sr, Sc, W, v = NULL) {
  p <- dim(Xc)[1]; q <- dim(Xc)[2]; n <- dim(Xc)[3]
  if (is.null(v)) v <- rep(TRUE, n); nv <- sum(v); keep <- which(v)
  Or <- solve(Sr); Oc <- solve(Sc)
  Xh <- Xc; info <- vector("list", n)
  for (i in keep) {
    Mi <- which(!W[, , i])
    if (length(Mi)) {
      Tm <- Or %*% Xc[, , i] %*% Oc
      K <- (Mi - 1) %% p + 1; L <- (Mi - 1) %/% p + 1
      C <- solve(Or[K, K, drop = FALSE] * Oc[L, L, drop = FALSE])
      Xh[, , i][Mi] <- Xc[, , i][Mi] - C %*% Tm[Mi]
      Ek <- matrix(0, length(Mi), p); Ek[cbind(seq_along(Mi), K)] <- 1
      El <- matrix(0, length(Mi), q); El[cbind(seq_along(Mi), L)] <- 1
      info[[i]] <- list(C = C, K = K, L = L, Ek = Ek, El = El)
    }
  }
  Sr_new <- matrix(0, p, p)
  for (i in keep) {
    Sr_new <- Sr_new + Xh[, , i] %*% Oc %*% t(Xh[, , i])
    f <- info[[i]]; if (!is.null(f)) Sr_new <- Sr_new + t(f$Ek) %*% (f$C * Oc[f$L, f$L]) %*% f$Ek
  }
  Sr_new <- Sr_new / (nv * q); Or_new <- solve(Sr_new)
  Sc_new <- matrix(0, q, q)
  for (i in keep) {
    Sc_new <- Sc_new + t(Xh[, , i]) %*% Or_new %*% Xh[, , i]
    f <- info[[i]]; if (!is.null(f)) Sc_new <- Sc_new + t(f$El) %*% (f$C * Or_new[f$K, f$K]) %*% f$El
  }
  list(Sr = Sr_new, Sc = Sc_new / (nv * p))
}

## ---- 5. Driver ---------------------------------------------------------------------
# calib = TRUE: lambda-calibration (journal, block 2). With c = chi2_{1,quant} and
# t1 = F_chi2_3(c) / F_chi2_1(c), lambda uses c / t1 + log t1 instead of c, and the final
# Sigma is divided by t1. In the 1-D proxy this makes the fixed point flag exactly 1 - quant
# of the clean cells and the returned scale Fisher-consistent; without it the fixed point
# sits at scale ~0.90 with ~1.4 % flagged (cm_proxy_fixpt).
# Observation trimming (block 2): alpha_obs < 1 allows up to n - ceiling(alpha_obs n) whole
# matrices to be trimmed; quant_obs = NULL trims exactly that many (hard quota), otherwise an
# observation is trimmed iff d_i > mu = log det Sigma_0 + pq log 2pi + c_obs with the same
# calibration (t_obs = F_{pq+2}(c) / F_pq(c), final Sigma also divided by t_obs).
cellmmcd <- function(X, alpha = 0.75, quant = 0.99, alpha_obs = 0.75, quant_obs = 0.99,
                     maxit = 100, tol = 1e-8, init = NULL, marg_cut = 3, calib = TRUE, trace = FALSE) {
  p <- dim(X)[1]; q <- dim(X)[2]; n <- dim(X)[3]
  M <- apply(X, c(1, 2), median); Xc <- sweep(X, c(1, 2), M)
  st <- if (is.null(init)) kpair(X, M) else init
  s <- if (!is.null(st$s)) st$s else apply(Xc, c(1, 2), mad)
  Wfix <- sweep(abs(Xc), c(1, 2), s, "/") <= marg_cut
  Or <- solve(st$Sr); Oc <- solve(st$Sc)
  cq <- qchisq(quant, 1); t1 <- if (calib) pchisq(cq, 3) / pchisq(cq, 1) else 1
  lambda <- cq / t1 + log(t1) + log(2 * pi) - log(outer(diag(Or), diag(Oc)))
  trim_obs <- !is.null(alpha_obs) && alpha_obs < 1
  h_obs <- if (trim_obs) ceiling(alpha_obs * n) else n
  t_obs <- 1; mu <- -Inf
  if (trim_obs && !is.null(quant_obs)) {
    co <- qchisq(quant_obs, p * q); if (calib) t_obs <- pchisq(co, p * q + 2) / pchisq(co, p * q)
    mu <- q * logdet(st$Sr) + p * logdet(st$Sc) + p * q * log(2 * pi) + co / t_obs + p * q * log(t_obs)
  }
  Sr <- st$Sr; Sc <- st$Sc; W <- Wfix; v <- rep(TRUE, n)
  obj <- numeric(0); prev <- Inf
  for (it in 1:maxit) {
    if (trim_obs) v <- cm_vstep(Xc, Sr, Sc, W, lambda, mu, h_obs)
    nv <- sum(v); h <- ceiling(alpha * nv)
    W <- cm_wstep(Xc, Sr, Sc, W, Wfix, lambda, h, v)
    f1 <- cm_obj(Xc, Sr, Sc, W, lambda, v, mu)
    ps <- cm_pstep(Xc, Sr, Sc, W, v); Sr <- ps$Sr; Sc <- ps$Sc
    f2 <- cm_obj(Xc, Sr, Sc, W, lambda, v, mu)
    obj <- c(obj, f1, f2)
    if (trace) cat(sprintf("it %d: F = %.4f (cells) %.4f (params), flagged %.3f, trimmed %d\n", it, f1, f2, mean(!W[, , v]), n - nv))
    if (abs(prev - f2) < tol * abs(f2)) break
    prev <- f2
  }
  fac <- t1 * t_obs
  sc <- scale_ident(Sr / fac, Sc, "det")
  list(Sr = sc$Sr, Sc = sc$Sc, M = M, W = W, v = v, obj = obj, iter = it, flagged = mean(!W[, , v]),
       trimmed = mean(!v), lambda = lambda, mu = mu, factor = fac, start = st, converged = it < maxit)
}

## ---- 6. Scale fixed point of the 1-D proxy (journal, block 2) ------------------------
# One cell, true conditional variance 1, working scale t, initial scale t0, threshold
# constant cc (= chi2_{1,quant}, or chi2_{1,quant} / t1 + log t1 when calibrated).
# Cell step flags iff x^2 > a(t) = t (cc + log(t0 / t)); parameter step t' = F3(a) + t (1 - F1(a)).
# Returns the fixed point t* and the flagged fraction at t*.
cm_proxy_fixpt <- function(quant = 0.99, calib = TRUE, t0 = 1, iter = 500) {
  cq <- qchisq(quant, 1); t1 <- if (calib) pchisq(cq, 3) / pchisq(cq, 1) else 1
  cc <- cq / t1 + log(t1); t <- 1
  for (k in 1:iter) { a <- t * (cc + log(t0 / t)); t <- pchisq(a, 3) + t * (1 - pchisq(a, 1)) }
  a <- t * (cc + log(t0 / t))
  c(t = t, flagged = 1 - pchisq(a, 1), t1 = t1)
}
