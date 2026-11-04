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
# mu: numeric constant per trimmed observation, a list (then cm_vstep's observed-cell threshold
# mu_i is charged, which depends on W_i), or -Inf (hard quota: nothing is charged).
cm_obj <- function(Xc, Sr, Sc, W, lambda, v = NULL, mu = 0) {
  d <- cm_obj_i(Xc, Sr, Sc, W, lambda)
  if (is.null(v) || all(v)) return(sum(d))
  if (is.list(mu)) return(sum(d[v]) + sum(cm_mu_i(Xc, W, lambda, mu)[!v]))
  sum(d[v]) + if (is.finite(mu)) mu * sum(!v) else 0
}
cm_mu_i <- function(Xc, W, lambda, mu) {
  p <- dim(Xc)[1]; q <- dim(Xc)[2]; n <- dim(Xc)[3]
  Or0 <- solve(mu$Sr0); Oc0 <- solve(mu$Sc0); base <- q * logdet(mu$Sr0) + p * logdet(mu$Sc0)
  thr <- numeric(n)
  for (i in 1:n) {
    Mi <- which(!W[, , i]); no <- p * q - length(Mi)
    ld <- if (length(Mi)) { K <- (Mi - 1) %% p + 1; L <- (Mi - 1) %/% p + 1; logdet(Or0[K, K, drop = FALSE] * Oc0[L, L, drop = FALSE]) } else 0
    thr[i] <- base + ld + no * log(2 * pi) + qchisq(mu$quant, no) / mu$t_obs + no * log(mu$t_obs) + sum(lambda[Mi])
  }
  thr
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
# shadow = TRUE (block 4): the cells of the *trimmed* observations are re-decided too, by the same rule
# without a quota (sequential coordinate descent on their own d_i). They do not enter F through d_i, but the
# observation step compares d_i - mu_i at the current W_i, so a matrix trimmed only because its cellwise
# outliers were still unflagged (W = Wfix at the first observation step) can be re-admitted once they are.
# engine = "cpp" (default when src/cm_wstep.cpp compiled, see R/load.R) runs the same loop in C++ (test 20: identical W).
cm_wstep <- function(Xc, Sr, Sc, W, Wfix, lambda, h, v = NULL, shadow = FALSE, engine = getOption("cellmmcd.engine", "R")) {
  p <- dim(Xc)[1]; q <- dim(Xc)[2]; n <- dim(Xc)[3]
  if (is.null(v)) v <- rep(TRUE, n); nv <- sum(v)
  Or <- solve(Sr); Oc <- solve(Sc)
  Tf <- array(0, dim(Xc)); for (i in if (shadow) 1:n else which(v)) Tf[, , i] <- Or %*% Xc[, , i] %*% Oc
  if (engine == "cpp") return(cm_wstep_cpp(Tf, Or, Oc, W, Wfix, lambda, as.integer(h), v, shadow))
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
  if (shadow && any(!v)) for (i in which(!v)) {
    Tfi <- Or %*% Xc[, , i] %*% Oc; Wi <- W[, , i] & Wfix[, , i]
    for (c in 1:(p * q)) {
      k <- (c - 1) %% p + 1; l <- (c - 1) %/% p + 1
      if (!Wfix[k, l, i]) next
      Wi[k, l] <- cm_delta1(Tfi, Or, Oc, which(!Wi), c, p) - lambda[k, l] <= 0
    }
    W[, , i] <- Wi
  }
  W
}

## ---- 3a. Full-column safeguard (block 3) --------------------------------------------------
# The objective is unbounded below when n - h >= n / max(p, q) (METHOD 4a): a singular Sigma^r fits one fully
# kept column and the quota allows a flag in every other column. Conjectured remedy: for every l at least p
# fully kept columns (i, l), and for every k at least q fully kept rows (i, k), among the kept observations.
# This step restores the condition after a cell step by unflagging, per violating column (row), the kept
# observations with the fewest flags in it (forced marginal flags are never unflagged; such observations are
# skipped). It can increase F; the driver records how many cells were unflagged (fit$unflagged).
cm_fullcol <- function(W, Wfix, v) {
  p <- dim(W)[1]; q <- dim(W)[2]; keep <- which(v); nun <- 0
  fix1 <- function(Wsl, Wfsl, need) {       # Wsl: cells x kept obs (logical kept); returns obs to unflag
    nf <- colSums(!Wsl); forced <- colSums(!Wfsl) > 0
    full <- sum(nf == 0); if (full >= need) return(integer(0))
    el <- which(nf > 0 & !forced); el <- el[order(nf[el])]
    el[seq_len(min(need - full, length(el)))]
  }
  for (l in 1:q) {
    j <- fix1(matrix(W[, l, keep], p), matrix(Wfix[, l, keep], p), p)
    for (i in keep[j]) { nun <- nun + sum(!W[, l, i]); W[, l, i] <- TRUE }
  }
  for (k in 1:p) {
    j <- fix1(matrix(W[k, , keep], q), matrix(Wfix[k, , keep], q), q)
    for (i in keep[j]) { nun <- nun + sum(!W[k, , i]); W[k, , i] <- TRUE }
  }
  list(W = W, unflagged = nun)
}

## ---- 3b. Observation step (block 2): exact minimisation over v given (Sigma, W) ------------
# Trim the observations with the largest d_i - mu_i > 0, at most n - h_obs of them.
# mu = -Inf: hard quota (classical trimmed likelihood, MMCD-type). Otherwise mu is a list
# (Sr0, Sc0, cobs, t_obs) and the threshold is evaluated on the observed cells of X_i only:
#   mu_i = log det Sigma0_{O_i O_i} + |O_i| log 2pi + c(|O_i|) / t_obs + |O_i| log t_obs + sum lambda[M_i],
# with c(d) = chi2_{d, quant_obs} and log det Sigma0_{OO} by the Schur identity at the initial estimate.
# Flagged cells therefore do not count against their matrix; a matrix is trimmed when its own kept
# cells are jointly beyond their quant_obs quantile (at the initial scale, like lambda).
# max_flag: a kept observation may have at most this fraction of its cells flagged; matrices beyond
# it are trimmed first (constraint on (v, W), see METHOD.md 4c).
# count_max (block 4): a kept observation may have at most this many flagged cells (adaptive majority rule, driver:
# the flag_quant quantile of Binomial(pq, eps_bar) at the current overall flag rate eps_bar); Inf = off.
cm_vstep <- function(Xc, Sr, Sc, W, lambda, mu, h_obs, max_flag = 1, count_max = Inf) {
  p <- dim(Xc)[1]; q <- dim(Xc)[2]; n <- dim(Xc)[3]
  d <- cm_obj_i(Xc, Sr, Sc, W, lambda)
  v <- rep(TRUE, n); budget <- n - h_obs
  if (budget <= 0) return(v)
  if (is.list(mu)) { d <- d - cm_mu_i(Xc, W, lambda, mu); cand <- which(d > 0) } else cand <- seq_len(n)
  nfl <- apply(!W, 3, sum); forced <- which(nfl > max_flag * p * q | nfl > count_max)
  if (length(forced)) { d[forced] <- Inf; cand <- union(forced, cand) }
  if (length(cand) > budget) cand <- cand[order(d[cand], decreasing = TRUE)[1:budget]]
  v[cand] <- FALSE
  v
}

## ---- 3c. Robust flag rate for the count rule (block 5) ------------------------------------
# The count rule (cm_vstep count_max) needs the cellwise flag rate eps_bar of the clean matrices. The
# overall rate among kept matrices is inflated by kept casewise outliers (scaled / structural matrices
# flag 15-40 % of their cells), which raises count_max until they pass (block 4: struct 20 % trimmed
# only 60 %). Robust version: mean flag count over the alpha fraction of kept matrices with the fewest
# flags, corrected for the trimming by the binomial model: under independent cellwise contamination the
# counts are Binomial(pq, eps), and g(eps) = E[mean of the lowest alpha fraction of Bin(pq, eps)] / pq is
# increasing, so eps_bar solves g(eps_bar) = observed lower-alpha mean / pq (uniroot). Fisher-consistent
# for eps under the binomial model (checked, test 21), and unaffected by up to (1 - alpha) matrices.
cm_eps_robust <- function(W, v, alpha = 0.75) {
  pq <- prod(dim(W)[1:2]); cnt <- apply(!W[, , v, drop = FALSE], 3, sum)
  m <- max(1, floor(alpha * length(cnt))); tm <- mean(sort(cnt)[1:m]) / pq
  if (tm <= 0) return(0)
  g <- function(e) { k <- qbinom(alpha, pq, e); ks <- seq_len(k) - 1
    (sum(ks * dbinom(ks, pq, e)) + k * (alpha - pbinom(k - 1, pq, e))) / alpha / pq }
  if (g(0.5) <= tm) return(0.5)
  uniroot(function(e) g(e) - tm, c(0, 0.5), tol = 1e-6)$root
}

## ---- 4. Parameter step: fixed point of the observed-data score (one flip-flop sweep) -----
# E[X Oc X' | obs] = Xh Oc Xh' + A, A[k,k'] = sum_{a,b in M} C_ab Oc[l_a,l_b] 1[k_a=k] 1[k_b=k'],
# with Xh the conditional mean (x_S - O_SS^{-1} T_S on S) and C = O_MM^{-1}. GEM step: F does not increase.
# Only the kept observations (v) enter.
# tau (block 3): multiplies the conditional covariance C of the flagged cells in the second-moment
# terms (not in the conditional mean). tau = 1 is EM for missing cells; tau = (1 - F3(c)) / (1 - F1(c))
# is the second moment of a N(0,1) residual beyond its cutoff c, i.e. flagged cells are treated as
# censored (known to be in the tail) rather than missing; in the 1-D proxy this makes the fixed point
# Fisher-consistent without the t1 factor (cm_proxy_fixpt(tau = ...)).
cm_pstep <- function(Xc, Sr, Sc, W, v = NULL, tau = 1) {
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
      info[[i]] <- list(C = tau * C, K = K, L = L, Ek = Ek, El = El)
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
# observation is trimmed iff d_i > mu_i, the quant_obs quantile of its own kept cells at the
# (hard quota, block 5: the scale is multiplied by the MCD consistency factor kappa_mcd(alpha_obs, pq), d = pq as above)
# initial estimate (cm_vstep), with the same calibration (t_obs = F_{pq+2}(c) / F_pq(c) at d = pq,
# an approximation for matrices with flagged cells; final Sigma also divided by t_obs).
# max_flag: matrices with more than this fraction of flagged cells are trimmed as a whole (they are
# casewise outliers by any standard; a shifted matrix keeps a few unflagged cells that look fine alone).
# tail = TRUE (block 3): censored treatment of flagged cells in the parameter step (cm_pstep tau), with
# calib switched off (the proxy fixed point is then at scale 1 with 1 - quant flagged).
# safeguard = TRUE (block 3): full-column/row constraint restored after each cell step (cm_fullcol, 3a).
# reweight = TRUE (block 3, the default since the bias finding of exp/cm_step6_bias.R): after convergence,
# re-flag at the quant_rw quantile and refit with the flags fixed (cm_reweight), the analogue of the
# reweighted MCD; the raw fit is returned as $raw. reweight = FALSE gives the raw estimator, which is
# asymptotically biased in shape at the Gaussian (JOURNAL, block 3).
# flag_quant (block 4): adaptive majority rule. With the shadow step, a casewise outlier whose worst cells are flagged
# can pass the kept-cell test (selection); so a kept matrix may have at most qbinom(flag_quant, pq, eps_bar) flagged
# cells, eps_bar = current overall flag rate among the kept matrices (at least 1 - quant). NULL = off (max_flag only).
# eps_robust (block 5): eps_bar from the alpha_obs fraction of kept matrices with the fewest flags, binomial-corrected
# (cm_eps_robust, 3c), instead of the plain mean, which kept casewise outliers inflate. FALSE = block-4 rule.
# shadow_it (block 4): during the first shadow_it iterations the cell step also re-decides the cells of the
# trimmed observations (cm_wstep shadow), so that the observation step can re-admit a matrix whose cellwise
# outliers are flagged by then. Those cells change F only through the trimmed charge mu_i(W_i), so F need
# not be monotone during these iterations; afterwards the algorithm is the block-3 one (monotone), which
# makes the convergence argument hold after at most shadow_it non-monotone iterations. 0 = off.
cellmmcd <- function(X, alpha = 0.75, quant = 0.99, alpha_obs = 0.75, quant_obs = 0.99, max_flag = 0.5,
                     maxit = 100, tol = 1e-8, init = NULL, marg_cut = 3, calib = TRUE, tail = FALSE,
                     reweight = TRUE, quant_rw = 0.999, safeguard = FALSE, shadow_it = 10, flag_quant = 0.999,
                     eps_robust = TRUE, trace = FALSE) {
  p <- dim(X)[1]; q <- dim(X)[2]; n <- dim(X)[3]
  M <- apply(X, c(1, 2), median); Xc <- sweep(X, c(1, 2), M)
  st <- if (is.null(init)) kpair(X, M) else init
  s <- if (!is.null(st$s)) st$s else apply(Xc, c(1, 2), mad)
  Wfix <- sweep(abs(Xc), c(1, 2), s, "/") <= marg_cut
  Or <- solve(st$Sr); Oc <- solve(st$Sc)
  cq <- qchisq(quant, 1); t1 <- if (calib && !tail) pchisq(cq, 3) / pchisq(cq, 1) else 1
  tau <- if (tail) (1 - pchisq(cq, 3)) / (1 - pchisq(cq, 1)) else 1
  lambda <- cq / t1 + log(t1) + log(2 * pi) - log(outer(diag(Or), diag(Oc)))
  trim_obs <- !is.null(alpha_obs) && alpha_obs < 1
  h_obs <- if (trim_obs) ceiling(alpha_obs * n) else n
  t_obs <- 1; mu <- -Inf
  if (trim_obs && is.null(quant_obs) && calib) t_obs <- 1 / kappa_mcd(alpha_obs, p * q)   # hard quota (block 5): MCD consistency factor at d = pq
  if (trim_obs && !is.null(quant_obs)) {
    co <- qchisq(quant_obs, p * q); if (calib) t_obs <- pchisq(co, p * q + 2) / pchisq(co, p * q)
    mu <- list(Sr0 = st$Sr, Sc0 = st$Sc, quant = quant_obs, t_obs = t_obs)
  }
  Sr <- st$Sr; Sc <- st$Sc; W <- Wfix; v <- rep(TRUE, n)
  obj <- numeric(0); prev <- Inf; unflagged <- 0; eps_bar <- NA; count_max <- Inf
  for (it in 1:maxit) {
    if (trim_obs) {
      eps_bar <- max(if (eps_robust) cm_eps_robust(W, v, alpha_obs) else mean(!W[, , v]), 1 - quant)
      count_max <- if (is.null(flag_quant)) Inf else qbinom(flag_quant, p * q, eps_bar)
      v <- cm_vstep(Xc, Sr, Sc, W, lambda, mu, h_obs, max_flag, count_max)
    }
    nv <- sum(v); h <- ceiling(alpha * nv)
    W <- cm_wstep(Xc, Sr, Sc, W, Wfix, lambda, h, v, shadow = trim_obs && it <= shadow_it)
    if (safeguard) { sg <- cm_fullcol(W, Wfix, v); W <- sg$W; unflagged <- unflagged + sg$unflagged }
    f1 <- cm_obj(Xc, Sr, Sc, W, lambda, v, mu)
    ps <- cm_pstep(Xc, Sr, Sc, W, v, tau); Sr <- ps$Sr; Sc <- ps$Sc
    f2 <- cm_obj(Xc, Sr, Sc, W, lambda, v, mu)
    obj <- c(obj, f1, f2)
    if (trace) cat(sprintf("it %d: F = %.4f (cells) %.4f (params), flagged %.3f, trimmed %d\n", it, f1, f2, mean(!W[, , v]), n - nv))
    if (abs(prev - f2) < tol * abs(f2)) break
    prev <- f2
  }
  fac <- t1 * t_obs
  sc <- scale_ident(Sr / fac, Sc, "det")
  fit <- list(Sr = sc$Sr, Sc = sc$Sc, M = M, W = W, Wfix = Wfix, v = v, obj = obj, iter = it,
              flagged = mean(!W[, , v]), trimmed = mean(!v), lambda = lambda, mu = mu, factor = fac,
              start = st, converged = it < maxit, tau = tau, unflagged = unflagged, shadow_it = shadow_it,
              eps_bar = eps_bar, count_max = count_max, t_obs = t_obs)
  if (reweight) cm_reweight(fit, X, quant_rw = quant_rw, maxit = maxit, tol = tol) else fit
}

## ---- 5b. Reweighting step (block 3) --------------------------------------------------------
# Squared standardized conditional residual of cell c given the kept cells of its matrix other than c.
cm_r2 <- function(Tf, Or, Oc, Mi, c, p) {
  S <- c(setdiff(Mi, c), c); jc <- length(S)
  K <- (S - 1) %% p + 1; L <- (S - 1) %/% p + 1
  inv <- solve(Or[K, K, drop = FALSE] * Oc[L, L, drop = FALSE]); g <- inv %*% Tf[S]
  as.numeric(g[jc]^2 / inv[jc, jc])
}
# Re-decide every cell of every kept matrix at the raw fit: flag iff r^2 > chi2_{1, quant_rw} (given the
# other kept cells of the raw fit) or forced by the marginal cut; keep the observation weights v of the raw
# fit; then run parameter steps with (v, W) fixed to the observed-cell MLE (EM, monotone); divide the scale
# by t_rw = F3(c_rw) / F1(c_rw) (proxy calibration of the ratio estimator at cutoff c_rw) and by the raw
# fit's t_obs. No quota is enforced in this step (as in the reweighted MCD).
cm_reweight <- function(fit, X, quant_rw = 0.999, maxit = 100, tol = 1e-8) {
  p <- dim(X)[1]; q <- dim(X)[2]; n <- dim(X)[3]
  Xc <- sweep(X, c(1, 2), fit$M); Or <- solve(fit$Sr); Oc <- solve(fit$Sc)
  crw <- qchisq(quant_rw, 1); v <- fit$v; W <- fit$W
  if (getOption("cellmmcd.engine", "R") == "cpp") {
    Tf <- array(0, dim(Xc)); for (i in which(v)) Tf[, , i] <- Or %*% Xc[, , i] %*% Oc
    W <- cm_reweight_cpp(Tf, Or, Oc, fit$W, fit$Wfix, crw, v)
  } else for (i in which(v)) {
    Tf <- Or %*% Xc[, , i] %*% Oc; Mi <- which(!fit$W[, , i])
    r2 <- vapply(1:(p * q), function(c) cm_r2(Tf, Or, Oc, Mi, c, p), 0)
    W[, , i] <- fit$Wfix[, , i] & (r2 <= crw)
  }
  Sr <- fit$Sr; Sc <- fit$Sc; zero <- matrix(0, p, q); prev <- Inf; obj <- numeric(0)
  for (it in 1:maxit) {
    ps <- cm_pstep(Xc, Sr, Sc, W, v); Sr <- ps$Sr; Sc <- ps$Sc
    f2 <- cm_obj(Xc, Sr, Sc, W, zero, v, 0); obj <- c(obj, f2)
    if (abs(prev - f2) < tol * abs(f2)) break
    prev <- f2
  }
  t_rw <- pchisq(crw, 3) / pchisq(crw, 1)
  t_obs <- if (is.list(fit$mu)) fit$mu$t_obs else fit$t_obs
  sc <- scale_ident(Sr / (t_rw * t_obs), Sc, "det")
  list(Sr = sc$Sr, Sc = sc$Sc, M = fit$M, W = W, Wfix = fit$Wfix, v = v, obj = obj, iter = it,
       flagged = mean(!W[, , v]), trimmed = mean(!v), factor = t_rw * t_obs, converged = it < maxit,
       raw = fit, quant_rw = quant_rw, lambda = fit$lambda, mu = fit$mu, start = fit$start, tau = fit$tau,
       unflagged = fit$unflagged)
}

## ---- 6. Scale fixed point of the 1-D proxy (journal, block 2) ------------------------
# One cell, true conditional variance 1, working scale t, initial scale t0, threshold
# constant cc (= chi2_{1,quant}, or chi2_{1,quant} / t1 + log t1 when calibrated).
# Cell step flags iff x^2 > a(t) = t (cc + log(t0 / t)); parameter step t' = F3(a) + t (1 - F1(a)).
# Returns the fixed point t* and the flagged fraction at t*.
cm_proxy_fixpt <- function(quant = 0.99, calib = TRUE, t0 = 1, iter = 500, tau = 1) {
  cq <- qchisq(quant, 1); t1 <- if (calib) pchisq(cq, 3) / pchisq(cq, 1) else 1
  cc <- cq / t1 + log(t1); t <- 1
  for (k in 1:iter) { a <- t * (cc + log(t0 / t)); t <- pchisq(a, 3) + tau * t * (1 - pchisq(a, 1)) }
  a <- t * (cc + log(t0 / t))
  c(t = t, flagged = 1 - pchisq(a, 1), t1 = t1)
}

## ---- 7. Influence functions in the 1-D proxy (block 4) -----------------------------------
# One cell with residual distribution F; the start is the squared standardized MAD t0(F); the raw
# estimator is the calibrated fixed point T_raw(F) = t*/t1 of t <- N_F(a(t)) + t (1 - D_F(a(t))),
# a(t) = t (c/t1 + log t1 + log(t0/t)), with N_F(a) = int x^2 1(x^2 <= a) dF and D_F(a) = int 1(x^2 <= a) dF
# (flagged cells are missing: they contribute the current scale t); the reweighted estimator is
# T_rw(F) = N_F(a_rw) / D_F(a_rw) / t_rw with a_rw = c_rw T_raw(F). Both equal 1 at F = Phi (Fisher-consistent).
# cm_proxy_T evaluates them exactly under F_eps = (1 - eps) Phi + eps delta_x (integrals in closed form);
# cm_proxy_if gives the influence functions derived in JOURNAL.md (block 4); test 19 compares them.
cm_proxy_T <- function(x, eps = 0, quant = 0.99, quant_rw = 0.999, iter = 2000) {
  c1 <- qchisq(quant, 1); t1 <- pchisq(c1, 3) / pchisq(c1, 1); cc <- c1 / t1 + log(t1)
  crw <- qchisq(quant_rw, 1); trw <- pchisq(crw, 3) / pchisq(crw, 1)
  N <- function(a) (1 - eps) * pchisq(a, 3) + eps * x^2 * (x^2 <= a)
  D <- function(a) (1 - eps) * pchisq(a, 1) + eps * (x^2 <= a)
  q3 <- qnorm(0.75); m <- qnorm((1 + 0.5 / (1 - eps)) / 2)          # median of |X| under F_eps if |x| > m
  if (abs(x) <= m) m <- qnorm((1 + (0.5 - eps) / (1 - eps)) / 2)
  t0 <- (m / q3)^2; t <- t0
  for (k in 1:iter) { a <- t * (cc + log(t0 / t)); t <- N(a) + t * (1 - D(a)) }
  raw <- t / t1; arw <- crw * raw
  c(raw = raw, rw = N(arw) / D(arw) / trw, t0 = t0)
}
cm_proxy_if <- function(x, quant = 0.99, quant_rw = 0.999) {
  c1 <- qchisq(quant, 1); t1 <- pchisq(c1, 3) / pchisq(c1, 1); s1 <- sqrt(c1); ph1 <- dnorm(s1)
  crw <- qchisq(quant_rw, 1); F1 <- pchisq(crw, 1); F3 <- pchisq(crw, 3); srw <- sqrt(crw); phrw <- dnorm(srw)
  q3 <- qnorm(0.75); if_t0 <- 2 * sign(abs(x) - q3) / (4 * q3 * dnorm(q3))   # IF of the squared standardized MAD
  num <- (x^2 <= c1) * (x^2 - t1) + t1 * ph1 * (c1 - t1) / s1 * if_t0
  den <- (c1 - t1)^2 * ph1 / (t1 * s1) - pchisq(c1, 1)
  if_raw <- -num / den / t1
  kept <- x^2 <= crw
  if_rw <- kept * (x^2 / F3 - 1 / F1) + srw * phrw * (crw * F1 - F3) / (F1 * F3) * if_raw
  c(raw = if_raw, rw = if_rw, t0 = if_t0)
}
