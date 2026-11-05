# Block 11, step 27: the large-p,q operating limit (STATE next-step 1). At 20x20, n = 100, random cells at rate eps
# (6 sd shifts as in step 26 (e)), the kept pattern falls into Prop. A's unbounded regime once n q (1-eps)^p < p.
# (a) Quantify: KL of the default cellMMCD vs the clean-data MLE (oracle flip-flop on the uncontaminated X) and vs MMCD,
# eps in {0, .05, .1, .15, .2, .25}, 3 reps; fully kept columns / rows of the kept pattern, Prop. A, condition numbers.
# (b) Remedy (i): the ridge parameter step (cm_ridge, rho = 0.05 and 0.2); bias at the Gaussian read off eps = 0.
# Seed 27. Output: output/cm_step27_limit.log, output/cm_step27_limit.rds
source("R/load.R"); set.seed(27); library(parallel)
p <- q <- 20; n <- 100; eps_grid <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25); reps <- 3
grid <- expand.grid(rep = 1:reps, eps = eps_grid); grid$seed <- 2700 + seq_len(nrow(grid))
one <- function(g) { set.seed(g$seed); S <- true_sigma(p, q, "full"); X0 <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); X <- X0
  if (g$eps > 0) { idx <- sample(p * q * n, round(g$eps * p * q * n)); ai <- arrayInd(idx, dim(X))
    X[idx] <- X[idx] + 6 * sample(c(-1, 1), length(idx), TRUE) * sqrt(diag(S$Sr)[ai[, 1]] * diag(S$Sc)[ai[, 2]]) }
  out <- list(eps = g$eps, rep = g$rep, exp_full = n * q * (1 - g$eps)^p)
  ff <- flipflop(X0, center = TRUE); out$kl_mle <- kl_mn(ff$Sr, ff$Sc, S$Sr, S$Sc)
  if (has_robustmatrix) { m <- tryCatch(robustmatrix::mmcd(X, alpha = 0.75), error = function(e) NULL)
    out$kl_mmcd <- if (is.null(m)) NA else kl_mn(m$cov_row, m$cov_col, S$Sr, S$Sc) } else out$kl_mmcd <- NA
  for (r in c(0, 0.05, 0.2)) { t0 <- proc.time()[3]; f <- tryCatch(cellmmcd(X, ridge = r), error = function(e) NULL); nm <- paste0("r", r)
    if (is.null(f)) { out[[paste0("kl_", nm)]] <- NA; next }
    A <- cm_check_pattern(f$raw); out[[paste0("kl_", nm)]] <- kl_mn(f$Sr, f$Sc, S$Sr, S$Sc)
    out[[paste0("cf_", nm)]] <- A$c_full; out[[paste0("df_", nm)]] <- A$d_full
    out[[paste0("propA_", nm)]] <- as.numeric((A$c_full >= 1 && A$c_full <= p - 1) || (A$d_full >= 1 && A$d_full <= q - 1))
    out[[paste0("kr_", nm)]] <- kappa(f$Sr, exact = TRUE); out[[paste0("kc_", nm)]] <- kappa(f$Sc, exact = TRUE)
    out[[paste0("flag_", nm)]] <- f$flagged; out[[paste0("trim_", nm)]] <- f$trimmed; out[[paste0("sec_", nm)]] <- proc.time()[3] - t0 }
  out$kr_true <- kappa(S$Sr, exact = TRUE); out$kc_true <- kappa(S$Sc, exact = TRUE); out }
res <- mclapply(seq_len(nrow(grid)), function(k) one(grid[k, ]), mc.cores = 4)
saveRDS(res, "output/cm_step27_limit.rds")
D <- do.call(rbind, lapply(res, function(o) as.data.frame(o[!vapply(o, is.null, TRUE)])))
sink("output/cm_step27_limit.log"); cat("cellMMCD step 27 (block 11): the large-p,q operating limit at 20x20, n = 100 (seed 27, 3 reps, 6 sd cells)\n")
cat(sprintf("true condition numbers kappa(Sr) %.1f kappa(Sc) %.1f\n", D$kr_true[1], D$kc_true[1]))
cat("\n(a)+(b) median over reps; exp_full = n q (1-eps)^p (Prop. A unbounded regime once < p = 20); r0 = default, r.05 / r.2 = ridge rho\n")
cat(sprintf("%5s %8s | %7s %7s | %7s %7s %7s | %5s %5s %5s | %5s %5s %5s | %6s %6s %6s | %5s %5s %5s\n", "eps", "exp_full", "KL_mle", "KL_mmcd", "KL_r0", "KL_r.05", "KL_r.2", "cf_r0", "cf.05", "cf.2", "A_r0", "A.05", "A.2", "kr_r0", "kr.05", "kr.2", "fl_r0", "fl.05", "fl.2"))
for (e in eps_grid) { d <- D[D$eps == e, ]; m <- function(x) median(x, na.rm = TRUE)
  cat(sprintf("%5.2f %8.1f | %7.3f %7.3f | %7.3f %7.3f %7.3f | %5.0f %5.0f %5.0f | %5.0f %5.0f %5.0f | %6.0f %6.0f %6.0f | %5.3f %5.3f %5.3f\n", e, d$exp_full[1], m(d$kl_mle), m(d$kl_mmcd),
    m(d$kl_r0), m(d$kl_r0.05), m(d$kl_r0.2), m(d$cf_r0), m(d$cf_r0.05), m(d$cf_r0.2), sum(d$propA_r0), sum(d$propA_r0.05), sum(d$propA_r0.2),
    m(d$kr_r0), m(d$kr_r0.05), m(d$kr_r0.2), m(d$flag_r0), m(d$flag_r0.05), m(d$flag_r0.2))) }
cat("\nper-rep KL (r0 / r.05 / r.2), fully kept cols/rows (r0), trimmed fraction (r0), seconds (r0):\n")
for (k in seq_len(nrow(D))) with(D[k, ], cat(sprintf("  eps %.2f rep %d: KL %.3f / %.3f / %.3f; cf %d df %d; trimmed %.2f; %.0f s\n", eps, rep, kl_r0, kl_r0.05, kl_r0.2, cf_r0, df_r0, trim_r0, sec_r0)))
sink(); cat(readLines("output/cm_step27_limit.log"), sep = "\n")
