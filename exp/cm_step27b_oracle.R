# Block 11, step 27b: what drives the KL rise at eps = .2 / .25 (20x20, n = 100, same seeds as step 27 rows 13-18)?
# Oracle: the observed-cell MLE with the TRUE cell flags (cm_pstep iterated from kpair, no trimming), i.e. the pure
# information loss of discarding the contaminated cells; plus precision / recall of the default fit's flags and the
# number of fully kept columns of the true pattern. Output: output/cm_step27b_oracle.log
source("R/load.R"); p <- q <- 20; n <- 100
sink("output/cm_step27b_oracle.log"); cat("step 27b (block 11): oracle observed-cell MLE with true flags vs default cellMMCD at 20x20, n = 100\n")
for (k in 13:18) { eps <- c(0.2, 0.25)[(k - 13) %/% 3 + 1]; set.seed(2700 + k); S <- true_sigma(p, q, "full"); X0 <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); X <- X0
  idx <- sample(p * q * n, round(eps * p * q * n)); ai <- arrayInd(idx, dim(X)); X[idx] <- X[idx] + 6 * sample(c(-1, 1), length(idx), TRUE) * sqrt(diag(S$Sr)[ai[, 1]] * diag(S$Sc)[ai[, 2]])
  Wt <- array(TRUE, dim(X)); Wt[idx] <- FALSE; Xc <- sweep(X, c(1, 2), apply(X, c(1, 2), median)); st <- kpair(X); Sr <- st$Sr; Sc <- st$Sc
  for (it in 1:60) { ps <- cm_pstep(Xc, Sr, Sc, Wt, NULL, 1, TRUE); Sr <- ps$Sr; Sc <- ps$Sc; Xc <- sweep(Xc, c(1, 2), ps$dM) }
  f <- cellmmcd(X); A <- cm_check_pattern(Wt); pr <- prec_recall(f$W, !Wt)
  cat(sprintf("  eps %.2f seed %d: oracle KL %.3f (true pattern: c_full %d, d_full %d); cellMMCD KL %.3f, flag precision %.3f recall %.3f, trimmed %.2f, kept-pattern c_full %d\n",
    eps, 2700 + k, kl_mn(Sr, Sc, S$Sr, S$Sc), A$c_full, A$d_full, kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), pr["precision"], pr["recall"], f$trimmed, cm_check_pattern(f$raw)$c_full)) }
sink(); cat(readLines("output/cm_step27b_oracle.log"), sep = "\n")
