# cellMMCD step 2: is the clean-data scale loss the 1-D tail-trimming factor? Seed 2.
# Prediction: kept-cell variance ratio F_chi2_3(c) / F_chi2_1(c) with c = chi2_{1,.99} (= 0.924).
source("R/load.R"); library(parallel); set.seed(2)
c99 <- qchisq(0.99, 1); pred <- pchisq(c99, 3) / pchisq(c99, 1)
sink("output/cm_step2_factor.log")
cat(sprintf("predicted scale ratio at quant .99: %.4f\n", pred))
for (pq in list(c(3, 4), c(5, 5))) for (n in c(100, 1000)) {
  p <- pq[1]; q <- pq[2]; S <- true_sigma(p, q, "full")
  r <- mclapply(1:8, function(k) { set.seed(100 * n + k); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); f <- cellmmcd(X)
    c(tr = tr_ratio(f$Sr, f$Sc, S$Sr, S$Sc), fl = f$flagged, skl = kl_shape(f$Sr, f$Sc, S$Sr, S$Sc), klc = kl_mn(f$Sr / pred, f$Sc, S$Sr, S$Sc)) }, mc.cores = 4)
  v <- do.call(rbind, r)
  cat(sprintf("%dx%d n=%4d: tr ratio %.3f (se %.3f), flagged %.3f, shape KL %.3f, KL after dividing by %.3f: %.3f\n",
              p, q, n, mean(v[, "tr"]), sd(v[, "tr"]) / sqrt(8), mean(v[, "fl"]), mean(v[, "skl"]), pred, mean(v[, "klc"])))
}
sink()
