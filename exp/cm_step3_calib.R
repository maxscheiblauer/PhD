# cellMMCD step 3: lambda-calibration (journal, block 2). Clean data, calib on/off, vs the MLE.
# Seed 3. Output: output/cm_step3_calib.log
source("R/load.R"); library(parallel); set.seed(3)
sink("output/cm_step3_calib.log")
fx <- cm_proxy_fixpt(calib = FALSE); cat(sprintf("proxy fixed point, uncalibrated: t* = %.4f, flagged %.4f; calibrated: t1 = %.4f, flagged .0100\n", fx["t"], fx["flagged"], cm_proxy_fixpt()["t1"]))
for (pq in list(c(3, 4), c(5, 5))) for (n in c(100, 1000)) {
  p <- pq[1]; q <- pq[2]; S <- true_sigma(p, q, "full")
  r <- mclapply(1:8, function(k) { set.seed(100 * n + k); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
    f0 <- cellmmcd(X, calib = FALSE); f1 <- cellmmcd(X, calib = TRUE); ff <- flipflop(X)
    c(tr0 = tr_ratio(f0$Sr, f0$Sc, S$Sr, S$Sc), fl0 = f0$flagged, kl0 = kl_mn(f0$Sr, f0$Sc, S$Sr, S$Sc),
      tr1 = tr_ratio(f1$Sr, f1$Sc, S$Sr, S$Sc), fl1 = f1$flagged, kl1 = kl_mn(f1$Sr, f1$Sc, S$Sr, S$Sc), skl1 = kl_shape(f1$Sr, f1$Sc, S$Sr, S$Sc),
      klff = kl_mn(ff$Sr, ff$Sc, S$Sr, S$Sc), it1 = f1$iter) }, mc.cores = 4)
  v <- colMeans(do.call(rbind, r))
  cat(sprintf("%dx%d n=%4d | uncalibrated: tr %.3f flagged %.3f KL %.3f | calibrated: tr %.3f flagged %.3f KL %.3f (shape %.3f) | MLE KL %.3f | iter %.0f\n",
              p, q, n, v["tr0"], v["fl0"], v["kl0"], v["tr1"], v["fl1"], v["kl1"], v["skl1"], v["klff"], v["it1"]))
}
sink()
