# Diagnostic for the efficiency gap (block 2): is it the noise of lambda from the kpair start?
# Clean 3x4 and 5x5, n = 1000, 8 seeds: cellMMCD (calibrated) with (a) default start, (b) lambda and
# start from the true Sigma (oracle), (c) default start but lambda recomputed from the final fit and
# one more run from that fit. Output: output/cm_step3b_lambda.log
source("R/load.R"); library(parallel); set.seed(33)
sink("output/cm_step3b_lambda.log")
for (pq in list(c(3, 4), c(5, 5))) { p <- pq[1]; q <- pq[2]; S <- true_sigma(p, q, "full"); n <- 1000
  r <- mclapply(1:8, function(k) { set.seed(3000 + k); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc)
    fa <- cellmmcd(X); fb <- cellmmcd(X, init = list(Sr = S$Sr, Sc = S$Sc))
    fc <- cellmmcd(X, init = list(Sr = fa$Sr * fa$factor, Sc = fa$Sc)); ff <- flipflop(X)
    c(a = kl_mn(fa$Sr, fa$Sc, S$Sr, S$Sc), b = kl_mn(fb$Sr, fb$Sc, S$Sr, S$Sc), c = kl_mn(fc$Sr, fc$Sc, S$Sr, S$Sc), mle = kl_mn(ff$Sr, ff$Sc, S$Sr, S$Sc),
      fla = fa$flagged, flb = fb$flagged, flc = fc$flagged, tra = tr_ratio(fa$Sr, fa$Sc, S$Sr, S$Sc), trb = tr_ratio(fb$Sr, fb$Sc, S$Sr, S$Sc)) }, mc.cores = 4)
  v <- colMeans(do.call(rbind, r))
  cat(sprintf("%dx%d n=1000 mean KL: default %.4f (flagged %.3f, tr %.3f) | oracle start+lambda %.4f (flagged %.3f, tr %.3f) | restart from fit %.4f (flagged %.3f) | MLE %.4f\n",
              p, q, v["a"], v["fla"], v["tra"], v["b"], v["flb"], v["trb"], v["c"], v["flc"], v["mle"]))
}
sink()
