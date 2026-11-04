# cellMMCD step 12 (block 5): why are structural casewise outliers at 20 % trimmed only 60 %? Seed 12.
# Setting of cm_step4_mixed (p = q = 5, n = 100, struct = Kronecker model with flipped correlation signs).
# For each matrix at the final raw fit: number of flagged cells, d_i - mu_i (kept-cell excess), trimmed or not;
# eps_bar and count_max used by the count rule, for the block-4 rule (overall mean) and the block-5 robust rate.
# Output: output/cm_step12_struct.log
source("R/load.R"); library(parallel); set.seed(12)
source_lines <- readLines("exp/cm_step4_mixed.R"); eval(parse(text = source_lines[grep("^contaminate <- ", source_lines):(grep("^one_rep <- ", source_lines) - 1)]))
diag_rep <- function(seed, type, eps_case, eps_cell, eps_robust) {
  set.seed(seed); p <- q <- 5; n <- 100; S <- true_sigma(p, q, "full")
  d <- contaminate(rmatnorm(n, Sr = S$Sr, Sc = S$Sc), S, type, eps_case, eps_cell); X <- d$X
  f <- cellmmcd(X, eps_robust = eps_robust); r <- f$raw
  Xc <- sweep(X, c(1, 2), r$M); Sr <- r$Sr * r$factor
  dd <- cm_obj_i(Xc, Sr, r$Sc, r$W, r$lambda) - cm_mu_i(Xc, r$W, r$lambda, r$mu)
  cnt <- apply(!r$W, 3, sum)
  list(kl = kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), bad = d$bad_case, cnt = cnt, dd = dd, v = r$v, eps_bar = r$eps_bar, count_max = r$count_max,
       trim_bad = if (any(d$bad_case)) mean(!r$v[d$bad_case]) else NA, trim_good = mean(!r$v[!d$bad_case]))
}
grid <- rbind(data.frame(type = "struct", eps_case = c(0.1, 0.2), eps_cell = 0), data.frame(type = "scale", eps_case = 0.2, eps_cell = 0),
              data.frame(type = "none", eps_case = 0, eps_cell = c(0, 0.05, 0.1)), data.frame(type = "struct", eps_case = 0.2, eps_cell = 0.05))
sink("output/cm_step12_struct.log"); cat("cellMMCD step 12: p = q = 5, n = 100, 10 reps; final raw fit. cnt = flagged cells per matrix, dd = d_i - mu_i.\n")
for (g in seq_len(nrow(grid))) for (er in c(FALSE, TRUE)) {
  reps <- mclapply(1:10, function(r) diag_rep(12000 + 100 * g + r, grid$type[g], grid$eps_case[g], grid$eps_cell[g], er), mc.cores = 4)
  bad <- unlist(lapply(reps, `[[`, "bad")); cnt <- unlist(lapply(reps, `[[`, "cnt")); dd <- unlist(lapply(reps, `[[`, "dd")); v <- unlist(lapply(reps, `[[`, "v"))
  cat(sprintf("\n== %s case %.2f cells %.2f, eps_robust = %s: median KL %.3f, trimmed bad %.2f good %.3f, eps_bar %.3f (range %.3f-%.3f), count_max %s\n",
              grid$type[g], grid$eps_case[g], grid$eps_cell[g], er, median(sapply(reps, `[[`, "kl")), median(sapply(reps, `[[`, "trim_bad")), median(sapply(reps, `[[`, "trim_good")),
              median(sapply(reps, `[[`, "eps_bar")), min(sapply(reps, `[[`, "eps_bar")), max(sapply(reps, `[[`, "eps_bar")), paste(sort(unique(sapply(reps, `[[`, "count_max"))), collapse = "/")))
  for (b in unique(bad)) { sel <- bad == b
    cat(sprintf("  %s matrices (%d): flagged count quantiles (0,.1,.25,.5,.75,.9,1): %s; dd quantiles (.1,.5,.9): %s; trimmed %.2f; trimmed among cnt<=count_max %.2f\n",
                if (b) "bad " else "good", sum(sel), paste(round(quantile(cnt[sel], c(0, .1, .25, .5, .75, .9, 1)), 1), collapse = " "),
                paste(round(quantile(dd[sel], c(.1, .5, .9)), 1), collapse = " "), mean(!v[sel]), mean(!v[sel & cnt <= median(sapply(reps, `[[`, "count_max"))]))) }
}
sink(); cat("done\n")
