# Step 2 (Q2): do denominator, kappa and scale identification only change the scale?
# Run from project root: Rscript exp/step2_normalization.R
source("R/load.R")
set.seed(202)
n <- 25
rules <- c("R0", "R1", "R3", "R2", "RJ")
grid <- expand.grid(denom = c("h", "n"), kappa = c(TRUE, FALSE),
                    scale_id = c("det", "none", "tr", "c11"), stringsAsFactors = FALSE)
ref_id <- which(grid$denom == "h" & grid$kappa & grid$scale_id == "det")

t0 <- proc.time()[3]
res <- list()
for (dims in list(c(3, 5), c(5, 5))) {
  p <- dims[1]; q <- dims[2]; S <- true_sigma(p, q, "full")
  X <- gen_data(n, S$Sr, S$Sc, "col", 0.1)$X
  st <- make_init(X, "ff")                        # same start for every setting
  for (rl in rules) {
    fits <- lapply(seq_len(nrow(grid)), function(g)
      tff(X, rule = rl, beta = 0.75, denom = grid$denom[g], kappa = grid$kappa[g],
          scale_id = grid$scale_id[g], init = st, tol = 1e-10, maxit = 300))
    ref <- fits[[ref_id]]
    key_tr <- function(f) sapply(seq_len(f$iter), function(k)
      paste(as.integer(c(f$a_trace[[k]], f$b_trace[[k]])), collapse = ""))
    kref <- key_tr(ref)
    for (g in seq_len(nrow(grid))) {
      f <- fits[[g]]; kf <- key_tr(f); m <- min(length(kf), length(kref))
      res[[length(res) + 1]] <- data.frame(
        p = p, q = q, rule = rl, grid[g, ], status = f$status, iter = f$iter,
        same_final_a = identical(f$a, ref$a), same_final_b = identical(f$b, ref$b),
        same_path = all(kf[seq_len(m)] == kref[seq_len(m)]),
        shape_prod = shape_err(f$Sc %x% f$Sr, ref$Sc %x% ref$Sr),
        shape_r = shape_err(f$Sr, ref$Sr), shape_c = shape_err(f$Sc, ref$Sc),
        scale_vs_ref = tr_ratio(f$Sr, f$Sc, ref$Sr, ref$Sc),
        logdet_Sc = logdet(f$Sc))
    }
  }
}
res <- do.call(rbind, res); rownames(res) <- NULL
elapsed <- proc.time()[3] - t0
saveRDS(list(res = res, elapsed = elapsed), "output/step2.rds")

summ <- aggregate(cbind(all_same_a = same_final_a & same_final_b, same_path,
                        max_shape_prod = shape_prod, max_shape_r = shape_r, max_shape_c = shape_c) ~ p + q + rule,
                  res, function(x) if (is.logical(x)) mean(x) else max(x))
cat(sprintf("Step 2 done in %.1f s\n", elapsed))
print(summ, digits = 3)
cat("\nScale of product relative to reference (denom h, kappa on, det):\n")
print(xtabs(scale_vs_ref ~ rule + interaction(denom, kappa), res[res$scale_id == "det" & res$p == 3, ]), digits = 4)
cat("\nScale drift of Sigma^c with scale_id = 'none' (log det Sc at the end, p = 3, q = 5):\n")
print(res[res$scale_id == "none" & res$p == 3, c("rule", "denom", "kappa", "iter", "logdet_Sc")], digits = 3)

png("output/step2_shape.png", 1000, 450, res = 110)
par(mfrow = c(1, 2), mar = c(4, 7, 2.5, 1))
for (pp in c(3, 5)) {
  r <- res[res$p == pp, ]
  stripchart(log10(pmax(shape_prod, 1e-16)) ~ rule, r, method = "jitter", pch = 19, las = 1,
             col = ifelse(r$scale_id == "none", "tomato", "grey30"),
             xlab = "log10 shape error of product vs reference", main = sprintf("p = %d, q = 5", pp))
}
dev.off()
