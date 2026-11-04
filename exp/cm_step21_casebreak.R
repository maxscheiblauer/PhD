# Block 7, step 21: casewise breakdown of the reweighted hard quota (STATE-6 item 4) against the raw hard quota and the
# soft default. p = q = 5, n = 100, 5 reps. Casewise fraction eps of the matrices, types: shift (+5 sd in every cell,
# as cm_step13) and struct (AR(-0.5) x AR(-0.7), as cm_step16), eps from 0.10 to 0.45. Seed 21.
# Output: output/cm_step21_casebreak.log / .rds
source("R/load.R"); library(parallel); set.seed(21)
p <- q <- 5; n <- 100; R <- 5; S <- true_sigma(p, q, "full"); sd_cell <- sqrt(outer(diag(S$Sr), diag(S$Sc)))
one <- function(seed, type, eps) { set.seed(seed); X <- rmatnorm(n, Sr = S$Sr, Sc = S$Sc); ib <- seq_len(round(eps * n))
  if (type == "shift") X[, , ib] <- X[, , ib] + 5 * sd_cell[rep(TRUE, p * q)] else X[, , ib] <- rmatnorm(length(ib), Sr = ar_cov(p, -0.5), Sc = ar_cov(q, -0.7))
  bad <- seq_len(n) %in% ib; out <- c()
  for (variant in c("soft", "hard_raw", "hard_rw")) {
    f <- tryCatch(switch(variant, soft = cellmmcd(X), hard_raw = cellmmcd(X, quant_obs = NULL, quant_rw_obs = FALSE), hard_rw = cellmmcd(X, quant_obs = NULL)), error = function(e) NULL)
    out <- c(out, if (is.null(f)) c(NA, NA, NA, 1) else c(kl_mn(f$Sr, f$Sc, S$Sr, S$Sc), mean(!f$v[bad]), mean(!f$v[!bad]), 0)) }
  names(out) <- paste(rep(c("soft", "hard_raw", "hard_rw"), each = 4), c("kl", "trim_bad", "trim_good", "err"), sep = "."); out }
grid <- expand.grid(eps = c(0.10, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45), type = c("shift", "struct"))
res <- list(); sink("output/cm_step21_casebreak.log")
cat("cellMMCD step 21 (block 7): casewise breakdown, p = q = 5, n = 100, 5 reps. Per variant: median KL [max] / trimmed bad / trimmed good matrices.\n")
for (g in seq_len(nrow(grid))) {
  r <- do.call(rbind, mclapply(1:R, function(s) one(21000 + 10 * g + s, as.character(grid$type[g]), grid$eps[g]), mc.cores = 4)); res[[g]] <- list(grid = grid[g, ], r = r)
  line <- sprintf("%-6s eps %.2f:", grid$type[g], grid$eps[g])
  for (v in c("soft", "hard_raw", "hard_rw")) line <- paste0(line, sprintf(" | %s KL %7.2f [%7.2f] tb %.2f tg %.2f", v, median(r[, paste0(v, ".kl")], na.rm = TRUE), max(r[, paste0(v, ".kl")], na.rm = TRUE), median(r[, paste0(v, ".trim_bad")], na.rm = TRUE), median(r[, paste0(v, ".trim_good")], na.rm = TRUE)))
  cat(line, "\n") }
sink(); saveRDS(res, "output/cm_step21_casebreak.rds"); cat("done\n")
