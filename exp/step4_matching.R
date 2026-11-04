# Step 4 (Q4/Q5): compare matching rules R0-R3 and the joint rule RJ.
# Needs output/step1.rds and output/step3.rds. Run from project root: Rscript exp/step4_matching.R
source("R/load.R")
s1 <- readRDS("output/step1.rds"); s3 <- readRDS("output/step3.rds")
runs <- s1$runs; groups <- s1$groups
rules <- c("R0", "R1", "R2", "R3", "RJ")

conv <- aggregate(cbind(converged = status == "converged", cycle = status == "cycle",
                        J_nonmono = n_Jinc > 0, iter = iter) ~ rule + beta, runs, mean)
stab <- aggregate(cbind(n_distinct_a, sKL_med_starts = sKL_med) ~ rule + beta, groups, median)
acc <- xtabs(KL_shape_best ~ rule + cond,
             aggregate(KL_shape_best ~ rule + cond, groups[groups$beta == 0.75, ], median))
root_spread <- aggregate(cbind(sKL_max_roots = sKL_max) ~ rule, s3$res[s3$res$beta == 0.75, ], median)
tab <- merge(merge(conv, stab), root_spread, all.x = TRUE)

## RJu (unconstrained C-step) on the Step 1 data: how often does J* run off to -infinity?
set.seed(101)                                   # regenerates exactly the Step 1 datasets
p <- q <- 5; n <- 25; S <- true_sigma(p, q, "full")
conds <- c("clean", "col", "obs", "wcol"); nds <- 5
starts <- c("identity", "ff", rep("elemental", 8))
data <- lapply(conds, function(cd) lapply(seq_len(nds), function(ds)
  gen_data(n, S$Sr, S$Sc, cd, eps = if (cd == "clean") 0 else 0.1)))
t0 <- proc.time()[3]
rju <- do.call(rbind, lapply(seq_along(conds), function(ci) do.call(rbind, lapply(seq_len(nds), function(ds)
  do.call(rbind, lapply(c(0.75, 0.9), function(be) do.call(rbind, lapply(seq_along(starts), function(s) {
    set.seed(10000 * ds + s)
    f <- tff(data[[ci]][[ds]]$X, rule = "RJu", beta = be, init = starts[s], trace = FALSE,
             kappa = FALSE, scale_id = "none")
    data.frame(cond = conds[ci], beta = be, status = f$status, iter = f$iter)
  }))))))))
elapsed <- proc.time()[3] - t0
rju_tab <- with(rju, table(paste(cond, beta), status))

saveRDS(list(tab = tab, acc = acc, rju = rju, rju_tab = rju_tab, elapsed = elapsed), "output/step4.rds")

png("output/step4_accuracy.png", 1300, 420, res = 110)
par(mfrow = c(1, 4), mar = c(4, 4.5, 2.5, 0.5))
for (cd in conds) {
  g <- groups[groups$cond == cd & groups$beta == 0.75, ]; g$rule <- factor(g$rule, rules)
  boxplot(KL_shape_best ~ rule, g, log = "y", col = "grey85", ylab = "shape KL (best of 10 starts)",
          main = sprintf("%s, beta = .75", cd))
}
dev.off()

cat(sprintf("Step 4 (RJu part) done in %.1f s\n", elapsed))
print(tab, digits = 3, row.names = FALSE)
cat("\nMedian shape-KL of best-J start, beta = 0.75 (Step 1 data):\n"); print(acc, digits = 3)
cat("\nRJu status (10 starts x 5 datasets per row):\n"); print(rju_tab)
