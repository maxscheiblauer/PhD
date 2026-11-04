# Step 1 (Q1): convergence, cycling, monotonicity of J, variation across starts.
# Run from project root: Rscript exp/step1_convergence.R
source("R/load.R")
set.seed(101)
p <- q <- 5; n <- 25; S <- true_sigma(p, q, "full")
conds <- c("clean", "col", "obs", "wcol")     # contaminated: eps = 0.1, delta = 5
rules <- c("R0", "R1", "R3", "R2", "RJ")       # R2, RJ are used again in Step 4
betas <- c(0.75, 0.9); nds <- 5
starts <- c("identity", "ff", rep("elemental", 8))

data <- lapply(conds, function(cd) lapply(seq_len(nds), function(ds)
  gen_data(n, S$Sr, S$Sc, cd, eps = if (cd == "clean") 0 else 0.1)))
names(data) <- conds

t0 <- proc.time()[3]
runs <- list(); groups <- list(); traces <- list()
for (cd in conds) for (ds in seq_len(nds)) for (rl in rules) for (be in betas) {
  X <- data[[cd]][[ds]]$X
  keep_trace <- cd == "col" && ds == 1 && be == 0.75
  fits <- lapply(seq_along(starts), function(s) {
    set.seed(10000 * ds + s)                    # same elemental subsets for every rule
    # kappa off, no rescaling: each half-step is then the exact minimizer of its half of J
    tff(X, rule = rl, beta = be, init = starts[s], trace = keep_trace,
        kappa = FALSE, scale_id = "none", denom = "h")
  })
  if (keep_trace) traces[[rl]] <- fits
  ok <- !sapply(fits, function(f) f$status %in% c("singular", "degenerate"))
  Jf <- sapply(fits, function(f) tail(f$J, 1))
  kl <- sapply(fits, function(f) kl_mn(f$Sr, f$Sc, S$Sr, S$Sc))
  kls <- sapply(fits, function(f) kl_shape(f$Sr, f$Sc, S$Sr, S$Sc))
  runs[[length(runs) + 1]] <- data.frame(
    cond = cd, ds = ds, rule = rl, beta = be, start = seq_along(starts), init = starts,
    status = sapply(fits, `[[`, "status"), iter = sapply(fits, `[[`, "iter"),
    revisit = sapply(fits, `[[`, "revisit"),
    n_Jinc = sapply(fits, function(f) sum(diff(f$J) > 1e-9 * abs(f$J[-1]), na.rm = TRUE)),
    J = Jf, KL = kl, KL_shape = kls)
  # spread across starts: symmetrized KL of every start to the best-J start
  best <- which.min(ifelse(ok, Jf, Inf))
  sk <- sapply(fits[ok], function(f) skl(f, fits[[best]]))
  groups[[length(groups) + 1]] <- data.frame(
    cond = cd, ds = ds, rule = rl, beta = be,
    n_distinct_a = length(unique(sapply(fits[ok], function(f) paste(as.integer(f$a), collapse = "")))),
    sKL_med = median(sk), sKL_max = max(sk), KL_best = kl[best], KL_shape_best = kls[best],
    KL_shape_ff_start = kls[2])
}
runs <- do.call(rbind, runs); groups <- do.call(rbind, groups)
elapsed <- proc.time()[3] - t0

# status table: fraction converged / cycle / maxit (and maxit with revisited weights)
runs$status2 <- ifelse(runs$status == "maxit" & runs$revisit, "maxit+revisit", runs$status)
tab_status <- with(runs, table(paste(rule, beta), status2))
tab_mono <- aggregate(cbind(nonmonotone = n_Jinc > 0) ~ rule + beta, runs, mean)
tab_spread <- aggregate(cbind(n_distinct_a, sKL_med, sKL_max, KL_best) ~ rule + beta + cond, groups, median)

saveRDS(list(runs = runs, groups = groups, tab_status = tab_status, tab_mono = tab_mono,
             tab_spread = tab_spread, elapsed = elapsed,
             traces = lapply(traces, function(fl) lapply(fl, `[`, c("J", "delta", "nchange", "status")))),
        "output/step1.rds")

## plots
png("output/step1_J_traces.png", 1500, 400, res = 110)
par(mfrow = c(1, length(rules)), mar = c(4, 4, 2.5, 0.5))
for (rl in rules) {
  Js <- lapply(traces[[rl]], `[[`, "J")
  Jmin <- min(unlist(Js), na.rm = TRUE)
  yl <- range(log10(unlist(Js) - Jmin + 1), na.rm = TRUE)
  plot(NA, xlim = c(1, max(lengths(Js))), ylim = yl, xlab = "iteration",
       ylab = "log10(J - min J + 1)", main = paste(rl, "(col, eps=.1, beta=.75)"))
  for (s in seq_along(Js)) lines(log10(Js[[s]] - Jmin + 1), col = s, lwd = 1.3)
}
dev.off()

png("output/step1_status.png", 1100, 500, res = 110)
par(mar = c(6, 4, 2, 9), xpd = TRUE)
tt <- prop.table(tab_status, 1)
barplot(t(tt), las = 2, col = c("grey80", "tomato", "orange", "gold")[seq_len(ncol(tt))],
        ylab = "fraction of runs", main = "Step 1: run status (all conditions, 10 starts)")
legend("topright", inset = c(-0.22, 0), legend = colnames(tt),
       fill = c("grey80", "tomato", "orange", "gold")[seq_len(ncol(tt))], bty = "n")
dev.off()

cat(sprintf("Step 1 done in %.1f s, %d runs\n", elapsed, nrow(runs)))
print(tab_status); print(tab_mono)
print(aggregate(cbind(n_distinct_a, sKL_med, sKL_max, KL_best) ~ rule + beta, groups, median))
print(tab_spread[tab_spread$beta == 0.75, ])
