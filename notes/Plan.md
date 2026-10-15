# Plan: exploratory R study of column-trimmed flip-flop estimation

> **Historical (5 Oct 2026).** This is the original plan for steps 0–6, and it has been carried out. Current status: `README.md` and `output/FINDINGS.md`. Paths have changed since: results go to `output/` (not `out/`).

## Context

Read `trimmed_flipflop_spec.md` first. It defines the model, the whitened vectors $z_{ij}$ and $u_{ik}$, the weighted updates, and the selection/matching rules R0–R3.

This is exploratory groundwork. The goal is to find out **whether and where the idea breaks**, not to produce a polished method.

Keep everything small:
- $p,q \in \{3,5\}$
- $n \in \{10, 25, 50\}$
- 50 replications at most
- No large grids, no parallelization

Each experiment should run in well under a minute.

## Ground rules

- Base R, plus `MASS` if needed. Check whether the CRAN package `robustmatrix` is installable. It should provide `mmcd()` and a matrix-normal MLE and sampler; use it as the MMCD reference if available. Otherwise implement MMCD yourself as rule R2 with C-steps and multiple starts, and say so in the report.
- Set a seed in every script. Keep functions in `R/` and experiments in `exp/`. Each experiment writes a short `.rds` result plus one or two plots to `out/`.
- Do not tune things to make results look good. Report failures plainly.
- Never compare $\hat\Sigma^r$ or $\hat\Sigma^c$ individually without first fixing scale. Compare the product $\hat\Sigma = \hat\Sigma^c\otimes\hat\Sigma^r$, or shape-normalized factors $\Sigma/|\Sigma|^{1/\dim}$.

## Step 0: Core functions (`R/`)

1. `rmatnorm(n, M, Sr, Sc)`: sampler. Validate it by checking that the empirical covariance of $\operatorname{vec}(X)$ approaches $\Sigma^c\otimes\Sigma^r$.
2. `flipflop(X, ...)`: plain MLE. It must match `robustmatrix` if available.
3. `root(S, type = c("sym", "chol", "randrot"))` returns $A$ with $AA' = S^{-1}$.
4. `tff(X, rule, hr, hc, denom, kappa, scale_id, root_type, init, maxit, tol, trace = TRUE)`: the trimmed flip-flop.
   - It must return the estimates, the iteration count, a convergence flag, and a per-iteration trace with:
     - the objective $J$ (spec §6),
     - $\|\hat\Sigma^{(k)}-\hat\Sigma^{(k-1)}\|_F$ for the Kronecker product,
     - the number of weights that changed,
     - and the full weight matrices $a$ and $b$, so cycles can be detected.
   - Inits:
     - (a) identity,
     - (b) plain flip-flop,
     - (c) several random elemental subsets.
5. `kappa(beta, m)` as in spec §4.
6. Metrics:
   - KL divergence between $\mathcal N(0,\hat\Sigma)$ and $\mathcal N(0,\Sigma)$,
   - shape errors of each factor,
   - and, for contaminated data, the recall and precision of flagged $z_{ij}$ versus truly contaminated columns. Precision/recall is only meaningful when $\Sigma^c$ is diagonal; note this wherever it is reported.

7. Unit tests:
   - With no trimming, `tff` equals `flipflop`.
   - The identity $\sum_j d_{ij} = \sum_k e_{ik}$ from spec §3 holds numerically.

## Step 1: Convergence (Q1)

Use clean data and contaminated data (see the contamination list below), rules R0, R1 and R3, and $\beta \in \{0.75, 0.9\}$.

- Record how often each run converges, cycles (a repeated weight configuration), or hits `maxit`.
- Check whether $J$ is monotone. Plot $J$ traces for a few runs.
- Use 10 random starts. Report how much the final estimates vary across starts.

## Step 2: Normalization (Q2)

Fix a dataset and seed. Vary the denominator ($h$ vs. $nq$), $\kappa$ (on/off), and the scale identification (all options in spec §4).

**Test the conjecture in spec §7:**
- Are the selected sets identical across settings?
- Are the shape-normalized factors identical up to tolerance?
- If not, find out why. For example, κ might enter before the ranking.

## Step 3: Root dependence (spec Q3)

Run each rule with $A^c$ symmetric, Cholesky, and 5 random rotations.

- Quantify the spread of $\hat\Sigma$ (KL between the variants) under no trimming, which should be zero, and under trimming.
- This is likely the most fundamental flaw, so give it a clear verdict.

## Step 4: Matching rule (Q4/Q5)

Compare R0, R1, R2 and R3 on convergence (Step 1 metrics), accuracy (KL), and stability across starts.

- If time allows, try one principled alternative. Choose $a$ to minimize the full trimmed criterion jointly, then derive the $\Sigma^c$ step from the *same* $a$ by re-expressing the retained information in the $u$-coordinates.
- If such a rule cannot be written down, explain why in one paragraph.

## Step 5: Consistency (Q6)

Use clean data with $n \in \{10, 25, 50\}$ and $\beta \in \{0.5, 0.75, 0.9\}$.

- Estimate $E[\hat\Sigma]$ via the ratio $\operatorname{tr}(\Sigma^{-1}\hat\Sigma)/(pq)$, with and without $\kappa$.
- Compare the oracle version (true $\Sigma^c$ fixed, only $\Sigma^r$ estimated) to the full plug-in version to isolate the plug-in bias.
- Check whether the empirical $d_{ij}$ at convergence look $\chi^2_p$ (QQ plot for one setting).

## Step 6: Comparison with MMCD (Q7)

Use $p = q = 5$, $n = 25$, and contamination fractions $\varepsilon \in \{0, 0.1, 0.2\}$.

Contamination types:
- **observation-wise:** whole $X_i$ drawn from a shifted or inflated distribution
- **column-wise:** one original column of $X_i$ replaced by an outlier
- **row-wise:** the analogue of column-wise for rows
- **cellwise:** random single cells
- **whitened-column:** one $z_{ij}$ replaced, generated in whitened coordinates. This is the best case for the idea.

Use each of the R rules (R2 = MMCD) at matched breakdown, i.e. comparable $\beta$. Report the median KL and IQR. Present this as a small table and one boxplot per contamination type.

## Deliverable

`out/FINDINGS.md`, at most about 2 pages:
- one short section per question Q1–Q7 in the spec, each with a verdict (works / fails / unclear), the key number or plot, and the suspected reason;
- a final list of "flaws confirmed", "flaws refuted" and "new questions".

Include exact commands to reproduce each experiment. Keep the prose short.