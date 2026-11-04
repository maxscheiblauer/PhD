# STATE: baton between overnight blocks

_Last updated: 2026-11-04 13:45 UTC by block 2 (session_01Fi3wVcgnq5CzuGXGF5MFJh). Each block rewrites this file at its end (PROTOCOL §6)._

## Current method: cellMMCD (status: implemented, tested 13/13 + core 31/31, benchmarked cellwise, casewise and mixed)
- **Definition** (`notes/cellmmcd/METHOD.md` §2): trimmed observed-cell likelihood with Kronecker Σ, penalty λ_kl per flagged cell, quota h = ⌈αn_v⌉ kept cells per position, **plus (block 2) observation weights v_i** with threshold μ_i on the kept cells of X_i (χ² quantile with |O_i| dof at the initial estimate), quota h_obs = ⌈α_obs n⌉, and the constraint that a kept matrix has at most half its cells flagged (`max_flag`). No imputation. Three exact/GEM blocks: v-step, cell step, parameter step. `quant_obs = NULL` gives the hard quota (MMCD-type), `alpha_obs = 1` switches observation trimming off.
- **Scale calibration (block 2, proven in the 1-D proxy, checked numerically):** λ built from c/t₁ + log t₁ with t₁ = F₃(c)/F₁(c) = 0.925, final Σ/(t₁ t_obs). Clean-data tr ratio 0.98–0.99 (was 0.895). `cm_proxy_fixpt()` gives the fixed points.
- **Existence (block 2, proven for data in general position):** F is unbounded below whenever n − h ≥ n/max(p,q) (so at α = .75 for max(p,q) ≥ 4). The estimator is the local minimiser from the `kpair` start; at n = 100 the descent never degenerated in any experiment, at n ≤ 12 it does (`output/cm_step5_exist.log`). Converse and a full-column remedy conjectured (JOURNAL, METHOD §4a).
- **Results (p = q = 5, n = 100, median KL; `output/cm_step4_mixed.log`, `cm_step1.log`):** clean 0.31 (MMCD 0.18, MLE 0.16); 5 % cells 0.29 (MMCD 45); 10 % / 20 % shifted matrices 0.24 / 0.37 (MMCD 0.18 / 0.21; cells-only 0.36 / 2.48); scale ×3 0.25 / 0.28; structural 0.20 / 0.30; 10 % shifted + 5 % cells 0.25 (MMCD 51). Trims 100 % of shifted/scaled, 85–95 % of structural matrices, 0–3 % of clean ones.
- **Known weakness:** efficiency at the Gaussian. KL ratio to the MLE 2.4 at n = 1000 (3×4) with quant = .99; 1.3 with quant = .999. Not caused by the start or by λ noise (oracle start gives the same KL). Unexplained (`output/cm_step3c_eff.log`).

## Code and tests (all passing at commit of this file)
- `R/cellmmcd.R`: `kpair`, `cm_obj_i`, `cm_obj`, `cm_mu_i`, `cm_delta1`, `cm_wstep`, `cm_vstep`, `cm_pstep`, `cellmmcd`, `cm_proxy_fixpt`.
- `tests/test_cellmmcd.R`: 13 tests (block 1's 8 + proxy fixed point, calibrated scale, three-block monotonicity, casewise detection, brute-force v-step). `tests/test_core.R`: 31/31.
- Experiments: `exp/cm_step1.R` (cellwise benchmark, block 1, uncalibrated), `cm_step2_factor.R`, `cm_step3_calib.R`, `cm_step3b_lambda.R`, `cm_step3c_eff.R`, `cm_step4_mixed.R`, `cm_step5_exist.R`; logs under `output/`. `cm_step1` has not been rerun with the calibrated default and the observation step.
- `tools/overnight/hb.sh <status>`: block number and session id are hard-coded (set to 2: edit for block 3).

## Prioritized next steps (concrete first actions)
1. **Rerun the cellwise benchmark with the current default** (`Rscript exp/cm_step1.R`, ~3 min) so that METHOD §9's first table reflects the calibrated estimator with the observation step; update the table.
2. **Efficiency at the Gaussian.** First action: in `exp/cm_step3c_eff.R` add a variant that runs the default and then one extra cell + parameter step at quant = .999 (a cellMCD-style final refinement) and a variant with `calib = FALSE` but final division by the empirical fixed point; compare the KL ratio to the MLE. Then look for the cause: compute the influence function of the proxy M-estimator, or check whether the loss is in the shape (it should not be) via `kl_shape`.
3. **Degenerate descents at small n.** Implement a safeguard: after each cell step, for every l check the number of fully kept columns (i,l); if it is below p, unflag the cells with the smallest Δ − λ until it is p (same for rows and q). Test on `exp/cm_step5_exist.R`; record whether monotonicity is lost. Alternatively a stopping rule on log det Ω_MM.
4. **Hard-quota variant:** apply κ(α, pq) (`kappa_mcd` in `R/matnorm.R`) to its scale and compare with the default on the step-4 grid; it is the natural bridge to MMCD and RJ.
5. **Location M** jointly from the kept cells (score fixed point); check the median plug-in effect under 20 % shifted matrices (the median moves).
6. **Speed**: `cm_wstep` loops in R over n·pq cells; vectorise the m = 0 case or Rcpp; then timings for p = q = 10, 20.
7. Plots for `output/` (none yet); then METHOD.md to paper form (the counting proof deserves a clean statement).

## Open questions and pitfalls
- Check the cellMCD paper (Raymaekers & Rousseeuw 2024, Section 2) for how they handle the unbounded objective; the counting argument applies to their objective too. arXiv/CRAN were not consulted in block 2.
- t_obs is evaluated at d = pq for all matrices (approximation for matrices with flagged cells).
- With `max_flag` the monotonicity of F holds only between transitions where a kept matrix crosses the bound (never observed to break in tests).
- λ and μ_i are fixed from the `kpair` start; the oracle-start experiment shows this costs nothing at n = 1000 on clean data.
- Don't compare raw factors; `scale_ident(..., "det")` is applied at the end of `cellmmcd`. `tff_multi(X, rule = "RJ", beta = .75)` already passes trace = FALSE.
- `pkill -f <script>` from a Bash tool call kills the tool's own shell (matched its command line); use `pkill -f "Rscript exp/..."` with care or TaskStop.

## Files to read first
`notes/overnight/RESEARCH_BRIEF.md`, this file, `notes/cellmmcd/METHOD.md`, the block-2 entries of `notes/cellmmcd/JOURNAL.md`, `R/cellmmcd.R` (sections 3b and 5), `output/cm_step4_mixed.log`, `output/cm_step3c_eff.log`.
