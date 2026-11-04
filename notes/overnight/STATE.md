# STATE: baton between overnight blocks

_Last updated: 2026-11-04 by block 6 (session_01Fi3wVfbEYSFfctE82FNt4U). Each block rewrites this file at its end (PROTOCOL §6)._

## Current method: cellMMCD, reweighted, shadow step, count rule (robust rate), joint location; fallback: reweighted hard quota (status: implemented, tested 25/25 + core 31/31, benchmarks rerun, paper-form summary in METHOD §0)
- **Definition** (`notes/cellmmcd/METHOD.md` §0 summary, §2, §4, §4b–4e): trimmed observed-cell likelihood with Kronecker Σ, calibrated penalty λ_kl, per-position quota h = ⌈αn⌉, forced marginal flags (3 MAD), observation weights v (kept-cell threshold μ_i, majority rule, count rule with the robust rate), shadow cell step, EM location step, then the reweighting (re-flag at .999, refit with flags fixed, scale factors t_rw·t_obs). **Unchanged default in block 6.**
- **Block-6 additions:** (1) `quant_rw_obs` (§4c): the hard quota (`quant_obs = NULL`) now re-decides the observation weights in the reweighting step by the kept-cell distance d_i² ≤ χ²_{|O_i|,.99} plus the raw count rule, no quota (`cm_d2_obs`, Schur form; test 24). Removes the raw hard quota's 12 % scale over-correction and its clean cost (0.275 → 0.200; soft 0.190); within 0.03 KL of the soft default on every casewise and cellwise case; **struct 20 % + 5 % cells 0.35 vs 0.79 (soft)**; recommended fallback. (2) C++ cell step with cached inverses (bordering / downdate of Ω_MM⁻¹, O(m²) per candidate; old routine kept as `cm_wstep_cpp0`, test 25): identical flags, 20×20 fit 8.6 → 2.4 s; faster than MMCD from 10×10.
- **Theory added (block 6, JOURNAL results 4 and 6; METHOD §7, §7b):** chain-rule IF of the reweighted estimator in the Kronecker model, IF = J⁻¹[v s(z_O) + K·IF_raw(z)], bounded (the marginal cut bounds every kept value; the start has bounded IF), finite limit for gross cells and gross matrices (proven under regularity; K not in closed form). Breakdown: explosion by gross cells impossible below 1/2 contaminated cells per position (proven via the marginal cut; checked: no breakdown at 45 % gross cells per position, `cm_step19_gross.log`); the block-5 "35–40 % for gross cells" was an artefact of random placement; the quota 1 − α governs a *bounded bias* blow-up for moderate cells, not explosion.
- **Checked numerically (block 6):** the full default is consistent at the Gaussian: n·KL 8.3 / 8.0 / 7.4 at n = 500 / 2000 / 8000 (MLE 8.9 / 7.4 / 5.9; `cm_step15_consistency.log`); raw hard quota unbiased but 2.1–2.8× the MLE, reweighted 0.98–1.41. Sensitivity curves (`cm_step18_sc.log`, `cm_fig_sc.png`, 3×3, n = 500): bounded, flat beyond the cutoff (one cell: 4.95 for x ≥ 3 vs MLE 824 at x = 30; shifted matrix: 5.5 for x ≥ 4 vs 3975), maximum just inside the cutoff.
- **Closed:** C++ parameter step (5 % of the run time; profile in JOURNAL).
- **Results with the default (median KL, n = 100), unchanged:** cellwise 3×4 / 5×5 clean 0.084 / 0.195, 5 / 10 / 20 % 6-sd cells 0.089 / 0.195, 0.092 / 0.183, 0.130 / 0.271; casewise 5×5 shift / scale / struct 20 % 0.17 / 0.22 / 0.23 (MMCD 0.21 / 0.27 / 0.23), shift + cells 0.22 (MMCD 51), struct + cells 0.79 (hard rw 0.35, MMCD 48).
- **Figures:** block 5 `cm_fig_breakdown.png`, `cm_fig_if.png`, `cm_fig_mixed.png`; block 6 `cm_fig_timing.png` (regenerated), `cm_fig_sc.png` (`exp/cm_plots_block6.R`).

## Code and tests (all passing at commit of this file)
- `R/cellmmcd.R`: driver `cellmmcd(..., quant_rw_obs = NULL)` (NULL = automatic: .99 for the hard quota, off otherwise); `cm_reweight(fit, X, quant_rw, maxit, tol, quant_rw_obs)` with the observation reweighting; `cm_d2_obs`. Fit carries `rw_obs`, `quant_rw_obs`, `t_obs`; `fit$raw$v` / `fit$raw$trimmed` keep the quota, `fit$v` is the reweighted set.
- `src/cm_wstep.cpp`: `cm_wstep_cpp` (cached inverses, block 6), `cm_wstep_cpp0` (block 4, test only), `cm_reweight_cpp`.
- `tests/test_cellmmcd.R`: 25 tests (tests 8 and 22 now check the quota on `fit$raw`; 24 reweighted hard quota; 25 cached cell step). `tests/test_core.R`: 31/31.
- Experiments added: `exp/cm_step15_consistency.R`, `cm_step16_hardrw.R`, `cm_step17_hardrw_cells.R`, `cm_step18_sc.R`, `cm_step19_gross.R`, `cm_plots_block6.R`; `cm_step10_timing.R` rerun (adds the hard rw column; earlier logs `_block4.log`, `_block6_before_cache.log`).
- `tools/overnight/hb.sh`: block number and session id hard-coded (set to 6: edit for block 7).

## Prioritized next steps (concrete first actions)
1. **Implosion / boundedness of F below the n/max(p,q) line (§4a converse)** is the main theoretical gap: the breakdown statement for moderate cells and the existence of the estimator both rest on it. First action: formulate the claim for the *local* minimiser from `kpair` on the feasible set with values in the marginal-cut box B (JOURNAL block 6, result 6): show that along a monotone descent from a start with bounded condition number, det Σ_OO terms keep the smallest eigenvalue away from 0 unless a kept configuration is degenerate; look at `cm_step11_smalln.log` for what the degenerate configurations look like. If a proof is out of reach, write the precise conjecture with the numerical evidence (n ≥ 3 max(p, q)).
2. **Paper draft from METHOD §0**: `notes/cellmmcd/PAPER.md` with the sections definition → algorithm → theory → simulation → limits, pulling the tables from §9 and the five figures; reorder §4a–4e (currently 4a, 4c, 4d, 4b, 4e). Cheap, owner-facing.
3. **Choice of default**: soft threshold vs reweighted hard quota. Evidence so far: equal within 0.03 KL everywhere except struct + cells (hard rw far better) and clean 5×5 (soft better by 0.01–0.03). A larger grid (3 seeds × n = 100, 200; 5×5, 8×8) with both would settle it; the decision belongs to the owner, so present both in the paper.
4. **Casewise breakdown of the reweighted hard quota** at > 25 % casewise contamination (the reweighting could re-admit outliers as the reweighted MCD does): rerun `cm_step13_breakdown.R`'s casewise rows with `quant_obs = NULL`.
5. **Speed, if pq > 20 matters:** `kpair` (9 %) and the R-level parameter step / objective are what remains; the cell step is no longer the bottleneck.
6. Timing beyond 20×20 (30×30, 40×40) now that the cell step is O(m²): one background run of `cm_step10_timing.R` with the grid extended.

## Open questions and pitfalls
- The observation reweighting of the hard quota uses `fit$count_max` from the raw fit; when the raw count rule is off (`flag_quant = NULL`) only the majority rule applies.
- The sensitivity curves are finite-sample (n = 500); the transmitted term K·IF_raw is invisible at that resolution, consistent with the proxy's 0.006 but not a measurement of K.
- Breakdown by moderate cells (inside 3 MAD) at the quota is a *bias* blow-up (KL 5 / 71 at 25 / 30 %), bounded by construction; say "bias breakdown" in the paper, not "explosion".
- Block-5 breakdown log (`cm_step13_breakdown.log`): the gross-cell rows used random placement; with position-wise placement there is no breakdown up to 45 %. Both facts are true; cite both.
- Inserting a row in the middle of an experiment grid changes the seeds of later rows: append new rows at the end.
- Each block's start-up reading costs ≈ $2; keep reads to this file, the brief, METHOD §0 and the block-6 JOURNAL entries.
- `pgrep -f` inside a waiting loop matches the loop's own shell; wait on a done marker. Foreground `sleep` is blocked: wait with a background `until` loop.

## Files to read first
`notes/overnight/RESEARCH_BRIEF.md`, this file, `notes/cellmmcd/METHOD.md` §0 (summary) and §7–8, the block-6 entries of `notes/cellmmcd/JOURNAL.md`, `R/cellmmcd.R` (sections 5–5b), `output/cm_step16_hardrw.log`, `output/cm_step15_consistency.log`.
