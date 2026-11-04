# STATE: baton between overnight blocks

_Last updated: 2026-11-04 by block 4 (session_01Fi3wVbjjQDXtd8HKUf1Bfi). Each block rewrites this file at its end (PROTOCOL §6)._

## Current method: cellMMCD, reweighted, with shadow step and count rule (status: implemented, tested 20/20 + core 31/31, benchmarked cellwise, casewise, mixed, timed to 20×20)
- **Definition** (`notes/cellmmcd/METHOD.md` §2, §4, §4b–4e): trimmed observed-cell likelihood with Kronecker Σ, penalty λ_kl per flagged cell (calibrated), quota h = ⌈αn_v⌉ per position, observation weights v_i with the kept-cell threshold μ_i, the majority rule (≤ 50 % flagged) **and (block 4) the adaptive count rule** (≤ qbinom(.999, pq, ε̄) flagged cells, ε̄ = current overall flag rate), **shadow cell step** (block 4: during the first 10 iterations the cells of trimmed matrices are re-decided too, so a matrix trimmed only because its cellwise outliers were unflagged is re-admitted), then the reweighting step (block 3: re-flag at .999 on the raw fit, refit with flags fixed, scale factors t_rw, t_obs).
- **Why the block-4 changes:** the observation step judged matrices at W = Wfix and locked them out (25 % of matrices trimmed under 10 % cells at 3–4 sd, 4 % of clean ones); the shadow step fixes that (trimmed 0 %, clean KL 0.198 → 0.182) but lets scaled/structural casewise outliers back in through selection (their worst cells flagged, the rest passes); the count rule restores their trimming (scale 100 %, struct 75 / 60 % at 10 / 20 %). All in JOURNAL (block 4) and METHOD §4e.
- **Results with the default (median KL, n = 100):** 3×4 / 5×5 clean 0.087 / 0.198 (MMCD 0.087 / 0.183, MLE 0.070 / 0.155); 5 / 10 / 20 % cells at 6 sd: 0.083 / 0.094 / 0.139 (3×4), 0.184 / 0.205 / 0.320 (5×5) vs cellMCD-on-vec 0.26 / 0.26 / 0.32 and 0.63 / 0.59 / 1.14, MMCD and MLE ≥ 14. Casewise 5×5: shift 10 / 20 % 0.17 / 0.22, scale 0.19 / 0.22, struct 0.18 / 0.32 (MMCD 0.18–0.27, 0.23 for struct 20 %), mixed 0.22 (MMCD 51). 10 % cells at 3 / 4 sd in every matrix: 0.33 / 0.26, trimmed 0, recall 0.92 / 1.00. Clean matrices trimmed 0 % everywhere.
- **Theory added (block 4):** influence functions of the raw and reweighted estimator in the 1-D proxy (METHOD §7a; proven in the proxy, checked to 1e-4 by exact contaminated functionals, test 19): bounded, gross-error sensitivity 7.4 (raw) / 10.0 (reweighted), ≈ 0 beyond the cutoff; MLE unbounded.
- **Speed (block 4):** cell step and reweighting pass in C++ (`src/cm_wstep.cpp`, RcppArmadillo, compiled by `R/load.R` into `output/cpp_cache/`, falls back to R; identical results, test 20). n = 100: 0.3 / 0.4 / 2.2 / 8.6 s at p = q = 5 / 10 / 15 / 20 (R loops: 0.8 / 2.7 / 12.4 / 38.8; MMCD 0.2–4.7 s). METHOD §5 has the table.
- **Operating range:** n ≥ 3·max(p, q) (degenerate descents in 20 % / 5 % of clean 4×4 runs at n = 8 / 10, none at n ≥ 12; `cm_step11_smalln.log`).

## Code and tests (all passing at commit of this file)
- `R/cellmmcd.R`: 1 kpair, 2 objective, 3 cell step `cm_wstep(..., shadow, engine)`, 3a `cm_fullcol` (off), 3b `cm_vstep(..., max_flag, count_max)`, 4 `cm_pstep(tau)`, 5 driver `cellmmcd(..., shadow_it = 10, flag_quant = 0.999, reweight = TRUE, ...)`, 5b `cm_r2`, `cm_reweight` (C++ path), 6 `cm_proxy_fixpt`, 7 `cm_proxy_T`, `cm_proxy_if`.
- `src/cm_wstep.cpp`: `cm_wstep_cpp`, `cm_reweight_cpp`. `R/load.R` compiles it (`has_cpp`, `options(cellmmcd.engine)`).
- `tests/test_cellmmcd.R`: 20 tests (block 3's 17 + shadow step, proxy IF, C++ = R). `tests/test_core.R`: 31/31.
- Experiments (logs in `output/`): `cm_step9_interplay.R` (lock-in; `_nocount.log` = before the count rule), `cm_step10_timing.R` (`_Rengine.log` = R loops), `cm_step11_smalln.R`; reruns of `cm_step1.R` (block-3 log kept as `cm_step1_block3.log`; rerun predates the count rule, which does not act there) and `cm_step4_mixed.R` (`_block3.log`, `_shadow_only.log` kept).
- `tools/overnight/hb.sh`: block number and session id hard-coded (set to 4: edit for block 5).

## Prioritized next steps (concrete first actions)
1. **Structural casewise outliers at 20 %** (KL 0.32 vs MMCD 0.23; only 60 % trimmed): first action: in `cm_step4_mixed.R` print, for struct matrices, the distribution of flagged counts and of d_i − μ_i at the final fit; then try the kept-cell test at a stricter quantile for matrices whose count exceeds qbinom(.9, pq, ε̄) (two-level rule), or charge flagged cells at δ = λ − E[Δ | clean] in μ_i. Keep the cellwise settings of `cm_step9` as the guard (trimmed must stay 0).
2. **Monotone formulation of the shadow phase** (METHOD §8 2b): F is monotone only after iteration 10 by construction; a trimmed matrix's charge μ_i(W_i) is the non-monotone part. Idea: charge trimmed matrices μ_i evaluated at W = Wfix (constant), and let the observation step compare d_i(W_i) − μ_i(W_i) as now (then the v-step is no longer the exact minimiser of F, but F has a fixed trimmed charge). Write it down in JOURNAL first.
3. **Monotone consistent formulation** (censored-cell objective, METHOD §8.2): only if the reweighting device is judged insufficient for the paper.
4. **C++ for `cm_pstep` and `cm_obj_i`** if pq > 20 matters (8.6 s at 20×20 now, mostly R loops over observations).
5. **Location M** jointly from the kept cells; the median plug-in is used (low priority: shifts are handled by trimming).
6. **Hard-quota variant:** apply κ(α, pq) (`kappa_mcd`) to its scale.
7. **Influence function in the Kronecker model** (beyond the proxy): the IF of the reweighted estimator given the flags is the IF of an observed-cell MLE; the proxy result (§7a) covers the cell-rule part. Breakdown under cellwise contamination still open.
8. METHOD.md to paper form; plots (`cm_step9` comparison, timing curve, IF curve: `cm_proxy_if` on a grid).

## Open questions and pitfalls
- The count rule's ε̄ is the *current* flag rate among kept matrices: under heavy cellwise contamination (≥ 20 %) it is large and the rule becomes inactive (by design); under heavy casewise contamination the casewise matrices inflate ε̄ until they are trimmed (iteration helps: trimmed matrices leave ε̄).
- The reweighted estimator keeps a κ(.999)²-small shape bias; t_obs is evaluated at d = pq for all matrices.
- Casewise shifts of 2 sd in every cell are missed (recall 0.25 in `cm_step8`); 5 sd are caught 100 %.
- A `pgrep -f <script>` inside a waiting loop matches the loop's own shell; wait on a "done" marker in the script's stdout instead. `pkill -f` kills the tool's own shell.
- Each block's start-up reading costs ≈ $2–3 of the cap; keep reads to this file, the brief, METHOD §4e/§5/§7a, and the block-4 journal entries.

## Files to read first
`notes/overnight/RESEARCH_BRIEF.md`, this file, `notes/cellmmcd/METHOD.md` (§4c–4e, §5, §7–9), the block-4 entries of `notes/cellmmcd/JOURNAL.md`, `R/cellmmcd.R` (sections 3b and 5), `output/cm_step4_mixed.log`, `output/cm_step9_interplay.log`.
