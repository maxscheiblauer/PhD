# STATE: baton between overnight blocks

_Last updated: 2026-11-04 by the planning session (block 0). Each block rewrites this file at its end (PROTOCOL §6)._

## Status
- No cellwise Kronecker method exists yet. The starting points are:
  - the owner's slice-trimmed likelihood **RJ** (`R/tff.R`, `exp/RJ_simple.R`, `output/FINDINGS.md`): robust to column-type contamination, weak against row-wise contamination, many local optima;
  - `robustmatrix::mmcd` (casewise only).
- Environment: `bash tools/overnight/setup.sh` works from a clean container (R 4.3.3, robustmatrix 0.1.5, cellWise 2.5.7), and `Rscript tests/test_core.R` passes 31/31.

## Next steps for block 1 (orientation, then decide)
1. Read `RESEARCH_BRIEF.md`, then `output/FINDINGS.md` and `notes/rj_objective.md`. Read the MMCD C-step and flip-flop in the robustmatrix source (`src/functions.cpp`: `mmleCpp`, `cstepCpp`, `mmcd`).
2. Skim `notes/00_setting_matrix_cellwise.md` §3–§5 for the literature map. Keep the owner's guidelines in mind: no imputation, EM as a last resort.
3. Write `notes/cellmmcd/JOURNAL.md`: candidate formulations, their expected complexity, and a reasoned choice of what to try first.
4. Start implementing and testing the first candidate on small problems (p, q ≤ 5).

## Files to read first
`notes/overnight/RESEARCH_BRIEF.md`, `output/FINDINGS.md`, `notes/rj_objective.md`, `R/tff.R`
