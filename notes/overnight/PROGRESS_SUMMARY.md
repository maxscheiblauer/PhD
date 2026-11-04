# Overnight progress summary (owner-facing)

Each block appends one short dated entry: results with numbers, decisions, dead ends. Newest at the bottom.

---

**Block 0: planning (2026-11-04).** Infrastructure set up: setup script, cost meter, stop guard, protocol, and brief. The literature dossier is `notes/00_setting_matrix_cellwise.md`. No research results yet.

**Block 1 (2026-11-04, 12:54–13:16 UTC, ≈ $9).** Method defined and implemented: **cellMMCD**, the cellMCD objective (trimmed observed-cell likelihood, no imputation) with the Kronecker constraint. The Schur identity det Σ_OO = det Σ · det Ω_MM makes it O(p²q + pq² + m³) per observation instead of O((pq)³). Algorithm: exact cell step + score fixed-point parameter step; monotone and convergent (proven, tested; boundedness under the quota still conjectured). Pairwise Kronecker start `kpair`. 8 new tests pass, core 31/31. Benchmark (n = 100, 6-sd cells): recall 1.00, false flags 1–2 %; KL at 5×5, 10 % bad cells: cellMMCD 0.32 vs cellMCD-on-vec 0.59, RJ 89, MMCD 119, MLE 144; at 20 %: 0.44 vs 1.14 and > 280. Clean-data cost: KL 0.37 vs MLE 0.16, mostly a 10 % scale shrink (tr ratio 0.895) from the 1 % tail trimming; the naive 1-D factor (0.925) does not explain it. Run time 0.1–0.7 s at n = 100, 10 s at 10×10 (pure R). Dead ends: product-set deletion and a trimmed pseudo-likelihood (rejected on paper, journal). Next: consistency factor inside the iteration, existence condition, casewise + cellwise mix.
