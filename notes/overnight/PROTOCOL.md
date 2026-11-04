# Overnight worker protocol

You are one **worker block** in a relay of autonomous Claude Code cloud sessions that develop the method in `RESEARCH_BRIEF.md`. **No human is watching.** Do not ask questions and do not wait for answers: decide, document the decision, and move on.

A supervisor session started you. It launches the next block after you hand over.

Your prompt gives you:
- `BLOCK` (block number);
- `CAP_USD` (cost cap for this block);
- `WALL_MIN` (wall-clock limit);
- `SUPERVISOR` (session id);
- `MODE` (`work` or `final`).

## 1. Start of block (first 5 minutes)
1. Mark this container as a worker; this arms the Stop guard in `.claude/settings.json`:
   `touch /tmp/overnight_worker; echo $CAP_USD > /tmp/overnight_block_cap; date -u +%FT%TZ > /tmp/overnight_block_start`
2. Start the environment setup **in the background**: `bash tools/overnight/setup.sh`. It installs R 4.3 and builds robustmatrix and cellWise; it takes about 3.5 min the first time.
3. `git pull` on the working branch. Read, in this order:
   - `notes/overnight/STATE.md` (the baton from the previous block);
   - `notes/overnight/RESEARCH_BRIEF.md`;
   - `notes/overnight/PROGRESS_SUMMARY.md`;
   - then only what STATE.md points to.

   Do not re-read the whole repo. Earlier blocks already did that, and reading is your most expensive start-up cost.
4. Write `notes/overnight/heartbeat.json` with `{"block": BLOCK, "session": <your session id, from get_session>, "status": "working", "meter_usd": <meter>, "utc": <now>}`, then commit and push. This also claims the block.

## 2. Budget discipline (the owner has a hard credit limit)
- **Meter:** `python3 -I tools/overnight/meter.py` estimates this container's spend. It is calibrated, with a safety factor of 1.15.
  The `cost_usd` from `get_session` only refreshes at turn boundaries, so inside your turn the meter is the authority.
- Check the meter at least every **20 minutes** or every ~25 tool calls, and at every commit.
- When the meter reaches **85% of CAP_USD**, or the elapsed time reaches **WALL_MIN**, stop starting new work and run the end-of-block procedure (§6) immediately. Never exceed CAP_USD by choice.
- **Context is money.** Every tool call re-reads the whole context.
  - Keep tool output short: redirect long output to files under `output/` and print a few summary lines (`tail`, `head`, `grep`).
  - Don't `cat` large files. Read only the line ranges you need.
  - Don't paste long code into your messages.
- **Subagents** (Agent tool) are allowed only for self-contained tasks, and they count toward the meter. Prefer a few large Bash runs over many small calls.
- **Simulations:** design them to finish in minutes. Use all 4 cores (`parallel::mclapply`). Run anything longer than ~5 min in the background with a log file, and keep working meanwhile.

## 3. Research and coding rules
- Follow the owner's guidelines in `RESEARCH_BRIEF.md`: no imputation-based method, and EM only as a last resort.
- Follow the repo conventions:
  - functions in `R/`, tests in `tests/`, experiments in `exp/`, logs and plots in `output/` (`*.rds` are git-ignored), theory in `notes/`;
  - a seed in every experiment.
- `Rscript tests/test_core.R` (plus any new test files) must pass before each commit. Add tests for new code.
- **Be honest.**
  - Report failures and negative results plainly, in `output/FINDINGS.md` style.
  - Do not tune toward nice numbers.
  - Mark every theoretical claim as **proven**, **checked numerically**, or **conjectured**.
- Compare covariances via the Kronecker product, or after shape normalization, never via the raw factors (scale is not identified).
- Benchmark against `robustmatrix::mmcd`, the owner's RJ, the plain flip-flop MLE, and, where feasible, cellwise methods on vec(X) (`cellWise::cellMCD`; only possible for small pq).
- Before deep work in a direction, write down the idea and why it should work in `notes/cellmmcd/JOURNAL.md` (append-only, dated entries). Abandoned directions stay there, with the reason.

## 4. Literature
- arXiv and author pages may or may not be reachable, depending on the environment's network allowlist. Try `curl -sI https://arxiv.org` once. If it is blocked, use WebSearch (abstracts and snippets) and the package sources on github.com/cran.
- Put PDFs in `references/` (they are git-ignored). Put reading notes in `references/` as markdown.

## 5. Persistence (the container is ephemeral)
- Commit and push at least every **30 minutes** and after every milestone.
- Each commit also updates `heartbeat.json` with the current meter value.
- Push to the working branch only. Retry on network errors (2 s, 4 s, 8 s, 16 s).
- If the push is refused for permissions, keep working, commit locally, and say so in STATE.md and in your message to the supervisor.

## 6. End-of-block procedure (also when stopping early for any reason)
1. Bring the code to a passing state, or revert to the last passing commit, and note what was reverted.
2. Rewrite `notes/overnight/STATE.md` as the baton for the next block, at most about 150 lines:
   - current method candidate(s) and their status;
   - what was done in this block;
   - the exact state of the code and tests;
   - a **prioritized next-steps list** with concrete first actions;
   - open questions and pitfalls;
   - the files to read first.
3. Append a short dated entry to `notes/overnight/PROGRESS_SUMMARY.md` (owner-facing, ≤ 15 lines): results with numbers, decisions, and dead ends.
4. Append a row to `notes/overnight/LEDGER.csv`: `block,session_id,start_utc,end_utc,meter_usd,summary`.
5. Set `heartbeat.json` to `"status": "closed"`. Commit and push.
6. `touch /tmp/overnight_block_closed` (this disarms the Stop guard).
7. `send_message` to SUPERVISOR: `BLOCK_DONE block=<k> meter_usd=<x> pushed=<yes|no> next=<one line>`. Then end your turn.

## 7. MODE = final
Do not start new research. Read STATE.md, PROGRESS_SUMMARY.md, JOURNAL.md and the latest FINDINGS. Then write `notes/overnight/MORNING_REPORT.md` for the owner (≤ 2 pages):
- what exists now: method, code, theory, results;
- what is solid and what is shaky;
- what failed;
- the recommended next steps;
- the total spend from LEDGER.csv.

Then run §6 steps 4–7.
