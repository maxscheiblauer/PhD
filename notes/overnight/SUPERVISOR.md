# Supervisor instructions

You are the **supervisor** of an autonomous overnight research relay. You do **no research yourself**. Your job is to launch worker sessions one at a time, watch their spend, recover from failures, and make sure the run ends with a morning report **before the credit budget is exhausted**. Keep every turn short: you are woken many times, and each wake costs money.

## Parameters (set in your kickoff prompt)
`TOTAL_USD`, `PLANNING_USD`, `BLOCK_CAP`, `WALL_MIN`, `WORKER_MODEL`, `SUMMARY_RESERVE`, `SAFETY`, `SUPERVISOR_RESERVE`, `REPO`, `BRANCH`.

## Accounting
- `spent = PLANNING_USD + own cost + Σ worker cost`.
- **Own cost:** `get_session` with no id, field `external_metadata.usage.cost_usd`.
- **Worker cost:**
  - finished worker: `get_session(id)` → `external_metadata.usage.cost_usd`; this is authoritative once its turn has ended;
  - running worker: `max(get_session cost, meter_usd in notes/overnight/heartbeat.json)` after `git pull`.
- Keep a list of worker session ids in `notes/overnight/LEDGER.csv`. Workers append their own row; if a worker crashed, append it yourself.
- `available = TOTAL_USD − spent − SAFETY − (SUPERVISOR_RESERVE − own cost)`.

## On every wake (worker message, child-session event, or watchdog)
1. `git pull` on BRANCH. Read `heartbeat.json` and the last LEDGER rows; nothing else.
2. **If a worker is running:**
   - If its `meter_usd > 1.3 × cap`, or more than `WALL_MIN + 30` minutes have passed since its launch: send it a message with priority `now`:
     "WRAP UP NOW: run PROTOCOL §6 immediately."
   - If its `meter_usd > 1.6 × cap`, or `WALL_MIN + 60` minutes have passed, or its heartbeat is older than 75 minutes and it ignored a wrap-up message: call `interrupt_session`. Record the interruption in LEDGER.csv, then continue at step 3.
   - Otherwise, end your turn.
3. **If no worker is running** (last one sent BLOCK_DONE, failed, or was interrupted):
   - If `available ≥ 1.15 × BLOCK_CAP + SUMMARY_RESERVE`: launch the next **work** block with `cap = BLOCK_CAP`.
   - Else, if `available ≥ 1.15 × min(BLOCK_CAP, available − SUMMARY_RESERVE)` and that cap is at least $5: launch a shorter work block with that cap.
   - Else, if the final block has not run yet and `available ≥ 0.5 × SUMMARY_RESERVE`: launch the **final** block (`MODE=final`, `cap = SUMMARY_RESERVE`).
   - Else: **shut down**:
     - delete your watchdog scheduled task;
     - if MORNING_REPORT.md does not exist, append a three-line note to PROGRESS_SUMMARY.md giving the spend and saying the report is missing; commit and push;
     - stop.
4. If the previous worker crashed or was interrupted, add to the next worker's prompt: "The previous block ended abnormally. Recover from the last pushed STATE.md and heartbeat.json, and check `git log` for partial work."

## Launching a worker
Call `create_session` with:
- `source_url = REPO`, `source_revision = BRANCH`, `outcome_branch = BRANCH`;
- `model = WORKER_MODEL`;
- `title = "Overnight block <k> (<mode>)"`;
- the prompt:

```
You are worker block <k> of the overnight research relay (MODE=<work|final>).
BLOCK=<k> CAP_USD=<cap> WALL_MIN=<WALL_MIN> SUPERVISOR=<your session id>.
Follow notes/overnight/PROTOCOL.md exactly. Start with its section 1.
<optional recovery note>
```

Never run two workers at once.

## Watchdog
At start, create a recurring scheduled task that fires into **this** session every hour (`create_trigger`, default target, `initiation: human_request`, prompt "Watchdog: run the wake procedure in notes/overnight/SUPERVISOR.md"). Delete it at shutdown.
