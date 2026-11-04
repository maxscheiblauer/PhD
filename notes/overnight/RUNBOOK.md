# Overnight run: runbook

## Architecture
```
planning session (this one)
   │  on "go": create_session(supervisor)
   ▼
SUPERVISOR  (Opus 5.5, slim context, no research)        ◄── hourly watchdog scheduled task
   │  launches one worker at a time; reads spend; interrupts runaways; ensures final report
   ▼
WORKER block 1 ──BLOCK_DONE──► WORKER block 2 ──► … ──► WORKER final (MORNING_REPORT.md)
   (smart model; fresh container; state handed over via git: STATE.md, heartbeat.json, LEDGER.csv)
```

**Why a relay of fresh sessions instead of one long session**
- *Bounded context.* Every tool call re-reads the whole context, so one session that runs all night gets more expensive per step as it goes. Fresh blocks start small.
- *Container loss.* A crash or timeout costs at most one block. Everything else is already pushed.
- *Real cost numbers.* `get_session` reports authoritative spend once a worker's turn has ended. The supervisor sums these numbers.
- *Forced documentation.* Every block must leave a baton (`STATE.md`) and a summary entry, so there is always a current summary, even if the run ends abruptly.

**Safety nets**
| Risk | Mechanism |
|---|---|
| Worker stops early or idles waiting for the user | Stop guard hook (`tools/overnight/stop_guard.sh`, armed only in worker containers) |
| Worker overspends inside its turn | Local meter (`tools/overnight/meter.py`, calibrated to slightly overestimate) → wrap-up at 85% of the cap; supervisor sends a wrap-up message at 1.3× the cap and interrupts at 1.6× |
| Worker crashes, or its container dies | Supervisor gets a child-session event or a watchdog wake → relaunches from the last pushed state |
| Credits run out before the summary | Supervisor reserves `SUMMARY_RESERVE` for a final block; `PROGRESS_SUMMARY.md` is updated after every block anyway |
| Fresh container has no R | `tools/overnight/setup.sh` (≈3.5 min, idempotent; tested 2026-11-04) |

## Budget (proposal; adjust before launch)
| Item | USD |
|---|---|
| Total credit | 100 |
| Planning session (this one, measured at launch) | ≈ 8 |
| Supervisor reserve | 5 |
| Final summary block | 6 |
| Safety margin | 5 |
| **Work blocks** | ≈ 76 → **5 blocks × $15** |

Prices (from the claude-api skill, cached 2026-10-06), USD per million tokens:

| Model | Input | Output | Cache read |
|---|---|---|---|
| Fable 5.1 | 10 | 50 | 0.25 |
| Opus 5.5 | 4 | 20 | 0.20 |

Cache writes are assumed to cost 2× input for the 1 h TTL. A $15 Fable block is very roughly 45–90 minutes of work, depending on how much it thinks versus runs code. The meter estimate was checked against `get_session` on 2026-11-04: $2.90 estimated vs $2.80 actual.

## Pre-flight checklist (owner)
1. **GitHub push access.** Done; pushes have worked since 2026-11-04.
2. **Network allowlist** (environment settings → Network access): `arxiv.org`, `export.arxiv.org`, `repositum.tuwien.at`, `wis.kuleuven.be`, `lirias.kuleuven.be`, `cran.r-project.org`, `cloud.r-project.org`.
   Optional, but literature work is much weaker without them.
3. Choose the **worker model** (the default proposal is `claude-fable-5-1`) and confirm the budget numbers.
4. Say "go" in the planning session. It measures its own cost, fills in `PLANNING_USD`, and creates the supervisor.

## Monitoring and abort
- **Progress:** `notes/overnight/PROGRESS_SUMMARY.md` and `heartbeat.json` on the branch, or the session list in the Claude app.
- **Stop everything:**
  - archive the supervisor session, or tell it "shut down";
  - delete its hourly scheduled task;
  - archive or interrupt the running worker.
- **Morning:** read `notes/overnight/MORNING_REPORT.md` first, then `PROGRESS_SUMMARY.md` and `notes/cellmmcd/`.
