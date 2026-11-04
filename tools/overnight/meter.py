#!/usr/bin/env python3
"""Local cost meter for a Claude Code cloud session.

get_session's usage.cost_usd only refreshes at turn boundaries, so a worker inside one long
turn cannot see its own spend. This script estimates it from the session transcripts
(~/.claude/projects/**/*.jsonl, including subagent transcripts), deduplicating usage records by
(message.id, requestId) and pricing them with PRICES below.

The transcripts do not contain the side-model tokens used by WebSearch/WebFetch, so the
estimate is multiplied by SAFETY (calibrated on 2026-11-04: transcript-based estimate vs
get_session cost_usd for the planning session; see notes/overnight/RUNBOOK.md).

Usage:  python3 -I tools/overnight/meter.py [--before ISO_TIMESTAMP] [--json]
Prints the estimated USD spent by every transcript in this container.
"""
import glob, json, os, sys

# USD per million tokens. Input/output/cache-read from the claude-api skill (cached 2026-10-06).
# Cache writes assumed 1.25x (5 min TTL) and 2x (1 h TTL) of the input price.
PRICES = {
    "claude-fable-5-1":  dict(inp=10.0, out=50.0, cr=0.25, cw5=12.5, cw1h=20.0),
    "claude-opus-5-5":   dict(inp=4.0,  out=20.0, cr=0.20, cw5=5.0,  cw1h=8.0),
    "claude-sonnet-5-5": dict(inp=2.0,  out=10.0, cr=0.20, cw5=2.5,  cw1h=4.0),
    "claude-haiku-5-5":  dict(inp=0.10, out=0.50, cr=0.01, cw5=0.125, cw1h=0.20),
}
FALLBACK = PRICES["claude-fable-5-1"]   # unknown model: price as the most expensive one
SAFETY = 1.15


def estimate(before=None):
    seen = {}
    pattern = os.path.expanduser("~/.claude/projects/**/*.jsonl")
    for path in glob.glob(pattern, recursive=True):
        try:
            fh = open(path, encoding="utf-8", errors="replace")
        except OSError:
            continue
        with fh:
            for line in fh:
                try:
                    d = json.loads(line)
                except ValueError:
                    continue
                if d.get("type") != "assistant":
                    continue
                m = d.get("message")
                if not isinstance(m, dict) or not m.get("usage"):
                    continue
                ts = d.get("timestamp") or ""
                if before and ts and ts >= before:
                    continue
                seen[(path, m.get("id"), d.get("requestId"))] = (m.get("model"), m["usage"])
    by_model, total = {}, 0.0
    for model, u in seen.values():
        p = PRICES.get(model, FALLBACK)
        cc = u.get("cache_creation") or {}
        cw1h = cc.get("ephemeral_1h_input_tokens")
        cw5 = cc.get("ephemeral_5m_input_tokens")
        cw_all = u.get("cache_creation_input_tokens") or 0
        if cw1h is None and cw5 is None:          # no breakdown: assume the dearer 1 h TTL
            cw1h, cw5 = cw_all, 0
        cw1h, cw5 = cw1h or 0, cw5 or 0
        usd = ((u.get("input_tokens") or 0) * p["inp"] + (u.get("output_tokens") or 0) * p["out"]
               + (u.get("cache_read_input_tokens") or 0) * p["cr"] + cw5 * p["cw5"] + cw1h * p["cw1h"]) / 1e6
        b = by_model.setdefault(model, dict(requests=0, usd=0.0, out=0, cache_read=0))
        b["requests"] += 1; b["usd"] += usd
        b["out"] += u.get("output_tokens") or 0; b["cache_read"] += u.get("cache_read_input_tokens") or 0
        total += usd
    return dict(usd_raw=round(total, 3), usd_est=round(total * SAFETY, 3), safety=SAFETY,
                requests=len(seen), by_model={k: {**v, "usd": round(v["usd"], 3)} for k, v in by_model.items()})


if __name__ == "__main__":
    before = sys.argv[sys.argv.index("--before") + 1] if "--before" in sys.argv else None
    r = estimate(before)
    if "--json" in sys.argv:
        print(json.dumps(r))
    else:
        print(f"estimated spend in this container: ${r['usd_est']:.2f} "
              f"(raw ${r['usd_raw']:.2f} x safety {r['safety']}, {r['requests']} requests)")
