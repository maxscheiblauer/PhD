#!/bin/bash
# usage: bash tools/overnight/hb.sh <status>   -- rewrites heartbeat.json with current meter
cd "$(dirname "$0")/../.."
M=$(python3 -I tools/overnight/meter.py --json | python3 -I -c 'import json,sys;print(round(json.load(sys.stdin)["usd_est"],2))')
printf '{"block": 11, "session": "session_01Fi3wVb8YpPmehxM4suTyC6", "status": "%s", "meter_usd": %s, "utc": "%s"}\n' "$1" "$M" "$(date -u +%FT%TZ)" > notes/overnight/heartbeat.json
echo "meter_usd=$M"
