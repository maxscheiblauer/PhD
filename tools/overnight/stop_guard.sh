#!/usr/bin/env bash
# Stop hook for overnight worker sessions. Inert unless /tmp/overnight_worker exists
# (created by the worker at the start of its block, see notes/overnight/PROTOCOL.md).
# While the block is open, it refuses to let the worker end its turn and tells it to continue.
# Safety valves: the block is closed (/tmp/overnight_block_closed), or 60 refusals in total.
[ -f /tmp/overnight_worker ] || exit 0
[ -f /tmp/overnight_block_closed ] && exit 0
n=$(( $(cat /tmp/overnight_stop_refusals 2>/dev/null || echo 0) + 1 ))
echo "$n" > /tmp/overnight_stop_refusals
[ "$n" -gt 60 ] && exit 0
here=$(cd "$(dirname "$0")" && pwd)
spent=$(python3 -I "$here/meter.py" 2>/dev/null | sed -n 's/.*\$\([0-9.]*\) (raw.*/\1/p')
cap=$(cat /tmp/overnight_block_cap 2>/dev/null || echo "?")
reason="Overnight worker protocol: your block is still open (meter: \$${spent:-?} of block cap \$${cap}). Do not stop and do not wait for the user: nobody is watching. If the meter is below 85% of the cap and the wall-clock limit is not reached, continue with the next item in notes/overnight/STATE.md. Otherwise run the end-of-block procedure in notes/overnight/PROTOCOL.md section 6 (it creates /tmp/overnight_block_closed); then you may stop."
python3 -I -c 'import json,sys; print(json.dumps({"decision":"block","reason":sys.argv[1]}))' "$reason"
