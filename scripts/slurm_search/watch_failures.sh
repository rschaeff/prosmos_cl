#!/bin/bash
# Scan a slurm_search OUT directory's parts/*.tsv for nonzero exit codes and
# print a summary + per-failure detail (query name, wall time, hit count, log
# tail). One-shot; loop externally with `watch` or a `while sleep`.
#
# Usage:
#   ./watch_failures.sh [OUT_DIR]      # default: ~/work/ecod_search_v3
#
# Exits 0 if no failures (or no parts yet), 1 if any rc!=0 are present.

set -u

OUT="${1:-$HOME/work/ecod_search_v3}"
PARTS="$OUT/parts"
LOGS="$OUT/logs"

[ -d "$PARTS" ] || { echo "no parts dir: $PARTS" >&2; exit 0; }

# parts/*.tsv format: <query>\t<runtime_sec>\t<exit_code>\t<hits>
TOTAL=$(find "$PARTS" -maxdepth 1 -name '*.tsv' -type f | wc -l)
[ "$TOTAL" -eq 0 ] && { echo "[$(date +%H:%M:%S)] 0 parts written yet"; exit 0; }

# Tally rc distribution. awk reads all tsvs at once.
SUMMARY=$(awk -F'\t' '
    { rc[$3]++ }
    END { for (r in rc) printf "  rc=%s: %d\n", r, rc[r] }
' "$PARTS"/*.tsv | sort)

FAILS=$(awk -F'\t' '$3 != 0 { print }' "$PARTS"/*.tsv)
NFAIL=$(printf '%s' "$FAILS" | grep -c . || true)

printf '[%s] %d parts written, %d nonzero-rc\n' "$(date +%H:%M:%S)" "$TOTAL" "$NFAIL"
printf '%s\n' "$SUMMARY"

if [ "$NFAIL" -eq 0 ]; then
    exit 0
fi

printf '\n=== failures ===\n'
printf 'query\truntime_s\trc\thits\n'
printf '%s\n' "$FAILS" | sort

# For each failure, show the last few lines of the searchmatrix log if it was
# preserved (array.sbatch keeps tail -200 at $OUT/logs/<query>.tail).
printf '\n=== last 5 log lines per failure ===\n'
while IFS=$'\t' read -r q rt rc hits; do
    [ -z "$q" ] && continue
    printf '\n--- %s (rc=%s, %ss, %s hits) ---\n' "$q" "$rc" "$rt" "$hits"
    if [ -f "$LOGS/${q}.tail" ]; then
        tail -5 "$LOGS/${q}.tail"
    else
        echo "(no log tail at $LOGS/${q}.tail)"
    fi
done <<< "$FAILS"

exit 1
