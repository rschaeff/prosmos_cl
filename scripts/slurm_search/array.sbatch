#!/bin/bash
#SBATCH --job-name=prosmos-search
#SBATCH --time=01:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
# Default to all partitions the 2G/1cpu task fits, so chained chunks spread
# across the cluster instead of starving on the busy 96GB default. submit.sh's
# PARTITION env still overrides this (command-line --partition wins).
#SBATCH --partition=192GB,All,96GB
#
# ProSMoS searchmatrix array job. One query per array task. Submit via
# scripts/slurm_search/submit.sh (which sets --array and --export).
#
# Required env (must be passed by submit.sh):
#   QUERY_LIST  - file with one absolute query.txt path per line
#   DB          - absolute path to metamatricesDB
#   OUT         - output root (NFS-visible, writable)
# Optional env:
#   SEARCH      - searchmatrix binary (default: repo's build/searchmatrix)

set -uo pipefail

: "${QUERY_LIST:?missing QUERY_LIST}"
: "${DB:?missing DB}"
: "${OUT:?missing OUT}"
: "${SEARCH:=/home/rschaeff/dev/prosmos_cl/searchMatrix/build/searchmatrix}"

# Pick this task's query from the manifest. submit.sh splits the manifest
# into chunks (SLURM's MaxArraySize is ~1001) and passes the chunk's base
# OFFSET; this task's global line is SLURM_ARRAY_TASK_ID + OFFSET.
: "${OFFSET:=0}"
LINE=$((SLURM_ARRAY_TASK_ID + OFFSET))
QUERY=$(sed -n "${LINE}p" "$QUERY_LIST")
if [ -z "$QUERY" ] || [ ! -f "$QUERY" ]; then
    echo "no query at line ${LINE} of $QUERY_LIST" >&2
    exit 2
fi
NAME=$(basename "$QUERY" .query)

# Per-task work area: use the compute node's local /tmp for the run dir AND
# for the stdout log. The earlier all-NFS layout strangled searchmatrix on
# I/O contention (logs are O(GB) of verbose per-entry chatter -- the v3 smoke
# took 60+ min crawling on the head node before we moved logs to local /tmp,
# then ran in 18-30 min). searchmatrix also writes to <cwd>/../sheetbug/ and
# silently exit(0)s if it can't, so the local-/tmp parent must have that dir.
LOCAL="/tmp/sm_${LINE}_$$"
mkdir -p "$LOCAL/run" "$LOCAL/sheetbug"
HITS="$OUT/hits/${NAME}"
mkdir -p "$HITS" "$OUT/parts" "$OUT/logs"
cd "$LOCAL/run"

START=$(date +%s)
"$SEARCH" "$QUERY" "$DB" "$HITS" > "$LOCAL/stdout.log" 2>&1
RC=$?
END=$(date +%s)
ELAPSED=$((END - START))

# Hit count is the only reliable success signal -- searchmatrix exits 0 even on
# parse errors and on its silent sheetbug-failure exit path.
HIT_COUNT=$(find "$HITS" -maxdepth 1 -type f | wc -l)

printf '%s\t%d\t%d\t%d\n' "$NAME" "$ELAPSED" "$RC" "$HIT_COUNT" \
    > "$OUT/parts/${LINE}.tsv"

# Save the last 200 lines of the searchmatrix log for forensics, drop the rest.
tail -200 "$LOCAL/stdout.log" > "$OUT/logs/${NAME}.tail" 2>/dev/null || true

# Clean local scratch (avoid filling /tmp on compute nodes with many tasks).
rm -rf "$LOCAL"

echo "[$(date +%H:%M:%S)] ${NAME}: ${ELAPSED}s rc=${RC} hits=${HIT_COUNT}" >&2
