#!/bin/bash
#SBATCH --job-name=prosmos-search
#SBATCH --time=01:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
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

# Pick this task's query from the manifest
QUERY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$QUERY_LIST")
if [ -z "$QUERY" ] || [ ! -f "$QUERY" ]; then
    echo "no query at line ${SLURM_ARRAY_TASK_ID} of $QUERY_LIST" >&2
    exit 2
fi
NAME=$(basename "$QUERY" .query)

# Per-task work area. searchmatrix unconditionally writes to <cwd>/../sheetbug/
# and silently exit(0)s if it can't open that file, so each task needs its own
# writable parent. Using a per-task subdir avoids cross-task collisions.
WORK="$OUT/work/${SLURM_ARRAY_TASK_ID}"
HITS="$OUT/hits/${NAME}"
mkdir -p "$WORK/run" "$HITS" "$OUT/parts" "$OUT/logs"
cd "$WORK/run"

START=$(date +%s)
"$SEARCH" "$QUERY" "$DB" "$HITS" > "$OUT/logs/${NAME}.stdout" 2>&1
RC=$?
END=$(date +%s)
ELAPSED=$((END - START))

# Hit count is the only reliable success signal — searchmatrix exits 0 even on
# parse errors and on its silent sheetbug-failure exit path.
HIT_COUNT=$(find "$HITS" -maxdepth 1 -type f | wc -l)

printf '%s\t%d\t%d\t%d\n' "$NAME" "$ELAPSED" "$RC" "$HIT_COUNT" \
    > "$OUT/parts/${SLURM_ARRAY_TASK_ID}.tsv"

# Clean up the sheetbug debug output (large, useless after the run)
rm -rf "$WORK/sheetbug"

echo "[$(date +%H:%M:%S)] ${NAME}: ${ELAPSED}s rc=${RC} hits=${HIT_COUNT}" >&2
