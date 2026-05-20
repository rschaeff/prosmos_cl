#!/bin/bash
# Submit a ProSMoS searchmatrix SLURM array for every *.query in QDIR.
# Usage:
#   QDIR=path/to/queries OUT=path/to/results [CONCURRENCY=64] ./submit.sh
#
# Required env:
#   QDIR  - directory of *.query files
#   OUT   - NFS-visible output root (will be created)
# Optional env:
#   DB           - metamatricesDB path (default: 2010 DB)
#   CONCURRENCY  - max concurrent array tasks (default: 64). Cap this on a
#                  busy cluster; per leda etiquette do not run all 1488 at once.
#   TIMELIMIT    - per-task SLURM --time. Default 01:00:00 from array.sbatch;
#                  override e.g. TIMELIMIT=04:00:00 for searches that exceed
#                  the default (we saw 5-of-9 design-target searches finish
#                  in 18-51 min on the 306MB 2010 DB; 4 of 9 needed >1h).

set -euo pipefail

: "${QDIR:?missing QDIR}"
: "${OUT:?missing OUT}"
: "${DB:=/home/rschaeff/src/Prosmos/ProSMoS/metamatrixdb/metamatricesDB}"
: "${CONCURRENCY:=64}"
: "${TIMELIMIT:=}"

[ -d "$QDIR" ] || { echo "QDIR not a directory: $QDIR" >&2; exit 1; }
[ -f "$DB" ]   || { echo "DB not found: $DB" >&2; exit 1; }

mkdir -p "$OUT/logs" "$OUT/parts" "$OUT/work" "$OUT/hits"

# Build the query manifest (one path per line, sorted for reproducibility).
QUERY_LIST="$OUT/queries.list"
find "$QDIR" -maxdepth 1 -name '*.query' -type f | sort > "$QUERY_LIST"
N=$(wc -l < "$QUERY_LIST")
[ "$N" -gt 0 ] || { echo "no *.query files in $QDIR" >&2; exit 1; }

HERE=$(cd "$(dirname "$0")" && pwd)

echo "submitting array of $N tasks (concurrency cap %${CONCURRENCY}${TIMELIMIT:+, --time=$TIMELIMIT})" >&2
TIME_ARG=()
[ -n "$TIMELIMIT" ] && TIME_ARG+=(--time="$TIMELIMIT")
ARRAY_ID=$(sbatch --parsable \
    "${TIME_ARG[@]}" \
    --array=1-${N}%${CONCURRENCY} \
    --chdir="$OUT/logs" \
    --export=QUERY_LIST="$QUERY_LIST",DB="$DB",OUT="$OUT" \
    "$HERE/array.sbatch")
echo "array job id: $ARRAY_ID" >&2

# Merge per-task TSV parts into one summary after the array finishes.
MERGE_ID=$(sbatch --parsable \
    --job-name=prosmos-merge \
    --time=00:05:00 \
    --mem=512M \
    --dependency=afterany:${ARRAY_ID} \
    --chdir="$OUT/logs" \
    --wrap="{ printf 'query\truntime_sec\texit_code\thits\n'; cat $OUT/parts/*.tsv 2>/dev/null | sort; } > $OUT/summary.tsv && echo 'summary at $OUT/summary.tsv'")
echo "merge job id:  $MERGE_ID (afterany)" >&2

echo "watch with: squeue -j ${ARRAY_ID} ; tail -f ${OUT}/logs/prosmos-search.${ARRAY_ID}_*.err" >&2
