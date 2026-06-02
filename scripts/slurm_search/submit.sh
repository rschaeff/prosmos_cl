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
# Recurses into subdirectories so callers can point QDIR at a parent like
# `example/ssp_enumerated/queries_typed/` that holds s3/, s4/, s5/ subdirs.
QUERY_LIST="$OUT/queries.list"
find "$QDIR" -name '*.query' -type f | sort > "$QUERY_LIST"
N=$(wc -l < "$QUERY_LIST")
[ "$N" -gt 0 ] || { echo "no *.query files in $QDIR" >&2; exit 1; }

HERE=$(cd "$(dirname "$0")" && pwd)

# SLURM caps array size (MaxArraySize, typically 1001 on this cluster). If we
# have more queries than that, split into chunks. Each chunk gets its own
# array job, with OFFSET set so its task IDs map back into the global manifest.
MAX_ARRAY=$(scontrol show config 2>/dev/null | awk '/MaxArraySize/{print $3}')
: "${MAX_ARRAY:=1000}"
# Use chunk size = MaxArraySize - 1 to leave headroom; SLURM rejects 1-N where
# N >= MaxArraySize.
CHUNK=$((MAX_ARRAY - 1))

echo "submitting $N tasks in chunks of <=${CHUNK} (concurrency %${CONCURRENCY}${TIMELIMIT:+, --time=$TIMELIMIT})" >&2
TIME_ARG=()
[ -n "$TIMELIMIT" ] && TIME_ARG+=(--time="$TIMELIMIT")
EXCLUDE_ARG=()
[ -n "${EXCLUDE_NODES:-}" ] && EXCLUDE_ARG+=(--exclude="$EXCLUDE_NODES")

# Chain the chunks: each chunk waits for the previous one to finish so the
# %${CONCURRENCY} cap is the TOTAL concurrent task count (not per-chunk ×
# n_chunks). On the leda cluster (~50 nodes) running 8 chunks × 64 concurrent
# would mean 512 simultaneous tasks, which would hog the whole cluster and
# starve other users.
ARRAY_IDS=()
PREV=""
offset=0
while [ "$offset" -lt "$N" ]; do
    remaining=$((N - offset))
    this_chunk=$((remaining < CHUNK ? remaining : CHUNK))
    DEP_ARG=()
    [ -n "$PREV" ] && DEP_ARG+=(--dependency=afterany:"$PREV")
    JOB_ID=$(sbatch --parsable \
        "${TIME_ARG[@]}" \
        "${EXCLUDE_ARG[@]}" \
        "${DEP_ARG[@]}" \
        --array=1-${this_chunk}%${CONCURRENCY} \
        --chdir="$OUT/logs" \
        --export=QUERY_LIST="$QUERY_LIST",DB="$DB",OUT="$OUT",OFFSET="$offset" \
        "$HERE/array.sbatch")
    ARRAY_IDS+=("$JOB_ID")
    echo "  array job id: $JOB_ID  (offset $offset, tasks 1..$this_chunk${PREV:+, afterany:$PREV})" >&2
    PREV="$JOB_ID"
    offset=$((offset + this_chunk))
done

# Merge per-task TSV parts after ALL array jobs finish.
DEPS=$(IFS=:; echo "${ARRAY_IDS[*]}")
MERGE_ID=$(sbatch --parsable \
    --job-name=prosmos-merge \
    --time=00:05:00 \
    --mem=512M \
    --dependency=afterany:${DEPS} \
    --chdir="$OUT/logs" \
    --wrap="{ printf 'query\truntime_sec\texit_code\thits\n'; cat $OUT/parts/*.tsv 2>/dev/null | sort; } > $OUT/summary.tsv && echo 'summary at $OUT/summary.tsv'")
echo "merge job id: $MERGE_ID (afterany on ${#ARRAY_IDS[@]} array job(s))" >&2

echo "watch with: squeue -j ${ARRAY_IDS[*]}" >&2
