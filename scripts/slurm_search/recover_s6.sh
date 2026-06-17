#!/bin/bash
# Throttled recovery launcher for the S6 sweep. Submits 1000-task chunks
# one at a time, waiting for the user job count to drop below a threshold
# before each submission. MaxJobCount=10000 on this cluster, so chaining all
# 142 remaining chunks via afterany at once exceeded the per-user limit.
#
# Usage: bash recover_s6.sh START_OFFSET
#   START_OFFSET — first chunk's offset (e.g. 10000 to resume after 10
#                   completed chunks of 1000)

set -uo pipefail

START_OFFSET="${1:?missing START_OFFSET}"
OUT=$HOME/work/prosmos_2026/ecod_search_v4_s6
QLIST=$OUT/queries.list
DB=$HOME/work/prosmos_2026/ecod_db_manual_v4/metamatricesDB
HERE=$(cd "$(dirname "$0")" && pwd)
N=151808
CHUNK=1000
CONC=128
HEADROOM=5000  # pause if queue has more than this many jobs

EXCLUDE=$(sinfo -h -N -o '%N' 2>&1 | sort -u | awk '/^leda[0-9]+$/ {n=substr($0,5); if (n+0 < 100) print}' | paste -sd,)

# Find the most recent already-submitted job we can chain off, so the first
# new chunk waits for it. Pull from squeue; fall back to none (will start
# immediately).
PREV=$(squeue -u "$USER" -h -t PD,R -o '%i' 2>/dev/null | awk -F_ '{print $1}' | sort -u | tail -1)
[ -z "$PREV" ] && echo "no previous job to chain off" >&2

ARRAY_IDS=()
offset=$START_OFFSET
while [ "$offset" -lt "$N" ]; do
    # Throttle: wait for the queue to drop below HEADROOM
    while true; do
        nq=$(squeue -u "$USER" -h | wc -l)
        if [ "$nq" -lt "$HEADROOM" ]; then break; fi
        echo "[$(date +%H:%M:%S)] queue=$nq >= $HEADROOM, waiting..."
        sleep 60
    done

    remaining=$((N - offset))
    this_chunk=$((remaining < CHUNK ? remaining : CHUNK))
    DEP_ARG=()
    [ -n "$PREV" ] && DEP_ARG+=(--dependency=afterany:"$PREV")
    for attempt in 1 2 3 4 5; do
        JOB_ID=$(sbatch --parsable \
            --exclude="$EXCLUDE" \
            "${DEP_ARG[@]}" \
            --array=1-${this_chunk}%${CONC} \
            --chdir="$OUT/logs" \
            --export=QUERY_LIST="$QLIST",DB="$DB",OUT="$OUT",OFFSET="$offset" \
            "$HERE/array.sbatch" 2>/dev/null)
        [ -n "$JOB_ID" ] && break
        echo "[$(date +%H:%M:%S)] sbatch failed (attempt $attempt), retrying in 30s..." >&2
        sleep 30
    done
    if [ -z "$JOB_ID" ]; then
        echo "[$(date +%H:%M:%S)] FAILED offset=$offset, aborting" >&2
        exit 1
    fi
    ARRAY_IDS+=("$JOB_ID")
    echo "[$(date +%H:%M:%S)] ok offset=$offset id=$JOB_ID (queue was $nq)"
    PREV=$JOB_ID
    offset=$((offset + this_chunk))
    sleep 2  # small per-sbatch sleep regardless
done

# Submit final merge.
DEPS=$(IFS=:; echo "${ARRAY_IDS[*]}")
MERGE_ID=$(sbatch --parsable \
    --job-name=prosmos-merge \
    --time=00:05:00 \
    --mem=512M \
    --dependency=afterany:"$DEPS" \
    --chdir="$OUT/logs" \
    --wrap="{ printf 'query\truntime_sec\texit_code\thits\n'; cat $OUT/parts/*.tsv 2>/dev/null | sort; } > $OUT/summary.tsv && echo 'summary at $OUT/summary.tsv'")
echo "[$(date +%H:%M:%S)] merge job id: $MERGE_ID"
echo "[$(date +%H:%M:%S)] submitted ${#ARRAY_IDS[@]} chunks + 1 merge"
