#!/bin/bash
#SBATCH --job-name=instr-agree
#SBATCH --time=02:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=8
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --partition=192GB,All,96GB
#
# Sharded instrument-agreement search. Each array task owns ~48 pairs; it builds
# a small exp DB and a small pred DB from ITS shard and searches both against the
# full 6,336-query manifest. Both halves of a pair are in the same shard so the
# task is self-contained, and both halves run on one node so binary/hardware/env
# are held fixed -- structure provenance is the only variable, which is the point.
set -uo pipefail

R=/home/rschaeff/work/prosmos_2026/instr_agree
MANIFEST=/home/rschaeff/work/prosmos_2026/s5_pdb_inv/manifest.txt
SEARCH=/home/rschaeff/dev/prosmos_cl/searchMatrix/build/searchmatrix
PROC=/home/rschaeff/dev/prosmos_cl/scripts/ecod_db_build/process_chunk.sh
ARCHIVE=/home/rschaeff/dev/prosmos_cl/scripts/slurm_search/archive_hits.py
PY=/sw/apps/Anaconda3-2023.09-0/bin/python
I=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")

# Manifest guard: searchmatrix treats an unreadable query path as "no query" and
# exits 0, so a bad manifest would yield a clean-looking sweep reporting zero hits.
NQ=0
while read -r q; do [ -n "$q" ] || continue; [ -f "$q" ] || { echo "bad query: $q" >&2; exit 3; }; NQ=$((NQ+1)); done < "$MANIFEST"
[ "$NQ" = "6336" ] || { echo "expected 6336 queries, got $NQ" >&2; exit 3; }

LOCAL=/tmp/instr_${SLURM_ARRAY_JOB_ID}_${I}
mkdir -p "$LOCAL" "$R/hitparts" "$R/dbs"

for SIDE in exp pred; do
    SHARD="$R/shards/${SIDE}_${I}.txt"
    NIN=$(wc -l < "$SHARD")
    echo "=== [$(date +%H:%M:%S)] task $I $SIDE: $NIN input structures"
    PARALLEL=8 bash "$PROC" "$SHARD" "$LOCAL/work_$SIDE" "$R/dbs/${SIDE}_${I}.db"
    NREC=$(grep -c '\.ssd' "$R/dbs/${SIDE}_${I}.db")
    echo "=== [$(date +%H:%M:%S)] task $I $SIDE: DB has $NREC records"
    # A shard yielding almost nothing means PALSSE/genmat failed wholesale on this
    # node; searching it would fabricate a low-agreement result. Demand most of them.
    [ "$NREC" -ge $((NIN / 2)) ] || { echo "$SIDE shard too small ($NREC/$NIN)" >&2; exit 4; }

    H="$LOCAL/hits_$SIDE"; mkdir -p "$H"
    while read -r q; do mkdir -p "$H/$(basename "$q" .query)"; done < "$MANIFEST"
    mkdir -p "$LOCAL/run_$SIDE/sheetbug" "$LOCAL/run_$SIDE/cwd"; cd "$LOCAL/run_$SIDE/cwd"
    "$SEARCH" "$MANIFEST" "$R/dbs/${SIDE}_${I}.db" "$H/" > "$LOCAL/${SIDE}.search.log" 2>&1
    echo "=== [$(date +%H:%M:%S)] task $I $SIDE: search rc=$? -- archiving"
    "$PY" "$ARCHIVE" chunk "$H" "$R/hitparts/${SIDE}_${I}.tsv.zst"
done

rm -rf "$LOCAL"
echo "[$(date +%H:%M:%S)] task $I done"
