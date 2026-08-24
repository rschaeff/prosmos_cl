#!/bin/bash
# Build the AFDB metamatricesDB from grey's pre-staged non-singleton
# structure set at ~grey/afdb.200m/non_singleton_4p9m_structures/.
#
# Usage:
#   OUT=$HOME/work/prosmos_2026/afdb_db [N_CHUNKS=210] [PARALLEL=8] \
#     [TIMELIMIT=02:00:00] [CONCURRENCY=64] ./submit.sh
#
# Manifest construction is local (no DB pull) -- we enumerate the three
# AFDB source subtrees directly:
#   - non_singleton_4p9m_structures/new/<shard>/<XX>.tar.gz  (~91k tarballs)
#   - non_singleton_4p9m_structures/dpam/<shard>/<XX>.tar.gz (~7k tarballs)
#   - non_singleton_4p9m_structures/ecod/ecod_*.pdb.gz       (~7.7k flat gz)
#
# Work units = paths. process_chunk.sh routes by extension.

set -euo pipefail

: "${OUT:?missing OUT}"
: "${N_CHUNKS:=210}"
: "${PARALLEL:=8}"
: "${TIMELIMIT:=02:00:00}"
: "${CONCURRENCY:=64}"
: "${AFDB_ROOT:=/home/grey/afdb.200m/non_singleton_4p9m_structures}"

HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT/manifest" "$OUT/chunks" "$OUT/fragments" "$OUT/logs"

# Stage 1: build the manifest of work units.
MANIFEST="$OUT/manifest/paths.txt"
echo "[$(date)] Enumerating AFDB work units from $AFDB_ROOT..."
{
    find "$AFDB_ROOT/new" -name '*.tar.gz' -type f
    find "$AFDB_ROOT/dpam" -name '*.tar.gz' -type f
    find "$AFDB_ROOT/ecod" -name '*.pdb.gz' -type f
} | sort > "$MANIFEST"
N_PATHS=$(wc -l < "$MANIFEST")
echo "[$(date)] Manifest: $N_PATHS work units"
[ "$N_PATHS" -gt 0 ] || { echo "empty manifest -- check AFDB_ROOT" >&2; exit 1; }

# Stage 2: chunk it.
CHUNK_SIZE=$(( (N_PATHS + N_CHUNKS - 1) / N_CHUNKS ))
echo "[$(date)] Splitting $N_PATHS units into $N_CHUNKS chunks of $CHUNK_SIZE"
find "$OUT/chunks" -name 'chunk_*.txt' -delete
awk -v cs="$CHUNK_SIZE" -v dir="$OUT/chunks" '
  { idx = int((NR - 1) / cs) + 1
    out = dir "/chunk_" idx ".txt"
    print >> out }
' "$MANIFEST"
ACTUAL=$(find "$OUT/chunks" -name 'chunk_*.txt' | wc -l)
if [ "$ACTUAL" -ne "$N_CHUNKS" ]; then
    echo "[$(date)] WARN: produced $ACTUAL chunks (expected $N_CHUNKS)" >&2
    N_CHUNKS="$ACTUAL"
fi

# Stage 3: submit the array. Same bad-MPI exclude pattern as the ECOD build.
: "${EXCLUDE_NODES:=leda23,leda25,leda27,leda29,leda31,leda33,leda35,leda37,leda39,leda41,leda43,leda45}"
echo "[$(date)] Submitting array of $N_CHUNKS tasks (concurrency $CONCURRENCY, --time $TIMELIMIT, --cpus-per-task $PARALLEL, --exclude=$EXCLUDE_NODES)"
ARRAY_ID=$(sbatch --parsable \
    --time="$TIMELIMIT" \
    --cpus-per-task="$PARALLEL" \
    --array="1-${N_CHUNKS}%${CONCURRENCY}" \
    --exclude="$EXCLUDE_NODES" \
    --chdir="$OUT/logs" \
    --export=CHUNK_DIR="$OUT/chunks",FRAG_DIR="$OUT/fragments",SCRIPT_DIR="$HERE",PARALLEL="$PARALLEL" \
    "$HERE/array.sbatch")
echo "[$(date)] Array job id: $ARRAY_ID"

# Stage 4: dependent merge.
MERGE_ID=$(sbatch --parsable \
    --job-name=afdb-db-merge \
    --time=00:30:00 \
    --mem=2G \
    --dependency=afterany:"$ARRAY_ID" \
    --chdir="$OUT/logs" \
    --wrap="bash -c '
set -euo pipefail
echo \"[\$(date)] Merging fragments...\"
find \"$OUT/fragments\" -name fragment_\\*.frag | sort -V | xargs cat > \"$OUT/metamatricesDB\"
echo \"[\$(date)] metamatricesDB: \$(stat -c %s \"$OUT/metamatricesDB\") bytes\"
echo \"[\$(date)] Entries: \$(grep -c \\.ssd \"$OUT/metamatricesDB\")\"
'")
echo "[$(date)] Merge job id: $MERGE_ID (afterany)"

echo
echo "Watch: squeue -j $ARRAY_ID,$MERGE_ID"
echo "Final DB: $OUT/metamatricesDB"
echo "Remember to post-filter for malformed entries (>4-digit-coord rows)"
echo "before the negspace sweep -- see leda v3 .clean pattern."
