#!/bin/bash
# Orchestrate ECOD metamatricesDB build as a SLURM array.
#
# Pulls a manifest of complete PDB-source ECOD domain paths from
# ecod_commons.derived_files, splits into N chunks, submits one array
# task per chunk, then a dependent merge job that concatenates all
# per-chunk fragments into the final metamatricesDB.
#
# Usage:
#   OUT=$HOME/work/ecod_db [N_CHUNKS=200] [PARALLEL=8] [TIMELIMIT=02:00:00] \
#     [DB_FILTER="status='complete' AND domain_source_type='pdb'"] \
#     ./submit.sh
#
# Required env:
#   OUT          — output root (NFS-visible)
#
# Optional env:
#   N_CHUNKS     — number of array tasks (default 200)
#   PARALLEL     — concurrent workers per task (default 8)
#   TIMELIMIT    — per-task --time (default 02:00:00)
#   DB_FILTER    — SQL WHERE clause for the manifest query. Default selects
#                  complete PDB-source domains. Override e.g. for F40:
#                    DB_FILTER="status='complete' AND domain_source_type='pdb'
#                               AND ecod_uid IN (SELECT ... F40)"
#   CONCURRENCY  — max concurrent array tasks (default 64)
#   PGHOST/PGPORT/PGUSER/PGDATABASE  — DB connection (defaults to dione/45000)

set -euo pipefail

: "${OUT:?missing OUT}"
: "${N_CHUNKS:=200}"
: "${PARALLEL:=8}"
: "${TIMELIMIT:=02:00:00}"
: "${DB_FILTER:=file_type_id=2 AND status='complete' AND domain_source_type='pdb'}"
: "${CONCURRENCY:=64}"

: "${PGHOST:=dione}"
: "${PGPORT:=45000}"
: "${PGUSER:=ecod}"
: "${PGDATABASE:=ecod_protein}"
: "${PGPASSWORD:=***REMOVED***}"
export PGHOST PGPORT PGUSER PGDATABASE PGPASSWORD

HERE=$(cd "$(dirname "$0")" && pwd)

mkdir -p "$OUT/manifest" "$OUT/chunks" "$OUT/work" "$OUT/fragments" "$OUT/logs"

# Stage 1: pull manifest of paths from the DB
MANIFEST="$OUT/manifest/paths.txt"
echo "[$(date)] Pulling manifest from ecod_commons.derived_files..."
psql -At -c "
  SELECT internal_path FROM ecod_commons.derived_files
  WHERE ${DB_FILTER}
  ORDER BY ecod_uid;
" > "$MANIFEST"
N_PATHS=$(wc -l < "$MANIFEST")
echo "[$(date)] Manifest: $N_PATHS paths"
[ "$N_PATHS" -gt 0 ] || { echo "empty manifest — refusing to submit" >&2; exit 1; }

# Stage 2: split into chunks. Use ceiling division so the last chunk
# may be smaller; numbered 1..N_CHUNKS to match SLURM array indices.
CHUNK_SIZE=$(( (N_PATHS + N_CHUNKS - 1) / N_CHUNKS ))
echo "[$(date)] Splitting $N_PATHS paths into $N_CHUNKS chunks of $CHUNK_SIZE"
# `split` would name files alphabetically; use awk to write numerically.
awk -v cs="$CHUNK_SIZE" -v dir="$OUT/chunks" '
  { idx = int((NR - 1) / cs) + 1
    out = dir "/chunk_" idx ".txt"
    print >> out }
' "$MANIFEST"

# Verify we got N_CHUNKS chunk files (the last may be missing if N_PATHS
# was an exact multiple of N_CHUNKS — handle by checking the last file).
ACTUAL_CHUNKS=$(find "$OUT/chunks" -name "chunk_*.txt" | wc -l)
if [ "$ACTUAL_CHUNKS" -ne "$N_CHUNKS" ]; then
    echo "[$(date)] WARNING: produced $ACTUAL_CHUNKS chunks (expected $N_CHUNKS)" >&2
    N_CHUNKS="$ACTUAL_CHUNKS"
fi

# Stage 3: submit the array
echo "[$(date)] Submitting array of $N_CHUNKS tasks (concurrency $CONCURRENCY, --time $TIMELIMIT, --cpus-per-task $PARALLEL)"
ARRAY_ID=$(sbatch --parsable \
    --time="$TIMELIMIT" \
    --cpus-per-task="$PARALLEL" \
    --array="1-${N_CHUNKS}%${CONCURRENCY}" \
    --chdir="$OUT/logs" \
    --export=CHUNK_DIR="$OUT/chunks",WORK_ROOT="$OUT/work",FRAG_DIR="$OUT/fragments",PARALLEL="$PARALLEL" \
    "$HERE/array.sbatch")
echo "[$(date)] Array job id: $ARRAY_ID"

# Stage 4: dependent merge — cat all fragments into the final DB.
# `--dependency=afterany` so partial fragments still get assembled if
# some tasks fail (the operator can re-run only the failed array indices).
MERGE_ID=$(sbatch --parsable \
    --job-name=ecod-db-merge \
    --time=00:30:00 \
    --mem=2G \
    --dependency=afterany:"$ARRAY_ID" \
    --chdir="$OUT/logs" \
    --wrap="
set -euo pipefail
echo \"[\$(date)] Merging fragments...\"
find \"$OUT/fragments\" -name 'fragment_*.frag' | sort -V | xargs cat > \"$OUT/metamatricesDB\"
echo \"[\$(date)] metamatricesDB: \$(stat -c %s \"$OUT/metamatricesDB\") bytes\"
echo \"[\$(date)] Entries: \$(grep -c '\.ssd' \"$OUT/metamatricesDB\")\"
")
echo "[$(date)] Merge job id: $MERGE_ID (afterany)"

echo
echo "Watch with:"
echo "  squeue -j $ARRAY_ID,$MERGE_ID"
echo "  tail -f $OUT/logs/ecod-db-build.${ARRAY_ID}_1.err"
echo
echo "Final DB will land at: $OUT/metamatricesDB"
