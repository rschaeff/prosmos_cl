#!/bin/bash
#SBATCH --job-name=prosmos-inv
#SBATCH --time=08:00:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=1
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#SBATCH --partition=192GB,All,96GB
#
# Loop-inverted ProSMoS search: ONE task per DB CHUNK, each task runs ALL queries
# (a manifest) against its chunk in a single DB parse. The DB is thus parsed once
# total (across chunks) instead of once per query -- ~2x less compute than the
# query-per-task model. Submit via submit_inverted.sh.
#
# Required env (from submit_inverted.sh):
#   CHUNK_LIST - file with one absolute DB-chunk path per line
#   MANIFEST   - file with one absolute query.txt path per line (all queries)
#   HITS_ROOT  - output root; hits land in HITS_ROOT/<queryname>/pdb<id>.txt
#   OUT        - run root (for parts/logs)
# Optional:
#   SEARCH     - searchmatrix binary (default: repo build)

set -uo pipefail
: "${CHUNK_LIST:?}"; : "${MANIFEST:?}"; : "${HITS_ROOT:?}"; : "${OUT:?}"
: "${SEARCH:=/home/rschaeff/dev/prosmos_cl/searchMatrix/build/searchmatrix}"
: "${OFFSET:=0}"

LINE=$((SLURM_ARRAY_TASK_ID + OFFSET))
CHUNK=$(sed -n "${LINE}p" "$CHUNK_LIST")
if [ -z "$CHUNK" ] || [ ! -f "$CHUNK" ]; then
    echo "no chunk at line ${LINE} of $CHUNK_LIST" >&2; exit 2
fi
NAME=$(basename "$CHUNK")

# Node-local scratch: searchmatrix writes <cwd>/../sheetbug and a DB-named log.
LOCAL="/tmp/smi_${LINE}_$$"
mkdir -p "$LOCAL/run" "$LOCAL/sheetbug"
mkdir -p "$OUT/parts" "$OUT/logs"
cd "$LOCAL/run"

# Chunks are small (~a few MB); reading from NFS is fine, no staging needed.
START=$(date +%s)
"$SEARCH" "$MANIFEST" "$CHUNK" "$HITS_ROOT" > "$LOCAL/stdout.log" 2>&1
RC=$?
END=$(date +%s)
ELAPSED=$((END - START))

# Hits produced by THIS chunk (files whose mtime is >= START would be exact, but
# a simple count of this chunk's records + rc is the health signal).
NREC=$(grep -c '\.ssd' "$CHUNK" 2>/dev/null)
printf '%s\t%d\t%d\t%d\n' "$NAME" "$ELAPSED" "$RC" "$NREC" > "$OUT/parts/${LINE}.tsv"
tail -50 "$LOCAL/stdout.log" > "$OUT/logs/${NAME}.tail" 2>/dev/null || true
rm -rf "$LOCAL"
echo "[$(date +%H:%M:%S)] chunk ${NAME}: ${ELAPSED}s rc=${RC} recs=${NREC}" >&2
