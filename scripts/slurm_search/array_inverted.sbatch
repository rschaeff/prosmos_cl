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

# Validate the manifest BEFORE searching. searchmatrix treats an unreadable query
# path as "no query" and exits 0, so a malformed manifest yields a full sweep that
# reports zero hits and looks like a scientific result. (This happened: a manifest
# built with `ls` -- aliased to long format -- put "-rw-r--r-- 1 rschaeff lab ..."
# on every line, so all 11 templates silently vanished and 800 chunks searched an
# empty query set.) Fail loudly on a bad path instead.
NQ=0
while read -r q; do
    [ -n "$q" ] || continue
    if [ ! -f "$q" ]; then
        echo "manifest $MANIFEST: not a file: '$q' (built with \`ls\`? use find)" >&2; exit 3
    fi
    NQ=$((NQ + 1))
done < "$MANIFEST"
if [ "$NQ" = "0" ]; then
    echo "manifest $MANIFEST has no queries" >&2; exit 3
fi

# Node-local scratch: searchmatrix writes <cwd>/../sheetbug and a DB-named log.
LOCAL="/tmp/smi_${LINE}_$$"
mkdir -p "$LOCAL/run" "$LOCAL/sheetbug"
mkdir -p "$OUT/parts" "$OUT/logs"
cd "$LOCAL/run"

# HITS_MODE=local: searchmatrix writes its one-file-per-hit output to the node's
# own disk, and we fold it into a single compressed TSV on NFS before exiting.
# The tree-mode alternative has every task writing ~30k tiny files straight into
# 6336 NFS directories shared with every other task -- the full AFDB sweep left
# ~24M inodes behind that way, and the write storm is in searchmatrix's hot loop.
# One TSV per chunk means ~800 files per sweep instead of tens of millions.
# Default is tree, so runs already in flight (and their retries) keep the layout
# they started with; submit_inverted.sh opts new runs into local.
: "${HITS_MODE:=tree}"
if [ "$HITS_MODE" = "local" ]; then
    HR="$LOCAL/hits"
    # searchmatrix does not create the per-query dirs, it only writes into them.
    mkdir -p "$HR"
    while read -r q; do mkdir -p "$HR/$(basename "$q" .query)"; done < "$MANIFEST"
    HR="$HR/"
else
    HR="$HITS_ROOT"
fi

# Chunks are small (~a few MB); reading from NFS is fine, no staging needed.
START=$(date +%s)
"$SEARCH" "$MANIFEST" "$CHUNK" "$HR" > "$LOCAL/stdout.log" 2>&1
RC=$?
END=$(date +%s)
ELAPSED=$((END - START))

# Fold this chunk's node-local hits into one zstd TSV on NFS. Records live in
# exactly one chunk, so the per-chunk TSVs concatenate with no dedup needed.
if [ "$HITS_MODE" = "local" ] && [ "$RC" = "0" ]; then
    mkdir -p "$OUT/hitparts"
    "${PYTHON:-/sw/apps/Anaconda3-2023.09-0/bin/python}" \
        "${ARCHIVE_PY:-/home/rschaeff/dev/prosmos_cl/scripts/slurm_search/archive_hits.py}" \
        chunk "$LOCAL/hits" "$OUT/hitparts/${NAME}.tsv.zst" >> "$LOCAL/stdout.log" 2>&1
    RC=$?
fi

# Hits produced by THIS chunk (files whose mtime is >= START would be exact, but
# a simple count of this chunk's records + rc is the health signal).
NREC=$(grep -c '\.ssd' "$CHUNK" 2>/dev/null)
# retry_missing.sh sets PARTS_PREFIX so a retry array's line numbers (1..n_missing)
# cannot overwrite the original run's parts, which are keyed by global chunk line.
printf '%s\t%d\t%d\t%d\n' "$NAME" "$ELAPSED" "$RC" "$NREC" > "$OUT/parts/${PARTS_PREFIX:+${PARTS_PREFIX}_}${LINE}.tsv"
tail -50 "$LOCAL/stdout.log" > "$OUT/logs/${NAME}.tail" 2>/dev/null || true
rm -rf "$LOCAL"
echo "[$(date +%H:%M:%S)] chunk ${NAME}: ${ELAPSED}s rc=${RC} recs=${NREC}" >&2
