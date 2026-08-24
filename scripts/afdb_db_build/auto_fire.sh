#!/bin/bash
# Post-build orchestrator. Runs AFTER the AFDB metamatricesDB merge job
# completes (--dependency=afterok:<merge_id>). Does, in sequence:
#
#   1. Filter the merged DB for malformed entries (>4-digit-coord lines)
#      -> metamatricesDB.clean. Pattern verified against leda v3 .clean
#      (same 707567 entries kept of 707651 total).
#
#   2. Sanity-check the .clean DB against two known-positive queries
#      (FHB s4-0026-0000, RLM s5-0098-0021). If either returns 0 hits,
#      something is broken (DB structure mismatch, encoding bug) -- abort
#      before launching the expensive sweep.
#
#   3. If sanity passes, submit the 800-query negspace sweep against the
#      .clean DB via scripts/slurm_search/submit.sh.
#
# Usage (typically invoked by sbatch --wrap; can also be run directly):
#   bash auto_fire.sh <DB_BUILD_OUT> <SWEEP_OUT>
#
# Exits 0 on success, nonzero if filter or sanity fails.

set -euo pipefail

DB_OUT="${1:?missing DB_BUILD_OUT (e.g. \$HOME/work/prosmos_2026/afdb_db)}"
SWEEP_OUT="${2:?missing SWEEP_OUT (e.g. \$HOME/work/prosmos_2026/afdb_negspace_sweep)}"

PROSMOS=/home/rschaeff/dev/prosmos_cl
SEARCH=$PROSMOS/searchMatrix/build/searchmatrix
QDIR=$PROSMOS/enum/negspace_queries
SANITY_FHB=$PROSMOS/example/ssp_enumerated/queries_typed/s4/s4-0026-0000.query
SANITY_RLM=$PROSMOS/example/ssp_enumerated/queries_typed/s5/s5-0098-0021.query

DB=$DB_OUT/metamatricesDB
CLEAN=$DB_OUT/metamatricesDB.clean

# ---------------- Stage 1: filter ----------------------------------------
echo "[$(date)] Stage 1: filtering malformed entries"
[ -s "$DB" ] || { echo "$DB missing or empty -- did the build merge fail?" >&2; exit 2; }
RAW_COUNT=$(grep -c '\.ssd' "$DB")
awk '!/[0-9][0-9][0-9][0-9]\.[0-9]/' "$DB" > "$CLEAN"
CLEAN_COUNT=$(grep -c '\.ssd' "$CLEAN")
echo "[$(date)] Filter: $RAW_COUNT -> $CLEAN_COUNT entries (dropped $((RAW_COUNT - CLEAN_COUNT)))"
[ "$CLEAN_COUNT" -gt 0 ] || { echo "filter produced empty DB -- pattern mismatch?" >&2; exit 3; }

# ---------------- Stage 2: sanity ----------------------------------------
echo "[$(date)] Stage 2: sanity controls"
SANITY_DIR=/tmp/afdb_sanity_$$
rm -rf "$SANITY_DIR"
mkdir -p "$SANITY_DIR/fhb_hits" "$SANITY_DIR/rlm_hits" "$SANITY_DIR/sheetbug" "$SANITY_DIR/run"
cd "$SANITY_DIR/run"

run_sanity() {
    local name="$1" query="$2" hitdir="$3"
    [ -f "$query" ] || { echo "sanity query missing: $query" >&2; return 1; }
    echo "[$(date)] running $name ($query) ..."
    local start=$(date +%s)
    timeout 1800 "$SEARCH" "$query" "$CLEAN" "$hitdir" > /dev/null 2>&1 || true
    local elapsed=$(($(date +%s) - start))
    local hits=$(find "$hitdir" -maxdepth 1 -type f | wc -l)
    echo "[$(date)] $name: $hits hits in ${elapsed}s"
    if [ "$hits" -lt 100 ]; then
        echo "$name returned only $hits hits -- expected thousands; aborting" >&2
        return 1
    fi
}

run_sanity FHB "$SANITY_FHB" "$SANITY_DIR/fhb_hits" || { rm -rf "$SANITY_DIR"; exit 4; }
run_sanity RLM "$SANITY_RLM" "$SANITY_DIR/rlm_hits" || { rm -rf "$SANITY_DIR"; exit 5; }
rm -rf "$SANITY_DIR"
echo "[$(date)] Sanity controls passed."

# ---------------- Stage 3: launch sweep ----------------------------------
echo "[$(date)] Stage 3: launching 800-query negspace sweep"
cd "$PROSMOS"
SWEEP_LOG=$(mktemp)
QDIR="$QDIR" \
DB="$CLEAN" \
OUT="$SWEEP_OUT" \
CONCURRENCY=128 \
TIMELIMIT=12:00:00 \
bash scripts/slurm_search/submit.sh 2>&1 | tee "$SWEEP_LOG"

# Parse the sweep-merge job id from submit.sh stdout.
# submit.sh prints e.g. "merge job id: 553204 (afterany on 1 array job(s))".
SWEEP_MERGE_ID=$(awk '/merge job id:/{print $4; exit}' "$SWEEP_LOG")
rm -f "$SWEEP_LOG"
if [ -z "$SWEEP_MERGE_ID" ]; then
    echo "[$(date)] WARN: could not parse sweep merge job id; analysis will need manual launch" >&2
    exit 0
fi
echo "[$(date)] sweep merge job id: $SWEEP_MERGE_ID"

# ---------------- Stage 4: chain analysis after sweep merge --------------
ANALYSIS_OUT=$SWEEP_OUT/analysis
mkdir -p "$ANALYSIS_OUT"
ANALYSIS_ID=$(sbatch --parsable \
    --job-name=afdb-analyze \
    --time=00:30:00 \
    --mem=2G \
    --dependency=afterok:$SWEEP_MERGE_ID \
    --chdir=$SWEEP_OUT \
    --output=analysis.%j.out \
    --error=analysis.%j.err \
    --wrap="bash $PROSMOS/scripts/afdb_db_build/analyze_quadrants.sh '$SWEEP_OUT' '$ANALYSIS_OUT'")
echo "[$(date)] analysis job id: $ANALYSIS_ID (afterok:$SWEEP_MERGE_ID)"
echo "[$(date)] auto_fire complete -- chain: sweep array -> sweep merge -> analysis -> $ANALYSIS_OUT/"
