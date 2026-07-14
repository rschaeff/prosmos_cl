#!/bin/bash
# Fold a finished sweep's hits/ tree into one compressed table (see archive_hits.py).
# Run this only after the sweep AND its retry pass are complete -- archiving a live
# tree would miss the hits still being written.
#
# Usage:
#   OUT=path/to/run [SHARDS=64] [CONCURRENCY=32] ./submit_archive.sh
set -euo pipefail
: "${OUT:?missing OUT}"
: "${SHARDS:=64}"
: "${CONCURRENCY:=32}"
HERE=$(cd "$(dirname "$0")" && pwd)
PY=/sw/apps/Anaconda3-2023.09-0/bin/python

HITS="$OUT/hits"
ARC="$OUT/archive"
[ -d "$HITS" ] || { echo "no $HITS" >&2; exit 1; }
mkdir -p "$ARC" "$OUT/logs"

JID=$(sbatch --parsable --job-name=hits-archive --time=08:00:00 --mem=4G \
    --array=0-$((SHARDS-1))%${CONCURRENCY} --chdir="$OUT/logs" \
    --wrap="$PY $HERE/archive_hits.py worker $HITS $ARC \$SLURM_ARRAY_TASK_ID $SHARDS")
echo "archive array: $JID ($SHARDS shards)" >&2

FID=$(sbatch --parsable --job-name=hits-archive-fin --time=04:00:00 --mem=16G \
    --dependency=afterok:"$JID" --chdir="$OUT/logs" \
    --wrap="$PY $HERE/archive_hits.py finalize $HITS $ARC $SHARDS")
echo "finalize:      $FID" >&2
echo >&2
echo "when $FID prints VERIFIED, the tree can be removed with:" >&2
echo "  rm -rf $HITS" >&2
