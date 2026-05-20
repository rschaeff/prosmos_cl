#!/bin/bash
# Process one chunk of ECOD domain PDB paths through PALSSE + generateMatrix
# and emit a per-chunk metamatricesDB fragment.
#
# Usage:
#   process_chunk.sh <chunk_paths.txt> <work_dir> <out.frag>
#
# Where:
#   chunk_paths.txt  — one PDB path per line (e.g. /data/ecod/.../12345.pdb)
#   work_dir         — scratch dir for .ssd / .out files (NFS-visible)
#   out.frag         — final per-chunk concatenated matrices file
#
# Steps:
#   1. PALSSE each PDB → .ssd in work_dir/ssd/
#   2. generateMatrix each .ssd → .out in work_dir/per_file/
#   3. cat all .out → out.frag

set -euo pipefail

CHUNK="${1:?chunk_paths.txt required}"
WORK="${2:?work_dir required}"
OUT="${3:?out.frag required}"

: "${PARALLEL:=8}"
: "${PALSSE_REPO:=/home/rschaeff/dev/palsse_cl}"
: "${GENMAT:=/home/rschaeff/dev/prosmos_cl/generateMatrix/build/generateMatrix}"

[ -f "$CHUNK" ] || { echo "chunk file not found: $CHUNK" >&2; exit 1; }
[ -x "$GENMAT" ] || { echo "generateMatrix not executable: $GENMAT" >&2; exit 1; }
[ -d "$PALSSE_REPO" ] || { echo "PALSSE repo not found: $PALSSE_REPO" >&2; exit 1; }

mkdir -p "$WORK/ssd" "$WORK/per_file"

N=$(wc -l < "$CHUNK")
echo "Processing $N PDB paths with PARALLEL=$PARALLEL"

# Stage 1: PALSSE each PDB
echo "[$(date +%H:%M:%S)] PALSSE stage starting"
time xargs -P "$PARALLEL" -I {} bash -c '
  pdb="{}"
  uid=$(basename "$pdb" .pdb)
  PYTHONPATH="'"$PALSSE_REPO"'" python -m palsse.cli "$pdb" -o "'"$WORK"'/ssd/${uid}.ssd" 2>/dev/null
' < "$CHUNK"
SSD_COUNT=$(find "$WORK/ssd" -name "*.ssd" | wc -l)
echo "[$(date +%H:%M:%S)] PALSSE: $SSD_COUNT .ssd files produced"

# Stage 2: generateMatrix per .ssd file
echo "[$(date +%H:%M:%S)] generateMatrix stage starting"
time find "$WORK/ssd" -name "*.ssd" -printf "%f\n" | xargs -P "$PARALLEL" -I {} bash -c '
  ssd="{}"
  uid="${ssd%.ssd}"
  "'"$GENMAT"'" -os "$ssd" "'"$WORK"'/ssd/" "'"$WORK"'/per_file/${uid}.out" > /dev/null 2>&1
'
OUT_COUNT=$(find "$WORK/per_file" -name "*.out" | wc -l)
echo "[$(date +%H:%M:%S)] generateMatrix: $OUT_COUNT .out files produced"

# Stage 3: concatenate into per-chunk fragment
echo "[$(date +%H:%M:%S)] Concatenating to $OUT"
# Sort filenames so concatenation order is reproducible across runs.
find "$WORK/per_file" -name "*.out" | sort | xargs cat > "$OUT"
echo "[$(date +%H:%M:%S)] Fragment size: $(stat -c %s "$OUT") bytes"

echo "[$(date +%H:%M:%S)] Done."
