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
# sbatch's --export whitelist strips PATH, so `python` resolves to
# /usr/bin/python which lacks gemmi and other PALSSE deps. Use an absolute
# path to the conda python (override PYTHON if PALSSE lives elsewhere).
: "${PYTHON:=/sw/apps/Anaconda3-2023.09-0/bin/python}"
# Per-domain timeouts. PALSSE can hang on huge or pathological PDBs (e.g.
# ribosomes); some chunks in the previous run timed out at 2h because one
# single PDB never returned. generateMatrix also hangs occasionally.
: "${PALSSE_TIMEOUT:=60s}"
: "${GENMAT_TIMEOUT:=30s}"

[ -f "$CHUNK" ] || { echo "chunk file not found: $CHUNK" >&2; exit 1; }
[ -x "$GENMAT" ] || { echo "generateMatrix not executable: $GENMAT" >&2; exit 1; }
[ -d "$PALSSE_REPO" ] || { echo "PALSSE repo not found: $PALSSE_REPO" >&2; exit 1; }
[ -x "$PYTHON" ] || { echo "PYTHON not executable: $PYTHON" >&2; exit 1; }

mkdir -p "$WORK/ssd" "$WORK/per_file"

N=$(wc -l < "$CHUNK")
echo "Processing $N PDB paths with PARALLEL=$PARALLEL"

# Stage 1: PALSSE each PDB. Per-domain failures are tolerated (|| true) so a
# handful of bad PDBs don't kill the whole chunk under `set -e` — the final
# .out count tells us how many actually succeeded. `timeout` bounds runaway
# domains (one stuck PDB previously soaked an entire 2h chunk).
echo "[$(date +%H:%M:%S)] PALSSE stage starting (per-domain timeout=$PALSSE_TIMEOUT)"
time xargs -P "$PARALLEL" -I {} bash -c '
  pdb="{}"
  uid=$(basename "$pdb" .pdb)
  timeout "'"$PALSSE_TIMEOUT"'" \
    env PYTHONPATH="'"$PALSSE_REPO"'" "'"$PYTHON"'" -m palsse.cli "$pdb" -o "'"$WORK"'/ssd/${uid}.ssd" 2>/dev/null || true
' < "$CHUNK"
SSD_COUNT=$(find "$WORK/ssd" -name "*.ssd" | wc -l)
echo "[$(date +%H:%M:%S)] PALSSE: $SSD_COUNT .ssd files produced"

# Stage 2: generateMatrix per .ssd file.
# generateMatrix is OpenMPI-linked. Under a SLURM step it sees SLURM_*/PMI_*
# in the env and tries to bootstrap via PMI, which fails on this cluster
# (OpenMPI built without --with-pmi). Strip ALL SLURM/PMI/PMIX/OMPI vars via
# env-glob so MPI_Init falls back to standalone singleton — a whitelist was
# too brittle (different nodes set different PMIX_* variants).
SLURM_PMI_UNSET=$(env | awk -F= '/^(SLURM|PMI|PMIX|OMPI)_/{printf " -u %s", $1}')
echo "[$(date +%H:%M:%S)] generateMatrix stage starting (per-domain timeout=$GENMAT_TIMEOUT)"
time find "$WORK/ssd" -name "*.ssd" -printf "%f\n" | xargs -P "$PARALLEL" -I {} bash -c '
  ssd="{}"
  uid="${ssd%.ssd}"
  timeout "'"$GENMAT_TIMEOUT"'" env '"$SLURM_PMI_UNSET"' \
    "'"$GENMAT"'" -os "$ssd" "'"$WORK"'/ssd/" "'"$WORK"'/per_file/${uid}.out" > /dev/null 2>&1 || true
'
OUT_COUNT=$(find "$WORK/per_file" -name "*.out" | wc -l)
echo "[$(date +%H:%M:%S)] generateMatrix: $OUT_COUNT .out files produced"

# Stage 3: concatenate into per-chunk fragment
echo "[$(date +%H:%M:%S)] Concatenating to $OUT"
# Sort filenames so concatenation order is reproducible across runs.
find "$WORK/per_file" -name "*.out" | sort | xargs cat > "$OUT"
echo "[$(date +%H:%M:%S)] Fragment size: $(stat -c %s "$OUT") bytes"

echo "[$(date +%H:%M:%S)] Done."
