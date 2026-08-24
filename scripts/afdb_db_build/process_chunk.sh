#!/bin/bash
# Process one chunk of AFDB work units (tar.gz or .pdb.gz files) through
# PALSSE + generateMatrix, emit a per-chunk metamatricesDB fragment.
#
# Usage:
#   process_chunk.sh <chunk_paths.txt> <work_dir> <out.frag>
#
# Where:
#   chunk_paths.txt — one path per line. Each path is either:
#     - a .tar.gz tarball containing bare .pdb files (AFDB new/, dpam/), or
#     - a .pdb.gz file (AFDB ecod/)
#   work_dir        — scratch dir. SHOULD be node-local /tmp; the caller
#                     (array.sbatch) is responsible for picking a local path.
#                     NFS work dirs murder this script: 25k+ small files /chunk.
#   out.frag        — final per-chunk concatenated matrices file (NFS OK,
#                     it's one big write at the end).
#
# Adapted from ecod_db_build/process_chunk.sh. The main difference: input
# units are bundled (tarballs of ~50 PDBs each + flat gzipped PDBs) rather
# than already-unpacked bare PDB paths, so we add an extraction stage and
# route by extension.

set -uo pipefail   # NOT -e: per-PDB failures are tolerated

CHUNK="${1:?chunk_paths.txt required}"
WORK="${2:?work_dir required}"
OUT="${3:?out.frag required}"

: "${PARALLEL:=8}"
: "${PALSSE_REPO:=/home/rschaeff/dev/palsse_cl}"
: "${GENMAT:=/home/rschaeff/dev/prosmos_cl/generateMatrix/build/generateMatrix}"
: "${PYTHON:=/sw/apps/Anaconda3-2023.09-0/bin/python}"
: "${PALSSE_TIMEOUT:=60s}"
: "${GENMAT_TIMEOUT:=30s}"

[ -f "$CHUNK" ] || { echo "chunk file not found: $CHUNK" >&2; exit 1; }
[ -x "$GENMAT" ] || { echo "generateMatrix not executable: $GENMAT" >&2; exit 1; }
[ -d "$PALSSE_REPO" ] || { echo "PALSSE repo not found: $PALSSE_REPO" >&2; exit 1; }
[ -x "$PYTHON" ] || { echo "PYTHON not executable: $PYTHON" >&2; exit 1; }

mkdir -p "$WORK/pdb" "$WORK/ssd" "$WORK/per_file"

# Strip SLURM/PMI env so generateMatrix's OpenMPI singleton fallback works.
SLURM_PMI_UNSET=$(env | awk -F= '/^(SLURM|PMI|PMIX|OMPI)_/{printf " -u %s", $1}')

N_UNITS=$(wc -l < "$CHUNK")
echo "[$(date +%H:%M:%S)] Chunk has $N_UNITS work units"

# ---------------- Stage 0: extract all PDBs from work units ---------------
# For each tarball: extract all .pdb inside to $WORK/pdb/. For each flat
# .pdb.gz: gunzip to $WORK/pdb/. Sequential is fine here -- extraction is
# I/O bound and parallelism doesn't help on local /tmp.
echo "[$(date +%H:%M:%S)] Extracting PDBs..."
EXTRACT_FAILED=0
while IFS= read -r path; do
    [ -z "$path" ] && continue
    case "$path" in
        *.tar.gz)
            # Extract bare .pdb files; tarballs contain flat layout (no subdirs).
            tar -C "$WORK/pdb" -xzf "$path" 2>/dev/null || EXTRACT_FAILED=$((EXTRACT_FAILED+1))
            ;;
        *.pdb.gz)
            base=$(basename "$path" .gz)
            gunzip -c "$path" > "$WORK/pdb/$base" 2>/dev/null || EXTRACT_FAILED=$((EXTRACT_FAILED+1))
            ;;
        *)
            echo "  unknown extension: $path" >&2
            EXTRACT_FAILED=$((EXTRACT_FAILED+1))
            ;;
    esac
done < "$CHUNK"
N_PDB=$(find "$WORK/pdb" -maxdepth 1 -name '*.pdb' -type f | wc -l)
echo "[$(date +%H:%M:%S)] Extracted $N_PDB PDBs ($EXTRACT_FAILED units failed)"
[ "$N_PDB" -eq 0 ] && { echo "no PDBs to process" >&2; touch "$OUT"; exit 0; }

# ---------------- Stage 1: PALSSE -> .ssd ---------------------------------
echo "[$(date +%H:%M:%S)] PALSSE stage starting (PARALLEL=$PARALLEL, timeout=$PALSSE_TIMEOUT)"
time find "$WORK/pdb" -maxdepth 1 -name '*.pdb' -type f -printf '%f\n' | xargs -P "$PARALLEL" -I {} bash -c '
  pdb="{}"
  uid="${pdb%.pdb}"
  timeout "'"$PALSSE_TIMEOUT"'" \
    env PYTHONPATH="'"$PALSSE_REPO"'" "'"$PYTHON"'" -m palsse.cli "'"$WORK"'/pdb/$pdb" -o "'"$WORK"'/ssd/${uid}.ssd" 2>/dev/null || true
'
SSD_COUNT=$(find "$WORK/ssd" -name '*.ssd' | wc -l)
echo "[$(date +%H:%M:%S)] PALSSE: $SSD_COUNT .ssd files produced"

# ---------------- Stage 2: generateMatrix -> .out -------------------------
echo "[$(date +%H:%M:%S)] generateMatrix stage starting (PARALLEL=$PARALLEL, timeout=$GENMAT_TIMEOUT)"
time find "$WORK/ssd" -name '*.ssd' -printf '%f\n' | xargs -P "$PARALLEL" -I {} bash -c '
  ssd="{}"
  uid="${ssd%.ssd}"
  timeout "'"$GENMAT_TIMEOUT"'" env '"$SLURM_PMI_UNSET"' \
    "'"$GENMAT"'" -os "$ssd" "'"$WORK"'/ssd/" "'"$WORK"'/per_file/${uid}.out" > /dev/null 2>&1 || true
'
OUT_COUNT=$(find "$WORK/per_file" -name '*.out' | wc -l)
echo "[$(date +%H:%M:%S)] generateMatrix: $OUT_COUNT .out files produced"

# ---------------- Stage 3: concatenate ------------------------------------
echo "[$(date +%H:%M:%S)] Concatenating to $OUT"
find "$WORK/per_file" -name '*.out' | sort | xargs cat > "$OUT"
echo "[$(date +%H:%M:%S)] Fragment size: $(stat -c %s "$OUT") bytes; ${N_PDB} input PDBs -> ${OUT_COUNT} matrices"

# Stage 4: free local /tmp for the next task on this node
rm -rf "$WORK/pdb" "$WORK/ssd" "$WORK/per_file"
echo "[$(date +%H:%M:%S)] Done."
