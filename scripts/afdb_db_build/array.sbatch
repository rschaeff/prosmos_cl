#!/bin/bash
#SBATCH --job-name=afdb-db-build
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err
#
# Per-array-task entry point for the AFDB metamatricesDB build. Expects env
# from submit.sh:
#   CHUNK_DIR   - directory of chunk_NNN.txt files (NFS, read-only)
#   FRAG_DIR    - per-task fragment outputs (NFS, written once at end)
#   SCRIPT_DIR  - directory holding process_chunk.sh
#
# Work dir is forced to node-local /tmp (the AFDB build does 25k+ tiny
# I/Os per chunk; running that on NFS murders performance).

set -euo pipefail

: "${CHUNK_DIR:?missing CHUNK_DIR}"
: "${FRAG_DIR:?missing FRAG_DIR}"
: "${SCRIPT_DIR:?missing SCRIPT_DIR}"
: "${PARALLEL:=8}"

TASK_ID="$SLURM_ARRAY_TASK_ID"
CHUNK="$CHUNK_DIR/chunk_${TASK_ID}.txt"
TASK_WORK="/tmp/afdb_build_${SLURM_JOB_ID}_${TASK_ID}"
FRAG="$FRAG_DIR/fragment_${TASK_ID}.frag"

echo "[$(date)] Task $TASK_ID on $(hostname) -- chunk=$CHUNK work=$TASK_WORK"

mkdir -p "$TASK_WORK"
PARALLEL="$PARALLEL" "$SCRIPT_DIR/process_chunk.sh" "$CHUNK" "$TASK_WORK" "$FRAG"
# process_chunk.sh cleans up its subdirs; this is the safety net.
rm -rf "$TASK_WORK"

echo "[$(date)] Task $TASK_ID finished"
