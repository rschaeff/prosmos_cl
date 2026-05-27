#!/bin/bash
#SBATCH --job-name=ecod-db-build
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err

# Per-array-task entry point. Expects environment from submit.sh:
#   CHUNK_DIR  — directory of chunk_NNN.txt files
#   WORK_ROOT  — root for per-task work dirs (NFS-visible)
#   FRAG_DIR   — per-task fragment outputs
#   SCRIPT_DIR — directory holding process_chunk.sh (SLURM stages this script
#                into /var/spool/slurmd/jobNNN/ so $(dirname "$0") is useless)
#
# Each task processes chunk_${SLURM_ARRAY_TASK_ID}.txt and writes
# fragment_${SLURM_ARRAY_TASK_ID}.frag.

set -euo pipefail

: "${CHUNK_DIR:?missing CHUNK_DIR}"
: "${WORK_ROOT:?missing WORK_ROOT}"
: "${FRAG_DIR:?missing FRAG_DIR}"
: "${SCRIPT_DIR:?missing SCRIPT_DIR}"
: "${PARALLEL:=8}"

TASK_ID="$SLURM_ARRAY_TASK_ID"
CHUNK="$CHUNK_DIR/chunk_${TASK_ID}.txt"
TASK_WORK="$WORK_ROOT/task_${TASK_ID}"
FRAG="$FRAG_DIR/fragment_${TASK_ID}.frag"

echo "[$(date)] Task $TASK_ID starting on $(hostname); chunk=$CHUNK; work=$TASK_WORK"

PARALLEL="$PARALLEL" "$SCRIPT_DIR/process_chunk.sh" "$CHUNK" "$TASK_WORK" "$FRAG"

echo "[$(date)] Task $TASK_ID finished"

# Clean up intermediate .ssd / .out files (keep only the fragment)
rm -rf "$TASK_WORK/ssd" "$TASK_WORK/per_file"
rmdir "$TASK_WORK" 2>/dev/null || true
