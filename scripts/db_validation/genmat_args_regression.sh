#!/bin/bash
# generateMatrix front-end regression: the argument handling in generMatrix.cpp
# may change, the records it writes may not. Runs a BASELINE binary (e.g. the
# v1.0.0 build the census DBs were made with) and a CANDIDATE over every .ssd
# in a directory and requires byte-identical output, for the legacy call
# (-os <name> <dir/> <out>) and for each new spelling the candidate accepts.
#
# Usage:
#   genmat_args_regression.sh <baseline_bin> <candidate_bin> <ssd_dir> [parallel]
# Exit 0 iff every file matches in every spelling. Work files go to a mktemp
# dir. Each call costs ~2 s of MPI start-up and there are five per file, so
# shard large corpora over a SLURM array: NSHARD=8 SHARD=<0..7> takes every
# 8th file, e.g.
#   sbatch --array=0-7 -c 16 --wrap='NSHARD=8 SHARD=$SLURM_ARRAY_TASK_ID \
#     genmat_args_regression.sh <base> <cand> <ssd_dir> 16'
# (5,088 AFDB files took ~65 min unsharded on 16 cores.)
set -u
BASE="${1:?baseline binary}"; CAND="${2:?candidate binary}"; SSD="${3:?ssd dir}"
PAR="${4:-8}"
SSD=$(cd "$SSD" && pwd)
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
# generateMatrix is an OpenMPI singleton; SLURM/PMI variables confuse it.
UNSET=$(env | awk -F= '/^(SLURM|PMI|PMIX|OMPI)_/{printf " -u %s", $1}')
export BASE CAND SSD W UNSET

NSHARD="${NSHARD:-1}"; SHARD="${SHARD:-0}"
find "$SSD" -maxdepth 1 -name '*.ssd' -printf '%f\n' | sort \
    | awk -v n="$NSHARD" -v k="$SHARD" '(NR - 1) % n == k' > "$W/list"
N=$(wc -l < "$W/list")
[ "$N" -gt 0 ] || { echo "no .ssd files in $SSD" >&2; exit 2; }

xargs -P "$PAR" -I{} bash -c '
  f={}; u=${f%.ssd}; d="$W/$u"; mkdir -p "$d"
  env $UNSET "$BASE" -os "$f" "$SSD/" "$d/base" >/dev/null 2>&1
  env $UNSET "$CAND" -os "$f" "$SSD/" "$d/legacy"  >/dev/null 2>&1
  env $UNSET "$CAND" -os "$f" "$SSD"  "$d/noslash" >/dev/null 2>&1
  env $UNSET "$CAND" "$SSD/$f" "$d/path"           >/dev/null 2>&1
  env $UNSET "$CAND" -os "$SSD/$f" "$d/ospath"     >/dev/null 2>&1
  for v in legacy noslash path ospath; do
    cmp -s "$d/base" "$d/$v" || echo "DIFF $u $v"
  done
  [ -s "$d/base" ] || echo "EMPTYBASE $u"
  rm -rf "$d"
' < "$W/list" > "$W/result"

NDIFF=$(grep -c '^DIFF' "$W/result")
NEMPTY=$(grep -c '^EMPTYBASE' "$W/result")
echo "shard $SHARD/$NSHARD  files: $N   differing (file x spelling): $NDIFF   baseline wrote nothing: $NEMPTY"
grep '^DIFF' "$W/result" | head -20
[ "$NDIFF" -eq 0 ]
