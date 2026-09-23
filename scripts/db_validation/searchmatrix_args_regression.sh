#!/bin/bash
# searchmatrix front-end regression: the argument handling may change, the hits
# may not. Runs a BASELINE binary (e.g. the v1.0.0 build the census sweeps used)
# and a CANDIDATE over a query set and a DB, and requires byte-identical hit
# trees for
#   - manifest mode (all queries, one DB pass), as the SLURM sweeps run it;
#   - single-query mode for every query, in the legacy spelling (absolute
#     paths, output dir with trailing slash) and, for the candidate only,
#     relative paths and an output dir without the slash.
#
# Usage:
#   searchmatrix_args_regression.sh <baseline_bin> <candidate_bin> <db> <query_dir> [parallel]
# Every *.query under <query_dir> (recursively) is used. Exit 0 iff all trees
# match. Work files go to a mktemp dir on local disk. Cost is 2 manifest passes
# plus 3 single-query passes per query over the whole DB, and low-dimension
# queries (S2/S3) hit most records: use a DB of a few thousand to ~20k records.
set -u
BASE="${1:?baseline binary}"; CAND="${2:?candidate binary}"
DB=$(readlink -f "${3:?db}"); QDIR=$(readlink -f "${4:?query dir}"); PAR="${5:-8}"
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
export BASE CAND DB W

find "$QDIR" -name '*.query' -type f | sort > "$W/manifest"
NQ=$(wc -l < "$W/manifest")
[ "$NQ" -gt 0 ] || { echo "no .query files under $QDIR" >&2; exit 2; }
FAIL=0

# Both binaries write ../sheetbug relative to their CWD: give each its own.
mkdir -p "$W/b/run" "$W/c/run"
( cd "$W/b/run" && "$BASE" "$W/manifest" "$DB" "$W/b/man/" >/dev/null 2>&1 ) &
( cd "$W/c/run" && "$CAND" "$W/manifest" "$DB" "$W/c/man"  >/dev/null 2>&1 ) &
wait
if diff -rq "$W/b/man" "$W/c/man" > "$W/man.diff"; then
    echo "manifest mode: identical ($NQ queries, $(find "$W/b/man" -type f | wc -l) hit files)"
else
    echo "manifest mode: DIFFERS"; head -20 "$W/man.diff"; FAIL=1
fi

# Single-query mode, one query per task.
xargs -P "$PAR" -I{} bash -c '
  q={}; n=$(basename "$q" .query); d="$W/s/$n"; mkdir -p "$d/b/run" "$d/c/run" "$d/r/run"
  (cd "$d/b/run" && "$BASE" "$q" "$DB" "$d/b/hits/" >/dev/null 2>&1)
  (cd "$d/c/run" && "$CAND" "$q" "$DB" "$d/c/hits/" >/dev/null 2>&1)
  cp "$q" "$d/r/run/q.query"; ln -s "$DB" "$d/r/run/db"
  (cd "$d/r/run" && "$CAND" q.query db hits >/dev/null 2>&1)
  diff -rq "$d/b/hits" "$d/c/hits"      >/dev/null 2>&1 || echo "DIFF $n legacy"
  diff -rq "$d/b/hits" "$d/r/run/hits"  >/dev/null 2>&1 || echo "DIFF $n relative"
  echo "HITS $n $(find "$d/b/hits" -type f | wc -l)"
  rm -rf "$d"
' < "$W/manifest" > "$W/single"
NDIFF=$(grep -c '^DIFF' "$W/single")
NHIT=$(awk '/^HITS/{s+=$3} END{print s+0}' "$W/single")
echo "single-query mode: $NQ queries x 2 spellings, $NHIT baseline hit files, $NDIFF differing"
grep '^DIFF' "$W/single" | head -20
[ "$NDIFF" -eq 0 ] || FAIL=1
exit $FAIL
