#!/bin/bash
# Final 4-quadrant analysis combining ALL PDB + AFDB sweep results:
#
#   PDB data sources:
#     ~/work/prosmos_2026/ecod_search_v4_negspace/parts/        (leda partial 1)
#     ~/work/prosmos_2026/ecod_search_v4_negspace_part2/parts/  (leda partial 2)
#     ~/work/prosmos_2026/untested_366_pdb/parts/               (366 fill-in)
#
#   AFDB data sources:
#     ~/work/prosmos_2026/afdb_negspace_sweep_v2/parts/         (all 800, original)
#     ~/work/prosmos_2026/untested_366_afdb/parts/              (366 re-run, sanity)
#
# Where both AFDB sources have the same query, the original v2 value is used
# (the re-run was for parallel-methodology validation, not replacement).

set -euo pipefail

OUT="${1:-/home/rschaeff/work/prosmos_2026/afdb_negspace_final}"
mkdir -p "$OUT"

PDB_DIRS=(
    /home/rschaeff/work/prosmos_2026/ecod_search_v4_negspace/parts
    /home/rschaeff/work/prosmos_2026/ecod_search_v4_negspace_part2/parts
    /home/rschaeff/work/prosmos_2026/untested_366_pdb/parts
)
AFDB_DIRS=(
    /home/rschaeff/work/prosmos_2026/afdb_negspace_sweep_v2/parts
    /home/rschaeff/work/prosmos_2026/untested_366_afdb/parts
)

# Build merged maps. Format: <query>\t<hits>. Take union, prefer first-seen
# value (so original v2 AFDB wins over the re-run, and partial-1 wins over
# overlapping partials).
PDB_MAP=$OUT/pdb_hits.tsv
AFDB_MAP=$OUT/afdb_hits.tsv

build_map() {
    local out="$1"; shift
    : > "$out"
    for d in "$@"; do
        find "$d" -maxdepth 1 -name '*.tsv' -type f 2>/dev/null | xargs cat 2>/dev/null
    done | awk -F'\t' '$3==0 && !($1 in seen) { seen[$1]=$4; print $1"\t"$4 } NF==0 { next }' > "$out"
    wc -l "$out"
}

echo "[$(date)] Building PDB map..."
build_map "$PDB_MAP" "${PDB_DIRS[@]}"
echo "[$(date)] Building AFDB map..."
build_map "$AFDB_MAP" "${AFDB_DIRS[@]}"

# Combine and label.
QUADRANTS=$OUT/quadrants.tsv
awk -F'\t' -v pdb="$PDB_MAP" -v afdb="$AFDB_MAP" '
BEGIN {
    while ((getline line < pdb) > 0) { split(line, a, "\t"); pdb_h[a[1]] = a[2] }
    while ((getline line < afdb) > 0) { split(line, a, "\t"); afdb_h[a[1]] = a[2] }
    queries = ""
    for (q in pdb_h) all[q]=1
    for (q in afdb_h) all[q]=1
}
END {
    printf "query\tpdb_hits\tafdb_hits\tquadrant\n"
    for (q in all) {
        p = (q in pdb_h) ? pdb_h[q] : "NA"
        a = (q in afdb_h) ? afdb_h[q] : "NA"
        if (p == "NA" || a == "NA")    label = "incomplete"
        else if (p == 0 && a == 0)     label = "truly_absent"
        else if (p == 0 && a > 0)      label = "AFDB_only"
        else if (p > 0 && a == 0)      label = "AFDB_loss"
        else                            label = "common"
        printf "%s\t%s\t%s\t%s\n", q, p, a, label
    }
}' /dev/null | sort > "$QUADRANTS"

# Summary report.
SUMMARY=$OUT/summary.txt
{
    echo "AFDB negspace sweep -- FINAL 4-quadrant analysis"
    echo "================================================="
    echo
    echo "PDB sources:"
    for d in "${PDB_DIRS[@]}"; do
        c=$(find "$d" -maxdepth 1 -name '*.tsv' -type f 2>/dev/null | wc -l)
        echo "  $c queries from $d"
    done
    echo "AFDB sources:"
    for d in "${AFDB_DIRS[@]}"; do
        c=$(find "$d" -maxdepth 1 -name '*.tsv' -type f 2>/dev/null | wc -l)
        echo "  $c queries from $d"
    done
    echo
    PDB_COVERED=$(wc -l < "$PDB_MAP")
    AFDB_COVERED=$(wc -l < "$AFDB_MAP")
    echo "Unique queries with PDB hit count:  $PDB_COVERED / 800"
    echo "Unique queries with AFDB hit count: $AFDB_COVERED / 800"
    echo
    echo "Per-quadrant counts:"
    tail -n +2 "$QUADRANTS" | awk -F'\t' '{c[$4]++} END {for (k in c) printf "  %-15s %5d\n", k, c[k]}' | sort
    echo
    echo "AFDB_only candidates (PDB=0, AFDB>0) -- the headline result:"
    awk -F'\t' '$4=="AFDB_only"' "$QUADRANTS" | sort -t$'\t' -k3,3 -rn
    echo
    echo "Top truly_absent candidates by skeleton (PDB=0 AND AFDB=0 across N typings):"
    awk -F'\t' '$4=="truly_absent" { split($1,p,"-"); c[p[1]"-"p[2]]++ } END { for (k in c) printf "  %s\t%2d typings\n", k, c[k] }' "$QUADRANTS" | sort -t$'\t' -k2,2 -rn | head -20
    echo
    echo "AFDB_loss cases (PDB>0, AFDB=0 -- anomalies):"
    awk -F'\t' '$4=="AFDB_loss"' "$QUADRANTS" | sort -t$'\t' -k2,2 -rn | head -20
} > "$SUMMARY"

# Separate files for the key sets.
awk -F'\t' '$4=="AFDB_only"' "$QUADRANTS" > "$OUT/afdb_only.tsv"
awk -F'\t' '$4=="truly_absent"' "$QUADRANTS" > "$OUT/truly_absent.tsv"
awk -F'\t' '$4=="AFDB_loss"' "$QUADRANTS" > "$OUT/afdb_loss.tsv"
awk -F'\t' '$4=="incomplete" {print $1, $2, $3}' "$QUADRANTS" > "$OUT/incomplete.tsv"

echo "[$(date)] Final analysis complete."
echo "Outputs in $OUT/:"
ls -la "$OUT"/*.tsv "$OUT/summary.txt" 2>/dev/null
echo
cat "$SUMMARY"
