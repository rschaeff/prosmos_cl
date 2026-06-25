#!/bin/bash
# Clean 4-quadrant analysis combining the post-fix PDB-experimental sweep
# (pdb_exp_negspace_sweep_v2) and the post-fix AFDB-non-singleton sweep
# (afdb_negspace_sweep_v3). Both sources used the patched searchmatrix
# binary, so hit counts are real (no rc=139/rc=134 crash contamination).
#
# Usage: bash quadrant_v2.sh <OUT_DIR>

set -euo pipefail

OUT="${1:-/home/rschaeff/work/prosmos_2026/afdb_negspace_final_v2}"
PDB_SUM=/home/rschaeff/work/prosmos_2026/pdb_exp_negspace_sweep_v2/summary.tsv
AFDB_SUM=/home/rschaeff/work/prosmos_2026/afdb_negspace_sweep_v3/summary.tsv
mkdir -p "$OUT"

[ -s "$PDB_SUM" ]  || { echo "missing $PDB_SUM"  >&2; exit 2; }
[ -s "$AFDB_SUM" ] || { echo "missing $AFDB_SUM" >&2; exit 2; }

# Per-query maps (only successful, rc=0 rows).
awk -F'\t' 'NR>1 && $3==0 {print $1"\t"$4}' "$PDB_SUM"  | sort > "$OUT/pdb_hits.tsv"
awk -F'\t' 'NR>1 && $3==0 {print $1"\t"$4}' "$AFDB_SUM" | sort > "$OUT/afdb_hits.tsv"

python3 <<EOF > "$OUT/quadrants.tsv"
pdb = {}
afdb = {}
with open("$OUT/pdb_hits.tsv") as f:
    for line in f:
        q, h = line.rstrip().split("\t")
        pdb[q] = int(h)
with open("$OUT/afdb_hits.tsv") as f:
    for line in f:
        q, h = line.rstrip().split("\t")
        afdb[q] = int(h)

print("query\tpdb_exp_hits\tafdb_hits\tquadrant")
for q in sorted(set(pdb) | set(afdb)):
    p = pdb.get(q, "NA")
    a = afdb.get(q, "NA")
    if p == "NA" or a == "NA":
        label = "incomplete"
    elif p == 0 and a == 0:
        label = "truly_absent"
    elif p == 0 and a > 0:
        label = "AFDB_only"
    elif p > 0 and a == 0:
        label = "AFDB_loss"
    else:
        label = "common"
    print(f"{q}\t{p}\t{a}\t{label}")
EOF

# Per-quadrant counts.
{
    echo "Experimental-PDB vs AFDB-non-singleton 4-quadrant analysis (post-fix)"
    echo "===================================================================="
    echo
    echo "PDB source:  pdb_exp_negspace_sweep_v2 (496,359 experimental-PDB entries)"
    echo "AFDB source: afdb_negspace_sweep_v3    (4,921,931 non-singleton AFDB entries)"
    echo "Binary:      searchmatrix d3e1e83 (pid+matrix buffer overflows fixed)"
    echo
    PDB_OK=$(wc -l < "$OUT/pdb_hits.tsv")
    AFDB_OK=$(wc -l < "$OUT/afdb_hits.tsv")
    echo "Queries with PDB hit count (rc=0):  $PDB_OK / 800"
    echo "Queries with AFDB hit count (rc=0): $AFDB_OK / 800"
    echo
    echo "Per-quadrant counts:"
    tail -n +2 "$OUT/quadrants.tsv" | awk -F'\t' '{c[$4]++} END {for (k in c) printf "  %-15s %5d\n", k, c[k]}' | sort
    echo
    echo "AFDB-only candidates (PDB=0, AFDB>0) -- the headline:"
    awk -F'\t' '$4=="AFDB_only"' "$OUT/quadrants.tsv" | sort -t$'\t' -k3,3 -rn | head -20
    echo
    echo "AFDB-loss anomalies (PDB>0, AFDB=0):"
    awk -F'\t' '$4=="AFDB_loss"' "$OUT/quadrants.tsv" | sort -t$'\t' -k2,2 -rn | head -20
    echo
    echo "Common (lit up in both):"
    awk -F'\t' '$4=="common"' "$OUT/quadrants.tsv" | sort -t$'\t' -k2,2 -rn | head -10
} > "$OUT/summary.txt"

# Split into per-quadrant files
for q in truly_absent AFDB_only AFDB_loss common incomplete; do
    awk -F'\t' -v q="$q" '$4==q' "$OUT/quadrants.tsv" > "$OUT/${q}.tsv"
done

echo "Done. Output: $OUT/"
cat "$OUT/summary.txt"
