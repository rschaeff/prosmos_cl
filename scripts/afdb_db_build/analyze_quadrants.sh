#!/bin/bash
# Cross-reference the AFDB negspace sweep results against the leda PDB-side
# negspace sweep results and produce the 4-quadrant classification.
#
# Usage:
#   analyze_quadrants.sh <AFDB_SWEEP_OUT> <OUT_DIR>
#
# Inputs:
#   AFDB_SWEEP_OUT/summary.tsv          — the 800-query AFDB sweep
#   leda PDB negspace parts dirs        — union of the two partial runs
#
# Outputs (written to OUT_DIR):
#   quadrants.tsv      — per-query: PDB hits, AFDB hits, quadrant label
#   summary.txt        — counts per quadrant + top AFDB-only candidates
#   afdb_only.txt      — the headline result: queries with PDB=0 AND AFDB>0
#   afdb_only.queries  — same as above, just the query names (one per line)

set -euo pipefail

AFDB_OUT="${1:?missing AFDB sweep OUT}"
ANALYSIS="${2:?missing analysis OUT dir}"

PDB_PARTS_A=/home/rschaeff/work/prosmos_2026/ecod_search_v4_negspace/parts
PDB_PARTS_B=/home/rschaeff/work/prosmos_2026/ecod_search_v4_negspace_part2/parts

mkdir -p "$ANALYSIS"

AFDB_SUM="$AFDB_OUT/summary.tsv"
[ -s "$AFDB_SUM" ] || { echo "AFDB summary missing: $AFDB_SUM" >&2; exit 2; }

# Build PDB hit map: query -> hits, from union of two partial sweep dirs.
# Each parts/<line>.tsv has "<query>\t<elapsed>\t<rc>\t<hits>".
PDB_MAP=$(mktemp)
{
    cat "$PDB_PARTS_A"/*.tsv 2>/dev/null
    cat "$PDB_PARTS_B"/*.tsv 2>/dev/null
} | awk -F'\t' '$3==0 {print $1"\t"$4}' | sort -u > "$PDB_MAP"
PDB_COVERED=$(wc -l < "$PDB_MAP")
echo "[$(date)] PDB-side leda coverage: $PDB_COVERED / 800 queries"

awk -F'\t' -v pdb_map="$PDB_MAP" '
BEGIN {
    while ((getline line < pdb_map) > 0) {
        split(line, a, "\t")
        pdb[a[1]] = a[2]
    }
    close(pdb_map)
}
NR > 1 {
    q = $1; afdb = $4
    if (q in pdb) {
        p = pdb[q]
        if (p == 0 && afdb == 0)      label = "truly_absent"
        else if (p == 0 && afdb > 0)  label = "AFDB_only"
        else if (p > 0 && afdb == 0)  label = "AFDB_loss"
        else                           label = "common"
        printf "%s\t%s\t%s\t%s\n", q, p, afdb, label
    } else {
        printf "%s\t%s\t%s\t%s\n", q, "NA", afdb, "PDB_untested"
    }
}' "$AFDB_SUM" | sort > "$ANALYSIS/quadrants.tsv"

# Per-quadrant counts.
{
    echo "AFDB negspace sweep -- 4-quadrant analysis"
    echo "=========================================="
    echo "AFDB summary:   $AFDB_SUM"
    echo "PDB leda sweep: union of v4_negspace + v4_negspace_part2 ($PDB_COVERED / 800 covered)"
    echo
    echo "Per-quadrant counts:"
    awk -F'\t' '{c[$4]++} END {for (k in c) printf "  %-15s %5d\n", k, c[k]}' "$ANALYSIS/quadrants.tsv" | sort
    echo
    echo "Top AFDB-only candidates (PDB=0, AFDB>0) by AFDB hit count:"
    awk -F'\t' '$4=="AFDB_only"' "$ANALYSIS/quadrants.tsv" | sort -t$'\t' -k3,3 -rn | head -20
} > "$ANALYSIS/summary.txt"

awk -F'\t' '$4=="AFDB_only"' "$ANALYSIS/quadrants.tsv" | sort -t$'\t' -k3,3 -rn > "$ANALYSIS/afdb_only.txt"
awk -F'\t' '$4=="AFDB_only" {print $1}' "$ANALYSIS/quadrants.tsv" | sort > "$ANALYSIS/afdb_only.queries"

rm -f "$PDB_MAP"

echo "[$(date)] Analysis complete. See $ANALYSIS/{summary.txt,quadrants.tsv,afdb_only.{txt,queries}}"
cat "$ANALYSIS/summary.txt"
