#!/bin/bash
# generateMatrix regression: validate our rebuilt pipeline (palsse_cl .ssd +
# generateMatrix/build) against the ORIGINAL ProSMoS 2010 metamatricesDB — the
# frozen gold-standard build by the ProSMoS authors. For each structure we have
# a PALSSE .ssd for that also exists in the 2010 DB, run our generateMatrix and
# compare to the 2010 DB record on the BUILD-CORRECTNESS invariants:
#   (1) the contact-matrix line  (what searchmatrix consumes) -- must be byte-equal
#   (2) the SSE type+range+length sequence
# Our build appends a chain char to the type ("EA" vs legacy "E"); that and the
# raw endpoint COORDINATES (input PDB version/frame; contacts are invariant) are
# reported as INFO, not FAIL.
set -u
OURS=/home/rschaeff/dev/prosmos_cl/generateMatrix/build/generateMatrix
DB2010=/home/rschaeff/src/Prosmos/ProSMoS/metamatrixdb/metamatricesDB
SSD_DIR=/home/rschaeff/work/prosmos_2026/palsse_oracle/out
PY=/sw/apps/Anaconda3-2023.09-0/bin/python
W=$(mktemp -d); cd "$W"
UNSET=$(env | awk -F= '/^(SLURM|PMI|PMIX|OMPI)_/{printf " -u %s", $1}')

for ssd in "$SSD_DIR"/*.ssd; do
    uid=$(basename "$ssd" .ssd)
    grep -q "^${uid}[.]ssd" "$DB2010" || continue
    cp "$ssd" .
    env $UNSET "$OURS" -os "${uid}.ssd" ./ "${uid}.out" >/dev/null 2>&1
    awk -v u="^${uid}[.]ssd" '$0~u{p=1} p{print; if(/^\*/){exit}}' "$DB2010" > db.rec
    "$PY" - "${uid}.out" db.rec "$uid" <<'PY'
import re, sys
def parse(path):
    txt=open(path).read()
    mat=next((l for l in txt.splitlines() if l.startswith('*')), '')
    # SSE tokens: <type><optional chain><start> -- <end> <len> ; drop chain+coords
    sses=re.findall(r'([HEL])[A-Za-z0-9]?\s*(\S+?)\s*--\s*(\S+?)\s+(\d+)\s+-?\d+\.', txt)
    ranges=[(t,s,e,l) for t,s,e,l in sses]
    # raw coord signature (first coord of each SSE) for the info flag
    coords=re.findall(r'\s(-?\d+\.\d+)', txt)
    return mat, ranges, coords
mo,ro,co=parse(sys.argv[1]); md,rd,cd=parse(sys.argv[2]); uid=sys.argv[3]
mat_ok = mo==md
rng_ok = ro==rd
crd_ok = co==cd
if mat_ok and rng_ok:
    note="coords identical" if crd_ok else "coords differ (input PDB frame; contacts invariant)"
    print(f"PASS {uid}  (contact-matrix + {len(ro)} SSE ranges == 2010 DB; {note})")
else:
    print(f"FAIL {uid}  (matrix={'ok' if mat_ok else 'DIFF'} ranges={'ok' if rng_ok else 'DIFF'})")
PY
done | tee /tmp/gmreg.$$
pass=$(grep -c '^PASS' /tmp/gmreg.$$); fail=$(grep -c '^FAIL' /tmp/gmreg.$$); rm -f /tmp/gmreg.$$
echo "=== $pass pass, $fail fail ==="
rm -rf "$W"
