#!/usr/bin/env python3
"""Count PDB sequence families (ECOD F-groups) and folds (T-groups) per S5 cell.

The PDB is heavily redundant (~30x: one drug target crystallised hundreds of
times), so a raw domain count is not a scaffold-diversity count. Dedup to the
ECOD F-group -- ECOD's sequence-family level -- to get the number of distinct
sequence families the PDB offers for each local topology, the fair comparator to
AFDB's cluster-rep families. Also keep the T-group (fold) count for reference.
"""
import subprocess, re, glob
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool

R = Path("/home/rschaeff/work/prosmos_2026")
OUT = R / "design_scaffold" / "pdb_cell_families.tsv"
ZSTD = "/sw/apps/Anaconda3-2023.09-0/bin/zstd"
QP = re.compile(r'^s5-(\d{4})-(\d{4})$')

# uid -> (t_group, f_group). Fields: 0=uid 2=t 3=f 4=pdb.
uid_t, uid_f = {}, {}
for line in open(R / "pdb_exp_build" / "uid_class.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4 and f[0].isdigit():
        if f[2]:
            uid_t[int(f[0])] = f[2]
        if f[3]:
            uid_f[int(f[0])] = f[3]
print(f"uid->t {len(uid_t):,}  uid->f {len(uid_f):,}", flush=True)


def scan(fn):
    """(sk,ty) -> set of uids, for one hitpart."""
    d = defaultdict(set)
    p = subprocess.Popen([ZSTD, "-dc", fn], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        t = line.rstrip("\n").split("\t")
        if len(t) < 2:
            continue
        m = QP.match(t[0])
        if m and t[1].isdigit():
            d[(int(m.group(1)), int(m.group(2)))].add(int(t[1]))
    p.wait()
    return d


def main():
    parts = sorted(glob.glob(str(R / "s5_pdb_inv" / "hitparts" / "*.tsv.zst")))
    print(f"hitparts: {len(parts)}", flush=True)
    # records live in exactly one chunk, so per-part uid sets union cleanly
    cell_uids = defaultdict(set)
    with Pool(16) as pool:
        for d in pool.imap_unordered(scan, parts, chunksize=4):
            for c, us in d.items():
                cell_uids[c] |= us
    print(f"cells with PDB hits: {len(cell_uids):,}", flush=True)
    with open(OUT, "w") as fh:
        for (sk, ty), us in sorted(cell_uids.items()):
            ndom = len(us)
            nf = len({uid_f[u] for u in us if u in uid_f})
            nt = len({uid_t[u] for u in us if u in uid_t})
            fh.write(f"{sk}\t{ty}\t{ndom}\t{nf}\t{nt}\n")
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
