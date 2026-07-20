#!/usr/bin/env python3
"""Capture AFDB rep NAMES per cell (not just counts) so we can recount cells in
50%-cluster units after the matched clustering. Each file in hits/<cell>/ is
pdb<accession>.txt; the accession is the grey non-singleton cluster rep.
"""
import os, re
from pathlib import Path
from multiprocessing import Pool

HITS = Path("/home/rschaeff/work/prosmos_2026/s5_inv/hits")
OUT = Path("/home/rschaeff/work/prosmos_2026/design_scaffold/afdb_cell_reps.tsv")
QP = re.compile(r'^s5-(\d{4})-(\d{4})$')


def reps(cell):
    p = HITS / cell
    m = QP.match(cell)
    sk, ty = int(m.group(1)), int(m.group(2))
    out = []
    try:
        for e in os.scandir(p):
            n = e.name
            if n.endswith(".txt"):
                acc = n[3:-4] if n.startswith("pdb") else n[:-4]
                out.append(f"{sk}\t{ty}\t{acc}")
    except OSError:
        pass
    return out


def main():
    cells = [d for d in os.listdir(HITS) if QP.match(d)]
    with Pool(16) as pool, open(OUT, "w") as fh:
        for chunk in pool.imap_unordered(reps, cells, chunksize=16):
            if chunk:
                fh.write("\n".join(chunk) + "\n")
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
