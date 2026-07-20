#!/usr/bin/env python3
"""Count AFDB sequence families per S5 cell.

The AFDB sweep ran in tree mode: hits/<cell>/pdb<record>.txt, one file per
(cell, record). searchmatrix writes exactly one file per hit RECORD (multiple
motif blocks live inside it), so the file count in a cell dir is the number of
distinct AFDB records lighting that cell.

Those records are non-singleton MMseqs cluster representatives (one per cluster,
50% id / 90% cov), so **file count = number of distinct sequence families that
realise this local topology** -- exactly the design-relevant scaffold-diversity
quantity, already deduplicated by construction.

os.scandir over NFS is the cost; parallelise across cell dirs.
"""
import os, re, sys
from pathlib import Path
from multiprocessing import Pool

HITS = Path("/home/rschaeff/work/prosmos_2026/s5_inv/hits")
OUT = Path("/home/rschaeff/work/prosmos_2026/design_scaffold/afdb_cell_families.tsv")
OUT.parent.mkdir(parents=True, exist_ok=True)
QP = re.compile(r'^s5-(\d{4})-(\d{4})$')


def count(cell):
    p = HITS / cell
    try:
        n = sum(1 for e in os.scandir(p) if e.name.endswith(".txt"))
    except OSError:
        n = -1
    return cell, n


def main():
    cells = [d for d in os.listdir(HITS) if QP.match(d)]
    print(f"cell dirs: {len(cells):,}", flush=True)
    with Pool(16) as pool:
        res = pool.map(count, cells, chunksize=32)
    bad = [c for c, n in res if n < 0]
    with open(OUT, "w") as fh:
        for cell, n in sorted(res):
            m = QP.match(cell)
            fh.write(f"{int(m.group(1))}\t{int(m.group(2))}\t{n}\n")
    tot = sum(n for _, n in res if n >= 0)
    print(f"wrote {OUT}  cells={len(res):,}  total AFDB family-hits={tot:,}  scandir-failed={len(bad)}",
          flush=True)


if __name__ == "__main__":
    main()
