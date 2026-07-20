#!/usr/bin/env python3
"""The decisive control for the n=5 divergence: cross-database Jaccard vs the
WITHIN-database null, both at the same hitter count.

Matched cross-database Jaccard (0.607 at n=5) is low, but that alone cannot tell
"AFDB and PDB sample different topological vocabularies" from "same vocabulary,
both far below saturation at this depth." At n=5, ~50k records light only ~2,500
of 6,336 cells, so even two draws of the SAME database will disagree on the
margins.

So at a common hitter count T, measure three Jaccards:
  cross      : lit(T AFDB hitters) vs lit(T PDB hitters)
  within-PDB : lit(T PDB hitters)  vs lit(another disjoint T PDB hitters)
  within-AFDB: lit(T AFDB hitters) vs lit(another disjoint T AFDB hitters)

If cross ~ the within-database nulls, the divergence is sampling -- one alphabet,
undersampled. If cross sits well BELOW both nulls, it is real vocabulary
divergence and n=5 genuinely locates the boundary.
"""
import subprocess, re, glob, pickle, random
from pathlib import Path
from collections import defaultdict
import numpy as np

R = Path("/home/rschaeff/work/prosmos_2026")
G = R / "s5_grid"
ZSTD = "/sw/apps/Anaconda3-2023.09-0/bin/zstd"
QP = re.compile(r'^s([345])-(\d{4})-(\d{4})$')
DRAWS = 8

afdb_sub = pickle.load(open(G / "afdb_sub_records.pkl", "rb"))
pdb_sub = pickle.load(open(G / "pdb_sub_records.pkl", "rb"))


def scan(pattern, sub):
    rec_cells = defaultdict(set)
    for fn in sorted(glob.glob(pattern)):
        p = subprocess.Popen([ZSTD, "-dc", fn], stdout=subprocess.PIPE, text=True)
        for line in p.stdout:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 2:
                m = QP.match(t[0])
                if m and t[1] in sub:
                    rec_cells[t[1]].add((int(m.group(2)), int(m.group(3))))
        p.wait()
    return rec_cells


def afdb_n5(sub):
    rec_cells = defaultdict(set)
    for line in open(R / "design_scaffold" / "afdb_cell_reps.tsv"):
        sk, ty, acc = line.rstrip("\n").split("\t")
        if acc in sub:
            rec_cells[acc].add((int(sk), int(ty)))
    return rec_cells


def lit(sample, rc):
    s = set()
    for r in sample:
        s |= rc[r]
    return s


def jac(a, b):
    return len(a & b) / len(a | b) if (a or b) else float("nan")


print(f"{'n':>2} {'T':>7} {'cross':>14} {'within-PDB':>14} {'within-AFDB':>14}   verdict", flush=True)
rows = []
for n in (3, 4, 5):
    arc = afdb_n5(afdb_sub) if n == 5 else scan(str(R / f"n{n}_afdb" / "hitparts" / "*.tsv.zst"), afdb_sub)
    prc = scan(str(R / ("s5_pdb_inv" if n == 5 else f"n{n}_pdb") / "hitparts" / "*.tsv.zst"), pdb_sub)
    A, P = list(arc), list(prc)
    # T small enough to draw two DISJOINT samples from each database
    T = min(len(A), len(P)) // 2
    rng = random.Random(0)
    cross, wp, wa = [], [], []
    for d in range(DRAWS):
        sa = rng.sample(A, 2 * T); sp = rng.sample(P, 2 * T)
        a1, a2 = lit(sa[:T], arc), lit(sa[T:], arc)
        p1, p2 = lit(sp[:T], prc), lit(sp[T:], prc)
        cross.append(jac(a1, p1))
        wp.append(jac(p1, p2))
        wa.append(jac(a1, a2))
    mc, mp, ma = np.mean(cross), np.mean(wp), np.mean(wa)
    null = min(mp, ma)
    verdict = ("SAMPLING (cross within null)" if mc >= null - 0.03
               else "REAL DIVERGENCE (cross below null)")
    rows.append({"n": n, "T": T, "cross": float(mc), "within_pdb": float(mp),
                 "within_afdb": float(ma), "verdict": verdict})
    print(f"{n:>2} {T:>7,} {mc:>8.3f}±{np.std(cross):.3f} {mp:>8.3f}±{np.std(wp):.3f} "
          f"{ma:>8.3f}±{np.std(wa):.3f}   {verdict}", flush=True)

import json
json.dump(rows, open(G / "jaccard_null.json", "w"), indent=1)
print(f"\nwrote {G/'jaccard_null.json'}")
