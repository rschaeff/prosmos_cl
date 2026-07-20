#!/usr/bin/env python3
"""Size-matched Jaccard: is the n=5 AFDB<->PDB divergence real, or just depth?

The raw curve shows AFDB and PDB lighting near-identical cell sets at n=3,4
(Jaccard ~0.97) then diverging at n=5 (0.674). But the sub-record sets are not
depth-matched -- 491,963 AFDB vs 49,640 PDB, a ~10x gap -- and at n=5 the grid is
unsaturated, so the deeper AFDB lights rare cells the PDB simply has not reached.
That is the same AFDB-only asymmetry we already showed is mostly sampling depth.

To isolate true vocabulary divergence, thin the AFDB SEARCHED population to PDB's
depth (keep each AFDB hitting record with probability q = |pdb_sub|/|afdb_sub|),
recompute the AFDB-lit cell set, and re-measure Jaccard vs PDB. Same method as
afdb_rarefy.py. Multiple draws.

  raw Jaccard drops but matched Jaccard stays flat  => the n=5 divergence is DEPTH.
  matched Jaccard still falls                        => real alphabet divergence,
                                                        and n=5 locates the boundary.
"""
import subprocess, re, glob, pickle, random
from pathlib import Path
from collections import defaultdict
import numpy as np

R = Path("/home/rschaeff/work/prosmos_2026")
G = R / "s5_grid"
ZSTD = "/sw/apps/Anaconda3-2023.09-0/bin/zstd"
QP = re.compile(r'^s([345])-(\d{4})-(\d{4})$')
TOTAL = {3: 40, 4: 672, 5: 6336}
DRAWS = 5

afdb_sub = pickle.load(open(G / "afdb_sub_records.pkl", "rb"))
pdb_sub = pickle.load(open(G / "pdb_sub_records.pkl", "rb"))
q = len(pdb_sub) / len(afdb_sub)
print(f"sub-records: AFDB {len(afdb_sub):,}  PDB {len(pdb_sub):,}  ->  keep fraction q = {q:.4f}", flush=True)


def scan(pattern, sub):
    out = defaultdict(set)
    for fn in sorted(glob.glob(pattern)):
        p = subprocess.Popen([ZSTD, "-dc", fn], stdout=subprocess.PIPE, text=True)
        for line in p.stdout:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 2:
                m = QP.match(t[0])
                if m and t[1] in sub:
                    out[(int(m.group(2)), int(m.group(3)))].add(t[1])
        p.wait()
    return out


def afdb_n5(sub):
    out = defaultdict(set)
    for line in open(R / "design_scaffold" / "afdb_cell_reps.tsv"):
        sk, ty, acc = line.rstrip("\n").split("\t")
        if acc in sub:
            out[(int(sk), int(ty))].add(acc)
    return out


print(f"\n{'n':>2} {'raw J':>7} {'matched J':>18} {'shared':>7} {'AFDBonly':>9} {'->matched':>10}", flush=True)
rows = []
for n in (3, 4, 5):
    ac = afdb_n5(afdb_sub) if n == 5 else scan(str(R / f"n{n}_afdb" / "hitparts" / "*.tsv.zst"), afdb_sub)
    pc = scan(str(R / ("s5_pdb_inv" if n == 5 else f"n{n}_pdb") / "hitparts" / "*.tsv.zst"), pdb_sub)
    p_lit = {c for c, s in pc.items() if s}
    a_lit_full = {c for c, s in ac.items() if s}
    raw_j = len(a_lit_full & p_lit) / len(a_lit_full | p_lit)

    # record -> cells it lights (AFDB), for fast re-lighting under a record sample
    rec_cells = defaultdict(set)
    for c, recs in ac.items():
        for r in recs:
            rec_cells[r].add(c)
    hitters = list(rec_cells)
    rng = random.Random(0)
    js, aos, shs = [], [], []
    for d in range(DRAWS):
        k = max(1, int(len(hitters) * q))
        keep = rng.sample(hitters, k)
        a_lit = set()
        for r in keep:
            a_lit |= rec_cells[r]
        js.append(len(a_lit & p_lit) / len(a_lit | p_lit))
        aos.append(len(a_lit - p_lit))
        shs.append(len(a_lit & p_lit))
    mj, sj = float(np.mean(js)), float(np.std(js))
    rows.append({"n": n, "raw_jaccard": raw_j, "matched_jaccard": mj, "matched_sd": sj,
                 "raw_afdb_only": len(a_lit_full - p_lit),
                 "matched_afdb_only": float(np.mean(aos)),
                 "matched_shared": float(np.mean(shs)), "pdb_lit": len(p_lit)})
    print(f"{n:>2} {raw_j:>7.3f} {mj:>13.3f} ±{sj:.3f} {np.mean(shs):>7.0f} "
          f"{len(a_lit_full-p_lit):>9} {np.mean(aos):>10.0f}", flush=True)

print("\nread: matched J ~ raw J at n=3,4 (saturated, depth invisible).")
print("      at n=5, if matched J jumps back toward ~0.95 the divergence was depth;")
print("      if it stays low, it is real vocabulary divergence.")
import json
json.dump(rows, open(G / "size_matched_jaccard.json", "w"), indent=1)
print(f"\nwrote {G/'size_matched_jaccard.json'}")
