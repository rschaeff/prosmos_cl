#!/usr/bin/env python3
"""Recount the sequence-diversity multiplier with BOTH databases on an identical
50%/90% MMseqs2 footing -- the matched-clustering correction to the F-group
headline.

PDB families  = distinct 50% clusters among a cell's lit domains (was: F-groups,
                a coarser homology unit that under-counted PDB diversity).
AFDB families = distinct 50% clusters among a cell's grey reps (was: raw rep
                count; reps are structurally distinct so this should barely move,
                which is itself worth confirming).

Coverage caveat carried in the numbers: 90% of PDB lit domains and 99% of AFDB
reps had a sequence to cluster; domains without one drop out of both the cluster
map and the per-cell recount, so the effect is a near-uniform thinning, not a
directional bias.
"""
import pickle, re
from pathlib import Path
from collections import defaultdict
import numpy as np

D = Path("/home/rschaeff/work/prosmos_2026/design_scaffold")
G = Path("/home/rschaeff/work/prosmos_2026/s5_grid")
sat = ~np.load(G / "impossible_mask.npy")


def load_map(tsv):
    """member -> cluster-rep id, from an mmseqs *_cluster.tsv."""
    m = {}
    for line in open(tsv):
        rep, mem = line.rstrip("\n").split("\t")
        m[mem] = rep
    return m

pdb_clu = load_map(D / "pdb_clu50_cluster.tsv")
afdb_clu = load_map(D / "afdb_clu50_cluster.tsv")
print(f"pdb members mapped {len(pdb_clu):,} -> {len(set(pdb_clu.values())):,} clusters")
print(f"afdb members mapped {len(afdb_clu):,} -> {len(set(afdb_clu.values())):,} clusters")

# AFDB rep-merge check: how much did 50% clustering collapse the grey reps?
n_reps = len(afdb_clu)
n_clu = len(set(afdb_clu.values()))
print(f"\nAFDB grey-rep -> 50%-cluster collapse: {n_reps:,} -> {n_clu:,} "
      f"({100*n_clu/n_reps:.1f}% retained) -- reps {'barely merge' if n_clu/n_reps>0.9 else 'merge materially'}")

# ---- per-cell PDB clusters
pdb_cell_uids = pickle.load(open(D / "pdb_cell_uids.pkl", "rb"))
pdb_cell_clu = {}
for c, uids in pdb_cell_uids.items():
    pdb_cell_clu[c] = len({pdb_clu[str(u)] for u in uids if str(u) in pdb_clu})

# ---- per-cell AFDB clusters
afdb_cell_accs = defaultdict(set)
for line in open(D / "afdb_cell_reps.tsv"):
    sk, ty, acc = line.rstrip("\n").split("\t")
    afdb_cell_accs[(int(sk), int(ty))].add(acc)
afdb_cell_clu = {c: len({afdb_clu[a] for a in accs if a in afdb_clu})
                 for c, accs in afdb_cell_accs.items()}

# ---- matched multiplier over shared satisfiable cells
shared = [c for c in set(pdb_cell_clu) & set(afdb_cell_clu)
          if sat[c] and pdb_cell_clu[c] > 0 and afdb_cell_clu[c] > 0]
mult = np.array([afdb_cell_clu[c] / pdb_cell_clu[c] for c in shared])
afv = np.array([afdb_cell_clu[c] for c in shared])
pfv = np.array([pdb_cell_clu[c] for c in shared])
print(f"\n=== MATCHED 50%/90% MULTIPLIER  (shared satisfiable cells: {len(shared):,})")
print(f"  per-cell median 50%-clusters:  PDB {np.median(pfv):.0f}   AFDB {np.median(afv):.0f}")
print(f"  per-cell median multiplier: {np.median(mult):.1f}x   mean {mult.mean():.1f}x")
for q in [10, 25, 50, 75, 90]:
    print(f"    p{q:>2}: {np.percentile(mult, q):5.1f}x")
print(f"  cells where AFDB >= PDB: {100*(mult>=1).mean():.1f}%")

tot_af = int(afv.sum()); tot_pf = int(pfv.sum())
print(f"\n  total over shared space: PDB {tot_pf:,}  AFDB {tot_af:,}  = {tot_af/tot_pf:.1f}x")

# under-templated in matched units
under = [c for c in shared if pdb_cell_clu[c] <= 5 and afdb_cell_clu[c] >= 50]
under.sort(key=lambda c: afdb_cell_clu[c] / pdb_cell_clu[c], reverse=True)
print(f"\n  under-templated (PDB<=5, AFDB>=50 in 50%-clusters): {len(under)} cells")
for c in under[:12]:
    print(f"    sk{c[0]:03d} ty{c[1]:02d}  PDB {pdb_cell_clu[c]:>3}  AFDB {afdb_cell_clu[c]:>4}  "
          f"{afdb_cell_clu[c]/pdb_cell_clu[c]:>4.0f}x")

import json
json.dump({
    "shared_cells": len(shared),
    "afdb_rep_retain_at_50pct": n_clu / n_reps,
    "median_multiplier_matched": float(np.median(mult)),
    "mean_multiplier_matched": float(mult.mean()),
    "tot_afdb_clusters": tot_af, "tot_pdb_clusters": tot_pf,
    "total_ratio_matched": tot_af / tot_pf,
    "under_templated_matched": [{"sk": c[0], "ty": c[1],
                                 "pdb": pdb_cell_clu[c], "afdb": afdb_cell_clu[c]} for c in under],
}, open(D / "matched.json", "w"), indent=1)
print(f"\nwrote {D/'matched.json'}")
