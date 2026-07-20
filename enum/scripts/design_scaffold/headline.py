#!/usr/bin/env python3
"""The design headline: over the SHARED S5 topological vocabulary, how much
sequence-diverse scaffold does the 200M add on top of the PDB?

Units:
  AFDB families = non-singleton MMseqs cluster reps lighting a cell (50% id /
                  90% cov). One rep = one sequence family, deduplicated already.
  PDB families  = distinct ECOD F-groups lighting a cell (ECOD's sequence-family
                  level), which dedups the PDB's ~30-40x deposition redundancy.

Granularity caveat, stated up front: MMseqs 50% clusters are FINER than ECOD
F-groups (an F-group can span several 50% clusters), so the AFDB/PDB family ratio
below is an over-estimate of the like-for-like multiplier. The direction and
order of magnitude are robust; firming up the exact factor needs PDB sequences
re-clustered at 50% (the matched-clustering follow-up). Reported as such.

Three products:
  1. the per-cell sequence-diversity multiplier (median, distribution)
  2. the total non-redundant scaffold expansion over shared space
  3. the under-templated cells: real topologies where the PDB offers few scaffold
     families and AFDB offers many -- the design-relevant target list.
"""
import numpy as np
from pathlib import Path

D = Path("/home/rschaeff/work/prosmos_2026/design_scaffold")
G = Path("/home/rschaeff/work/prosmos_2026/s5_grid")
imp = np.load(G / "impossible_mask.npy")
sat = ~imp

af = {}   # (sk,ty) -> families
for line in open(D / "afdb_cell_families.tsv"):
    sk, ty, n = line.split("\t")
    if int(n) > 0:
        af[(int(sk), int(ty))] = int(n)
pf, pt, pd = {}, {}, {}   # F-groups, T-groups, domains
for line in open(D / "pdb_cell_families.tsv"):
    sk, ty, ndom, nF, nT = line.rstrip("\n").split("\t")
    pd[(int(sk), int(ty))] = int(ndom)
    pf[(int(sk), int(ty))] = int(nF)
    pt[(int(sk), int(ty))] = int(nT)

af_cells = set(af); pf_cells = {c for c in pf if pf[c] > 0}
# require PDB F-groups > 0: a handful of shared cells are lit only by
# ECOD-unassigned PDB domains (no F-group), so a family ratio is undefined there.
shared = [c for c in af_cells & pf_cells if sat[c]]
nodef = len([c for c in af_cells & set(pf) if sat[c] and pf[c] == 0])
print(f"cells lit (satisfiable): AFDB {len([c for c in af_cells if sat[c]]):,}  "
      f"PDB {len([c for c in pf_cells if sat[c]]):,}  shared {len(shared):,}  "
      f"(dropped {nodef} shared cells with 0 PDB F-groups)")

# ---- 1. per-cell multiplier over shared vocabulary
mult = np.array([af[c] / pf[c] for c in shared])
afv = np.array([af[c] for c in shared]); pfv = np.array([pf[c] for c in shared])
print("\n=== 1. SEQUENCE-DIVERSITY MULTIPLIER (shared cells, AFDB families / PDB F-groups)")
print(f"  per-cell median families:  PDB {np.median(pfv):.0f}   AFDB {np.median(afv):.0f}")
print(f"  per-cell median multiplier: {np.median(mult):.1f}x   mean {mult.mean():.1f}x")
for q in [10, 25, 50, 75, 90]:
    print(f"    p{q:>2}: {np.percentile(mult, q):5.1f}x")
print(f"  cells where AFDB >= PDB families: {100*(mult>=1).mean():.1f}%")

# ---- 2. total non-redundant scaffold expansion
tot_af = sum(af[c] for c in shared)
tot_pf = sum(pf[c] for c in shared)
tot_pd = sum(pd[c] for c in shared)
print("\n=== 2. TOTAL NON-REDUNDANT SCAFFOLD POOL over shared topological space")
print(f"  PDB domains (redundant):     {tot_pd:,}")
print(f"  PDB families (F-group):      {tot_pf:,}")
print(f"  AFDB families (cluster rep): {tot_af:,}")
print(f"  fold-over-PDB-families:      {tot_af/tot_pf:.1f}x")
print(f"  vs even the redundant PDB domain count: {tot_af/tot_pd:.1f}x")

# ---- 3. under-templated cells: PDB-thin, AFDB-rich, real shared topologies
print("\n=== 3. UNDER-TEMPLATED CELLS (PDB families <= 5, AFDB families >= 50)")
under = [c for c in shared if pf[c] <= 5 and af[c] >= 50]
under.sort(key=lambda c: af[c] / pf[c], reverse=True)
print(f"  count: {len(under)}   (real topologies the PDB barely templates, AFDB richly does)")
print(f"  {'cell':>12} {'PDB_fam':>7} {'PDB_dom':>7} {'AFDB_fam':>8} {'mult':>6}")
for c in under[:15]:
    print(f"  sk{c[0]:03d} ty{c[1]:02d}  {pf[c]:>7} {pd[c]:>7} {af[c]:>8} {af[c]/pf[c]:>5.0f}x")

# also the strict version: how many shared cells is the PDB "thin" on
thin = [c for c in shared if pf[c] <= 3]
print(f"\n  shared cells with <=3 PDB families: {len(thin):,} "
      f"({100*len(thin)/len(shared):.0f}% of shared); "
      f"they carry {sum(af[c] for c in thin):,} AFDB families")

import json
json.dump({
    "shared_cells": len(shared),
    "median_multiplier": float(np.median(mult)),
    "mean_multiplier": float(mult.mean()),
    "tot_afdb_families": tot_af, "tot_pdb_families": tot_pf, "tot_pdb_domains": tot_pd,
    "fold_expansion_vs_families": tot_af / tot_pf,
    "under_templated": [{"sk": c[0], "ty": c[1], "pdb_fam": pf[c], "pdb_dom": pd[c],
                          "afdb_fam": af[c]} for c in under],
}, open(D / "headline.json", "w"), indent=1)
print(f"\nwrote {D/'headline.json'}")
