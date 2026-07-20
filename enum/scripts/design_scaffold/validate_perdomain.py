#!/usr/bin/env python3
"""Validate the per-protein pLDDT gate against exact per-domain pLDDT.

Two questions:
  1. Agreement -- how well does the per-protein proxy track the true per-domain
     value? (correlation, and disagreement rate at the 70 threshold)
  2. Does the surviving-cell verdict hold under exact per-domain pLDDT? Recount
     the 31 surviving cells' confident-scaffold counts and re-test PDB<=5 &
     AFDB>=50.
"""
import json
from pathlib import Path
from collections import defaultdict
import numpy as np

D = Path("/home/rschaeff/work/prosmos_2026/design_scaffold")

# exact per-domain pLDDT
pd_pl = {}
for line in open(D / "surv_plddt.tsv"):
    did, pl = line.rstrip("\n").split("\t")
    pd_pl[did] = float(pl)
print(f"per-domain pLDDT extracted: {len(pd_pl):,}")

# per-protein proxy used earlier: rebuild scaffold->proxy from plddt_filter.json cells?
# simpler: reload the proxy by re-deriving from the filter's accession map is not
# stored, so compare only where we have both. Load the proxy from a cached dump.
proxy = {}
pj = json.loads((D / "proxy_scaffold_plddt.json").read_text()) if (D / "proxy_scaffold_plddt.json").exists() else {}
for k, v in pj.items():
    proxy[k] = v

# ---- agreement (only scaffolds with both)
both = [(pd_pl[d], proxy[d]) for d in pd_pl if d in proxy and proxy[d] is not None]
if both:
    a = np.array([x[0] for x in both]); b = np.array([x[1] for x in both])
    r = np.corrcoef(a, b)[0, 1]
    # threshold agreement at 70
    ta = (a >= 70); tb = (b >= 70)
    print(f"\nAGREEMENT (n={len(both):,}): pearson r={r:.3f}")
    print(f"  mean |per-domain - per-protein| = {np.abs(a-b).mean():.1f} pLDDT")
    print(f"  same side of 70 gate: {100*(ta==tb).mean():.1f}%")
    print(f"  per-domain median {np.median(a):.1f}  per-protein median {np.median(b):.1f}")

# ---- recount surviving cells with EXACT per-domain pLDDT
surv = {tuple(c) for c in json.loads((D / "surv_cells.json").read_text())}
cell_sc = defaultdict(set)
for line in open(D / "afdb_cell_reps.tsv"):
    sk, ty, acc = line.rstrip("\n").split("\t")
    k = (int(sk), int(ty))
    if k in surv:
        cell_sc[k].add(acc)

pf = {(c["sk"], c["ty"]): c for c in json.load(open(D / "plddt_filter.json"))["cells"]}
print(f"\n{'cell':>12} {'PDB':>4} {'AFDB>=70(prot)':>13} {'AFDB>=70(dom)':>13} {'>=80(dom)':>9} {'cov':>5}")
still70 = still80 = 0
rows = []
for k in sorted(surv, key=lambda k: -pf[k]["afdb_ge70"]):
    scs = cell_sc[k]
    covered = [s for s in scs if s in pd_pl]
    n70 = sum(1 for s in covered if pd_pl[s] >= 70)
    n80 = sum(1 for s in covered if pd_pl[s] >= 80)
    pdb = pf[k]["pdb"]
    cov = len(covered) / max(len(scs), 1)
    if pdb <= 5 and n70 >= 50:
        still70 += 1
    if pdb <= 5 and n80 >= 50:
        still80 += 1
    rows.append((k, pdb, pf[k]["afdb_ge70"], n70, n80, cov))
for k, pdb, p70, n70, n80, cov in rows[:31]:
    print(f"  sk{k[0]:03d} ty{k[1]:02d} {pdb:>4} {p70:>13} {n70:>13} {n80:>9} {cov:>4.0%}")

print(f"\nsurviving cells still design-grade under EXACT per-domain pLDDT>=70 (AFDB>=50): {still70}/31")
print(f"  ... at >=80: {still80}/31")

json.dump({"n_perdomain": len(pd_pl),
           "still70": still70, "still80": still80,
           "cells": [{"sk": k[0], "ty": k[1], "pdb": pdb, "afdb70_dom": n70,
                      "afdb80_dom": n80, "coverage": cov}
                     for k, pdb, p70, n70, n80, cov in rows]},
          open(D / "perdomain_validated.json", "w"), indent=1)
print(f"wrote {D/'perdomain_validated.json'}")
