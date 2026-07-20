#!/usr/bin/env python3
"""Do the 'AFDB-only' cells concentrate in low-confidence AlphaFold models?

pLDDT comes from the B-factor column of the PRE-CHOPPED DPAM domain structures we
extracted -- i.e. averaged over exactly the residues of the domain that matched.
(Not ted_domain.plddt: TED draws different boundaries than DPAM, so that mean is
over partly the wrong range.)

This is a WITHIN-AFDB contrast, so it is free of the cross-database depth
confound that dominates every AFDB-vs-PDB comparison: both strata are drawn from
the same sweep at the same depth. If AFDB-only cells are model artifacts their
domains should be systematically less confident than domains lighting shared cells.
"""
import json, glob, os
import numpy as np
from collections import defaultdict

G = "/home/rschaeff/work/prosmos_2026/s5_grid"
DOM = "/home/rschaeff/work/prosmos_2026/afdb_domain_struct"
CELLS = "/home/rschaeff/dev/prosmos_inspect/data/afdb_assigned/cell"

A = np.load(f"{G}/grid_afdb_nT_rebuilt.npy")
P = np.load(f"{G}/grid_pdb_nT.npy").astype(int)
IMP = np.load(f"{G}/impossible_mask.npy")
SAT = ~IMP
afdb_only = (A > 0) & (P == 0) & SAT
shared = (A > 0) & (P > 0) & SAT
print(f"cells: AFDB-only {afdb_only.sum():,}  shared {shared.sum():,}")

# per-domain mean CA pLDDT
plddt = {}
for f in glob.glob(f"{DOM}/*.pdb"):
    did = os.path.basename(f)[:-4]
    v = [float(l[60:66]) for l in open(f, errors="replace")
         if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    if v:
        plddt[did] = sum(v) / len(v)
print(f"domains with pLDDT: {len(plddt):,}")

# which cells does each exemplar domain light?
dom_cells = defaultdict(list)
for cf in sorted(glob.glob(f"{CELLS}/*.json")):
    d = json.load(open(cf))
    sk, ty = d["sk"], d["ty"]
    for e in d["exemplars"]:
        dom_cells[e["did"]].append((sk, ty))

only_v, shared_v, both_v = [], [], []
for did, cs in dom_cells.items():
    p = plddt.get(did)
    if p is None:
        continue
    o = any(afdb_only[sk, ty] for sk, ty in cs)
    s = any(shared[sk, ty] for sk, ty in cs)
    if o and not s:
        only_v.append(p)
    elif s and not o:
        shared_v.append(p)
    elif o and s:
        both_v.append(p)


def desc(name, v):
    v = np.array(v)
    if not v.size:
        print(f"  {name}: none"); return
    print(f"  {name:32s} n={v.size:>6,}  mean {v.mean():5.1f}  median {np.median(v):5.1f}  "
          f"<70 {100*(v<70).mean():4.1f}%  <50 {100*(v<50).mean():4.1f}%")


print("\nmean pLDDT of exemplar domains, by the cells they light:")
desc("lights AFDB-only cells ONLY", only_v)
desc("lights shared cells ONLY", shared_v)
desc("lights both", both_v)

from scipy.stats import mannwhitneyu
if only_v and shared_v:
    u, pv = mannwhitneyu(only_v, shared_v, alternative="two-sided")
    print(f"\n  Mann-Whitney AFDB-only vs shared: p={pv:.3g}")
    print(f"  difference in means: {np.mean(only_v)-np.mean(shared_v):+.1f} pLDDT")

# the cells that survived the depth control / both significance tests
pers = None
try:
    import pickle
    pers = pickle.load(open(f"{G}/afdb_persist.pkl", "rb"))
except Exception:
    pass
if pers:
    rob = set(pers["robust"]); frag = set(pers["fragile"])
    rv, fv = [], []
    for did, cs in dom_cells.items():
        p = plddt.get(did)
        if p is None:
            continue
        if any((sk, ty) in rob for sk, ty in cs):
            rv.append(p)
        if any((sk, ty) in frag for sk, ty in cs):
            fv.append(p)
    print("\nwithin the AFDB-only set, by depth-control outcome:")
    desc("ROBUST cells (survive thinning)", rv)
    desc("FRAGILE cells (depth artifacts)", fv)

four = [(16, 12), (34, 0), (158, 27), (171, 7)]
fv4 = [plddt[d] for d, cs in dom_cells.items() if plddt.get(d) and any(c in four for c in cs)]
if fv4:
    desc("the 4 doubly-significant cells", fv4)
json.dump({"only": only_v, "shared": shared_v, "both": both_v},
          open(f"{G}/plddt_strat.json", "w"))
print("\nwrote plddt_strat.json")
