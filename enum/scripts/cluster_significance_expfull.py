#!/usr/bin/env python3
"""Cluster-unit AFDB-only significance against the corrected 1.4M experimental sweep.

WHY A CLUSTER UNIT AT ALL. The per-structure / per-fold depth correction still
over-calls, because AFDB is massively sequence-redundant: one topology realized by
a 10,000-member sequence family counts as 10,000 chances to light a cell. The
matched-clustering result (Ruczinski 20->10, design-scaffold 8x->4.6x) showed the
honest unit collapses that redundancy. Here we redo the Fisher + Poisson-binomial
depth correction with a SEQUENCE CLUSTER as the count unit.

UNITS (asymmetric thresholds -- flagged, not hidden):
  AFDB : 50%-id / 90%-cov mmseqs clusters (afdb_clu50), the published unit.
         cell_clu_k carries per-cluster member counts, so the depth thinning
         deflates a big family to ~1 surviving cluster instead of thousands of
         structures. This side is UNCHANGED from cluster_significance.py.
  PDB  : 50%-id / 90%-cov mmseqs clusters of the new 1.4M lit domains
         (pdb_clu50_expfull_cluster.tsv, built by build_pdb_clu50_expfull.py),
         MATCHED to afdb_clu50's threshold. An earlier F-group attempt gave only
         10,358 families vs 1.18M AFDB clusters -- a 100x scale asymmetry that
         INFLATED the depth-sig count (74->185, wrong direction) because AFDB was
         counted in fine units and PDB in coarse ones. Matched thresholds fix it.

The depth correction (the defensible test) uses only AFDB cluster member-counts
and the PDB observed count, so the threshold asymmetry bears on the Fisher margins
more than on the depth q. Reports both.

Compares against significance_expfull.json (the fold-unit recompute) so the
over-call shrinkage is explicit -- the whole point.
"""
import glob
import json
import re
import subprocess
from collections import Counter, defaultdict
from math import erf, sqrt
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

R = Path("/home/rschaeff/work/prosmos_2026")
G = R / "s5_grid"
D = R / "design_scaffold"
ZSTD = "/sw/apps/Anaconda3-2023.09-0/bin/zstd"
QP = re.compile(r"^s5-(\d{4})-(\d{4})$")
NSK, NTY = 198, 32
SAT = ~np.load(G / "impossible_mask.npy")


def loadmap(path):
    m = {}
    for line in open(path):
        rep, mem = line.rstrip("\n").split("\t")
        m[mem] = rep
    return m


# ---- AFDB 50%-cluster machinery (UNCHANGED) ----------------------------------
afdb_clu = loadmap(D / "afdb_clu50_cluster.tsv")
NA = len({*afdb_clu.values()})
A_clu = np.zeros((NSK, NTY), int)
cell_clu_k = defaultdict(Counter)
tmp = defaultdict(set)
for line in open(D / "afdb_cell_reps.tsv"):
    sk, ty, acc = line.rstrip("\n").split("\t")
    c = afdb_clu.get(acc)
    if c:
        cell_clu_k[(int(sk), int(ty))][c] += 1
        tmp[(int(sk), int(ty))].add(c)
for c, s in tmp.items():
    A_clu[c] = len(s)
print(f"AFDB 50%-clusters {NA:,}; AFDB-lit cells {int((A_clu > 0).sum()):,}", flush=True)

# ---- PDB 50%-cluster (MATCHED to AFDB) from the NEW 1.4M sweep -----------------
pdb_clu = loadmap(D / "pdb_clu50_expfull_cluster.tsv")   # member(uid str) -> rep
P_clu = np.zeros((NSK, NTY), int)
pcells = defaultdict(set)
all_clu = set()
for fn in sorted(glob.glob(str(R / "s5_exp_full_g198" / "hitparts" / "*.tsv.zst"))):
    p = subprocess.Popen([ZSTD, "-dc", fn], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        t = line.rstrip("\n").split("\t")
        if len(t) >= 2:
            m = QP.match(t[0])
            if m and t[1].isdigit():
                cl = pdb_clu.get(str(int(t[1])))         # strip zero-pad to match fasta id
                if cl:
                    pcells[(int(m.group(1)), int(m.group(2)))].add(cl)
                    all_clu.add(cl)
    p.wait()
for c, st in pcells.items():
    P_clu[c] = len(st)
NP = len(all_clu)
print(f"PDB 50%-clusters (lit) {NP:,}; PDB-lit cells {int((P_clu > 0).sum()):,}", flush=True)

# depth ratio Q is PER-STRUCTURE (how deeply AFDB is sampled vs PDB), because the
# thinning acts on cluster MEMBER structures (ks below). It is NOT families/clusters
# -- that was an ~80x error that over-thinned everything to non-significance. Same
# per-structure ratio as the fold recompute: PDB lit domains / AFDB lit structures.
N_PDB_DOM = int((G / "pdb_ndistinct_expfull.txt").read_text())   # 826,903 (fold recompute)
N_AFDB_DOM = 1_191_270                                            # AFDB lit structures
Q = N_PDB_DOM / N_AFDB_DOM
assert 0 < Q < 1, f"Q={Q:.3f} not a thinning fraction"
print(f"depth ratio Q(per-structure) = {N_PDB_DOM:,} / {N_AFDB_DOM:,} = {Q:.4f}", flush=True)

# ---- Fisher on cluster counts ------------------------------------------------
cells, pv = [], []
for sk in range(NSK):
    for ty in range(NTY):
        if not SAT[sk, ty]:
            continue
        a, p = int(A_clu[sk, ty]), int(P_clu[sk, ty])
        if a == 0 and p == 0:
            continue
        _, pval = fisher_exact([[a, NA - a], [p, NP - p]])
        cells.append((sk, ty))
        pv.append(pval)
q = multipletests(np.array(pv), method="fdr_bh")[1]

# ---- depth correction: thin AFDB clusters to PDB depth -----------------------
def ncdf(x, mu, sd):
    return 0.5 * (1 + erf((x - mu) / (sd * sqrt(2)))) if sd > 0 else float(x >= mu)


p_depth = np.ones((NSK, NTY))
for (sk, ty) in cells:
    ks = np.fromiter(cell_clu_k[(sk, ty)].values(), float)
    if ks.size == 0:
        continue
    surv = 1.0 - (1.0 - Q) ** ks
    Pv = int(P_clu[sk, ty])
    if Pv == 0:
        p_depth[sk, ty] = 2 * min(float(np.prod(1.0 - surv)), 0.5)   # exact P(X<=0)
    else:
        mu = float(surv.sum())
        sd = sqrt(float((surv * (1 - surv)).sum()))
        p_depth[sk, ty] = 2 * min(ncdf(Pv + 0.5, mu, sd),
                                  1 - ncdf(Pv - 0.5, mu, sd), 0.5)
qd = multipletests(np.array([p_depth[c] for c in cells]), method="fdr_bh")[1]
idx = {c: i for i, c in enumerate(cells)}

# ---- the AFDB-only question, cluster unit ------------------------------------
afdb_only = [c for c in cells if P_clu[c] == 0 and A_clu[c] > 0]
f_sig = [c for c in afdb_only if q[idx[c]] < 0.05]
d_sig = [c for c in afdb_only if qd[idx[c]] < 0.05]
both = sorted(set(f_sig) & set(d_sig))

# compare to the fold-unit recompute
fold = json.load(open(G / "significance_expfull.json"))
print("\n=== cluster unit (50% matched both sides) vs fold unit, on the 1.4M sweep ===")
print(f"{'':32s}{'fold unit':>11}{'cluster unit':>14}")
print(f"{'AFDB-only cells tested':32s}{fold['afdb_only_tested']:>11}{len(afdb_only):>14}")
print(f"{'Fisher q<0.05':32s}{fold['fisher_sig']:>11}{len(f_sig):>14}")
print(f"{'depth-corrected q<0.05':32s}{fold['depth_sig']:>11}{len(d_sig):>14}")
print(f"{'BOTH (defensible)':32s}{fold['both']:>11}{len(both):>14}")

OLD5 = [(16, 12), (171, 7), (158, 27), (34, 0), (33, 7)]
print("\n=== the published 5 depth-significant cells, cluster unit ===")
tracked = []
for (sk, ty) in OLD5:
    pf = int(P_clu[sk, ty])
    i = idx.get((sk, ty))
    qd_v = float(qd[i]) if i is not None else None
    surv = (i is not None and pf == 0 and qd[i] < 0.05)
    stat = "still AFDB-only" if pf == 0 else f"PDB-lit ({pf} clusters)"
    print(f"  sk{sk:03d} ty{ty:02d}: {stat}"
          + (f", depth q={qd_v:.2e} {'SURVIVES' if surv else 'drops'}"
             if qd_v is not None else ""))
    tracked.append(dict(sk=sk, ty=ty, pdb_clusters=pf, depth_q=qd_v, survives=bool(surv)))

json.dump(dict(
    unit_afdb="seq50_cluster", unit_pdb="seq50_cluster_matched", depth_ratio=Q,
    NA=NA, NP=NP, afdb_only=len(afdb_only),
    fisher_sig=len(f_sig), depth_sig=len(d_sig), both=len(both),
    both_cells=[[int(sk), int(ty), int(A_clu[sk, ty]),
                 float(q[idx[(sk, ty)]]), float(qd[idx[(sk, ty)]])] for sk, ty in both],
    old_five=tracked,
), open(G / "cluster_significance_expfull.json", "w"), indent=1)
print(f"\nwrote {G/'cluster_significance_expfull.json'}")
