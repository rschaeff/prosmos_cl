#!/usr/bin/env python3
"""Per-cell PDB-vs-AFDB significance, two ways, so the difference is visible.

(1) FISHER (what was asked for). Per satisfiable cell, a 2x2 on FOLD counts:
        [[a, NA-a], [p, NP-p]]   a = AFDB folds lighting the cell, NA = 3,656
                                 p = PDB  folds lighting the cell, NP = 3,935
    Fold-level is already the right unit: counting structures would treat the
    ~30x redundancy within a T-group as independent evidence. BH-corrected.

(2) DEPTH-CORRECTED. Fisher's margins encode a sampling-depth difference that is
    not biological: AFDB searched 4.92M structures vs PDB 496k, so each fold is
    sampled ~10x deeper and gets ~10x more chances to be SEEN lighting a cell.
    We therefore ask a different question: if AFDB were thinned to PDB's depth,
    how many folds would still light this cell?
      fold f with k_f member structures in the cell survives thinning to fraction
      q with probability 1-(1-q)^k_f  -> the thinned fold count is Poisson-binomial
      over the folds of that cell. Exact, no resampling.
    A cell is depth-robust only if PDB's count is extreme against THAT distribution.

Outputs fisher_cells.npz + a 3-panel figure.
"""
import numpy as np, glob, pickle, json
from collections import defaultdict, Counter
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

G = "/home/rschaeff/work/prosmos_2026/s5_grid"
A = np.load(f"{G}/grid_afdb_nT_rebuilt.npy").astype(int)
P = np.load(f"{G}/grid_pdb_nT.npy").astype(int)
IMP = np.load(f"{G}/impossible_mask.npy")
NA, NP = 3656, 3935
SAT = ~IMP

# ---------- (1) Fisher on fold counts
cells, pv, odds = [], [], []
for sk in range(198):
    for ty in range(32):
        if not SAT[sk, ty]:
            continue
        a, p = int(A[sk, ty]), int(P[sk, ty])
        if a == 0 and p == 0:
            continue                      # no data either way
        orr, pval = fisher_exact([[a, NA - a], [p, NP - p]])
        cells.append((sk, ty)); pv.append(pval); odds.append(orr)
pv = np.array(pv)
q = multipletests(pv, method="fdr_bh")[1]
print(f"cells tested (satisfiable, lit in >=1 DB): {len(cells):,}")
print(f"  Fisher q<0.05: {(q < 0.05).sum():,}")

# ---------- (2) depth-corrected: AFDB thinned to PDB depth
tg = {}
for line in open(f"{G}/afdb_name_tgroup.tsv"):
    n, t = line.rstrip("\n").split("\t"); tg[n] = t
pairs = []
for f in glob.glob(f"{G}/afdb_rare/*.pkl"):
    pairs.extend(pickle.load(open(f, "rb")))
cell_fold_k = defaultdict(Counter)          # cell -> fold -> #member structures
for cid, r in pairs:
    t = tg.get(r)
    if t:
        cell_fold_k[cid][t] += 1
# depth ratio: AFDB lit structures vs PDB lit structures (the like-for-like unit)
Q = 295_360 / 1_191_270
print(f"  thinning AFDB to q={Q:.3f} (PDB lit / AFDB lit structures)")


def pois_binom_cdf(ps, x):
    """P(X <= x) for independent Bernoulli(ps) via DP."""
    dist = np.zeros(len(ps) + 1); dist[0] = 1.0
    for pr in ps:
        dist[1:] = dist[1:] * (1 - pr) + dist[:-1] * pr
        dist[0] *= (1 - pr)
    return dist[:x + 1].sum()


exp_thin = np.zeros_like(A, float)
p_depth = np.ones_like(A, float)
for (sk, ty) in cells:
    ks = np.array(list(cell_fold_k[sk * 32 + ty].values()), dtype=float)
    if ks.size == 0:
        continue
    surv = 1.0 - (1.0 - Q) ** ks            # P(fold still seen at PDB depth)
    exp_thin[sk, ty] = surv.sum()
    # is PDB's count low against the thinned-AFDB distribution? (one-sided both ways)
    lo = pois_binom_cdf(surv, int(P[sk, ty]))
    hi = 1.0 - pois_binom_cdf(surv, int(P[sk, ty]) - 1) if P[sk, ty] > 0 else 1.0
    p_depth[sk, ty] = 2 * min(lo, hi, 0.5)  # two-sided
pd_list = np.array([p_depth[sk, ty] for sk, ty in cells])
qd = multipletests(pd_list, method="fdr_bh")[1]
print(f"  depth-corrected q<0.05: {(qd < 0.05).sum():,}")

# ---------- the cells Qian wants: one DB zero, the other significant
res = {}
for tag, cond in [("PDB0_AFDBsig", lambda a, p: p == 0 and a > 0),
                  ("AFDB0_PDBsig", lambda a, p: a == 0 and p > 0)]:
    idx = [i for i, (sk, ty) in enumerate(cells) if cond(A[sk, ty], P[sk, ty])]
    f_sig = [i for i in idx if q[i] < 0.05]
    d_sig = [i for i in idx if qd[i] < 0.05]
    both = sorted(set(f_sig) & set(d_sig))
    res[tag] = dict(n=len(idx), fisher_sig=len(f_sig), depth_sig=len(d_sig), both=len(both),
                    cells=[[int(cells[i][0]), int(cells[i][1]), int(A[cells[i]]), int(P[cells[i]]),
                            float(q[i]), float(qd[i])] for i in both])
    print(f"\n{tag}: {len(idx):,} cells | Fisher q<0.05: {len(f_sig):,} | "
          f"depth-corrected q<0.05: {len(d_sig):,} | BOTH: {len(both):,}")

np.savez(f"{G}/fisher_cells.npz",
         cells=np.array(cells), pv=pv, q=q, odds=np.array(odds),
         p_depth=pd_list, q_depth=qd, exp_thin=exp_thin)
json.dump(res, open(f"{G}/fisher_zero_cells.json", "w"), indent=1)
print("\nwrote fisher_cells.npz + fisher_zero_cells.json")
