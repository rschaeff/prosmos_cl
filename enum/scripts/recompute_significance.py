#!/usr/bin/env python3
"""Recompute AFDB-only S5 significance against the corrected 1.4M experimental sweep.

WHY THIS SUPERSEDES fisher_cells.py + the published significant set. The published
result compared AFDB against a PDB occupancy grid built from the OLD experimental
sweep, which had two defects now fixed:

  DEPTH. The old sweep covered 496k experimental domains = 33.5% of ECOD's
  experimental set. The other two thirds had no structure file, so "PDB never
  samples this cell" often meant "absent from the searchable third". The
  regenerated set covers 1,400,111 domains (93.8%). Every PDB-zero cell must be
  retested against this depth.

  CORRECTNESS. The old sweep ran a searchmatrix whose SSE loop was bounded by line
  length rather than numberElment, injecting phantom SSEs from records whose
  coordinates overflow the 8-char field, and silently losing four crashed chunks.
  Fixed in 1e5546c; this reads s5_exp_full_g198, produced by the fixed binary.

QUERY SET MUST MATCH THE BASELINE. The published AFDB grid, depth-ratio constant,
and afdb_rare inputs all come from queries_graph198_alltypings (fully typed:
v/C/T/t/c interaction codes). It is NOT queries_adjonly, which is orientation-blind
(x only) and matches ~20x more structures. Sweeping the new DB with adjonly and
comparing to the typed AFDB grid inflated one cell from 0 to 229,424 domains. The
_check_query_set guard below refuses to run on the wrong sweep.

The AFDB grid is unchanged: the AFDB DB has zero overflow records (AlphaFold
models are origin-centred), so that sweep was never affected.

METHOD (identical to the published fisher_cells.py, so the two are comparable):
  (1) Fisher on FOLD counts per cell -- fold (ECOD T-group) is the unit, because
      counting domains would treat PDB's ~30x redundancy as independent evidence.
  (2) DEPTH-CORRECTED. Fisher's margins encode a sampling-depth difference that is
      not biological. Thin AFDB to PDB's per-structure depth q and ask whether PDB
      still lights fewer folds than expected: fold f with k_f member structures
      survives thinning with prob 1-(1-q)^k_f, so the thinned fold count is
      Poisson-binomial over that cell's folds. A cell is depth-robust only if PDB's
      observed fold count is extreme against that distribution. Both BH/FDR.

THE HEADLINE is how many of the old 619 PDB-zero / AFDB>0 candidate cells still
have zero PDB folds now, and how many survive depth correction -- and in
particular whether the two cells known to be contaminated by experimental domains
that lacked structure files at the time (sk016 ty12 and sk033 ty07)
survive now that those structures exist and were searched.

Outputs (all under s5_grid/):
  grid_pdb_nT_expfull.npy         rebuilt PDB fold grid
  grid_pdb_ndom_expfull.npy       rebuilt PDB domain grid (for depth ratio)
  significance_expfull.npz        per-cell Fisher + depth q-values
  significance_expfull.json       summary + the surviving cells + contamination check
"""
from __future__ import annotations

import argparse
import glob
import json
import pickle
import re
import subprocess
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

R = Path("/home/rschaeff/work/prosmos_2026")
G = R / "s5_grid"
ZSTD = "/sw/apps/Anaconda3-2023.09-0/bin/zstd"
QP = re.compile(r"^s5-(\d{4})-(\d{4})$")
NSK, NTY = 198, 32

# cells that carried an ecod_* (experimental) domain in the AFDB sweep despite
# reading PDB-zero in the old grid -- the false positives this recompute must
# resolve. From the /ruczinski contamination check.
CONTAMINATED = {(16, 12), (33, 7)}


def load_uid_class():
    """uid -> (t_group, f_group). Prefer the fresh full pull if present, else the
    build-time table (99.5% coverage of the new DB)."""
    fresh = G / "uid_class_expfull.tsv"
    src = fresh if fresh.exists() else (R / "pdb_exp_build" / "uid_class.tsv")
    uid_t, uid_f = {}, {}
    for line in open(src):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and f[0].isdigit():
            if f[2]:
                uid_t[int(f[0])] = f[2]
            if f[3]:
                uid_f[int(f[0])] = f[3]
    print(f"uid_class from {src.name}: uid->t {len(uid_t):,}  uid->f {len(uid_f):,}",
          flush=True)
    return uid_t, uid_f


def _scan(fn):
    """(sk,ty) -> set(uid) for one hitpart. Experimental DB: every record is a
    numeric ecod_uid, so no provenance split is needed here."""
    d = defaultdict(set)
    p = subprocess.Popen([ZSTD, "-dc", fn], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        t = line.rstrip("\n").split("\t")
        if len(t) >= 2:
            m = QP.match(t[0])
            if m and t[1].isdigit():
                d[(int(m.group(1)), int(m.group(2)))].add(int(t[1]))
    p.wait()
    return d


def build_pdb_grid(uid_t, uid_f, sweep="s5_exp_full_g198", nproc=16):
    """Rebuild PDB fold / family / domain grids from the corrected sweep."""
    parts = sorted(glob.glob(str(R / sweep / "hitparts" / "*.tsv.zst")))
    print(f"scanning {len(parts)} hitparts of {sweep} ...", flush=True)
    cell_uids = defaultdict(set)
    with Pool(nproc) as pool:
        for d in pool.imap_unordered(_scan, parts, chunksize=4):
            for c, us in d.items():
                cell_uids[c] |= us
    nT = np.zeros((NSK, NTY), int)
    nDom = np.zeros((NSK, NTY), int)
    all_lit = set()
    for (sk, ty), us in cell_uids.items():
        nDom[sk, ty] = len(us)
        nT[sk, ty] = len({uid_t[u] for u in us if u in uid_t})
        all_lit |= us
    n_distinct = len(all_lit)                 # DISTINCT domains lighting >=1 cell
    print(f"cells with PDB hits: {len(cell_uids):,}  distinct exp domains lit: "
          f"{n_distinct:,}  (sum-over-cells {int(nDom.sum()):,} double-counts)",
          flush=True)
    return nT, nDom, n_distinct


def pois_binom_cdf(ps, x):
    """P(X <= x) for independent Bernoulli(ps) via exact DP."""
    dist = np.zeros(len(ps) + 1)
    dist[0] = 1.0
    for pr in ps:
        dist[1:] = dist[1:] * (1 - pr) + dist[:-1] * pr
        dist[0] *= (1 - pr)
    return dist[: x + 1].sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", default="s5_exp_full_g198")
    ap.add_argument("--rebuild", action="store_true",
                    help="rebuild the PDB grid even if cached")
    a = ap.parse_args()

    A = np.load(G / "grid_afdb_nT_rebuilt.npy").astype(int)   # AFDB folds (clean)
    IMP = np.load(G / "impossible_mask.npy")
    SAT = ~IMP
    NA_FOLDS = 3656                                           # AFDB fold universe
    NP_FOLDS = 3935                                           # PDB  fold universe

    # ---- new PDB grid from the corrected sweep
    cache = G / "grid_pdb_nT_expfull.npy"
    if cache.exists() and not a.rebuild:
        P = np.load(cache).astype(int)
        Pd = np.load(G / "grid_pdb_ndom_expfull.npy").astype(int)
        n_distinct = int((G / "pdb_ndistinct_expfull.txt").read_text())
        print(f"loaded cached PDB grid ({int((P > 0).sum())} folded cells)")
    else:
        uid_t, uid_f = load_uid_class()
        P, Pd, n_distinct = build_pdb_grid(uid_t, uid_f, a.sweep)
        np.save(cache, P)
        np.save(G / "grid_pdb_ndom_expfull.npy", Pd)
        (G / "pdb_ndistinct_expfull.txt").write_text(str(n_distinct))

    # ---- depth ratio, recomputed for the new experimental depth
    N_PDB_DOM = n_distinct               # DISTINCT experimental domains lighting any cell
    N_AFDB_DOM = 1_191_270              # DISTINCT AFDB structures lighting any cell
    Q = N_PDB_DOM / N_AFDB_DOM
    assert 0 < Q < 1, (f"depth ratio Q={Q:.3f} is not a thinning fraction; PDB lit "
                       f"({N_PDB_DOM:,}) should be < AFDB lit ({N_AFDB_DOM:,}). "
                       f"Wrong query set? (adjonly over-matches)")
    print(f"depth ratio Q = {N_PDB_DOM:,} / {N_AFDB_DOM:,} = {Q:.4f}", flush=True)

    # ---- Fisher on fold counts
    cells, pv = [], []
    for sk in range(NSK):
        for ty in range(NTY):
            if not SAT[sk, ty]:
                continue
            af, pf = int(A[sk, ty]), int(P[sk, ty])
            if af == 0 and pf == 0:
                continue
            _, pval = fisher_exact([[af, NA_FOLDS - af], [pf, NP_FOLDS - pf]])
            cells.append((sk, ty))
            pv.append(pval)
    pv = np.array(pv)
    q = multipletests(pv, method="fdr_bh")[1]

    # ---- depth correction: AFDB thinned to PDB depth
    tg = {}
    for line in open(G / "afdb_name_tgroup.tsv"):
        n, t = line.rstrip("\n").split("\t")
        tg[n] = t
    pairs = []
    for f in glob.glob(str(G / "afdb_rare" / "*.pkl")):
        pairs.extend(pickle.load(open(f, "rb")))
    cell_fold_k = defaultdict(Counter)
    for cid, rec in pairs:
        t = tg.get(rec)
        if t:
            cell_fold_k[cid][t] += 1

    p_depth = np.ones((NSK, NTY))
    exp_thin = np.zeros((NSK, NTY))
    for (sk, ty) in cells:
        ks = np.array(list(cell_fold_k[sk * NTY + ty].values()), dtype=float)
        if ks.size == 0:
            continue
        surv = 1.0 - (1.0 - Q) ** ks
        exp_thin[sk, ty] = surv.sum()
        lo = pois_binom_cdf(surv, int(P[sk, ty]))
        p_depth[sk, ty] = 2 * min(lo, 0.5)         # one-sided low, two-sided guard
    pd_list = np.array([p_depth[sk, ty] for sk, ty in cells])
    qd = multipletests(pd_list, method="fdr_bh")[1]

    idx = {c: i for i, c in enumerate(cells)}

    # ---- the AFDB-only question, before vs after
    P_old = np.load(G / "grid_pdb_nT.npy").astype(int)
    old_pool = [(sk, ty) for sk in range(NSK) for ty in range(NTY)
                if SAT[sk, ty] and P_old[sk, ty] == 0 and A[sk, ty] > 0]
    now_zero = [c for c in old_pool if P[c] == 0]
    now_lit = [c for c in old_pool if P[c] > 0]

    afdb_only = [(sk, ty) for (sk, ty) in cells if P[sk, ty] == 0 and A[sk, ty] > 0]
    fisher_sig = [c for c in afdb_only if q[idx[c]] < 0.05]
    depth_sig = [c for c in afdb_only if qd[idx[c]] < 0.05]
    both = sorted(set(fisher_sig) & set(depth_sig))

    print("\n=== AFDB-only recompute (fold unit) ===")
    print(f"old PDB-zero / AFDB>0 candidate cells:        {len(old_pool)}")
    print(f"  now LIT by experimental (>=1 PDB fold):     {len(now_lit)}  "
          f"({100*len(now_lit)/len(old_pool):.0f}% -- newly sampled)")
    print(f"  still PDB-zero:                             {len(now_zero)}")
    print(f"still-AFDB-only cells tested:                 {len(afdb_only)}")
    print(f"  Fisher q<0.05:                              {len(fisher_sig)}")
    print(f"  depth-corrected q<0.05:                     {len(depth_sig)}")
    print(f"  BOTH (defensible):                          {len(both)}")

    # ---- the contamination check, first-class
    print("\n=== contaminated depth-significant cells (had exp domains lacking files) ===")
    contam = []
    for (sk, ty) in sorted(CONTAMINATED):
        pf = int(P[sk, ty])
        verdict = ("RESOLVED: now PDB-lit, no longer AFDB-only" if pf > 0
                   else "still PDB-zero")
        i = idx.get((sk, ty))
        qd_v = float(qd[i]) if i is not None else None
        print(f"  sk{sk:03d} ty{ty:02d}: PDB folds now = {pf}  -> {verdict}"
              + (f"  (depth q={qd_v:.2e})" if qd_v is not None else ""))
        contam.append(dict(sk=sk, ty=ty, pdb_folds_now=pf, resolved=pf > 0,
                           depth_q=qd_v))

    # ---- old depth-significant five, tracked
    print("\n=== the published 5 depth-significant cells, now ===")
    OLD5 = [(16, 12), (171, 7), (158, 27), (34, 0), (33, 7)]
    tracked = []
    for (sk, ty) in OLD5:
        pf = int(P[sk, ty])
        i = idx.get((sk, ty))
        stat = ("still AFDB-only" if pf == 0 else f"now PDB-lit ({pf} folds)")
        qd_v = float(qd[i]) if i is not None else None
        surv = (i is not None and pf == 0 and qd[i] < 0.05)
        print(f"  sk{sk:03d} ty{ty:02d}: {stat}"
              + (f", depth q={qd_v:.2e} {'SURVIVES' if surv else 'drops'}"
                 if qd_v is not None else ""))
        tracked.append(dict(sk=sk, ty=ty, pdb_folds_now=pf, depth_q=qd_v,
                            survives=bool(surv)))

    np.savez(G / "significance_expfull.npz",
             cells=np.array(cells), pv=pv, q=q, p_depth=pd_list, q_depth=qd,
             exp_thin=np.array([exp_thin[sk, ty] for sk, ty in cells]))
    json.dump(dict(
        depth_ratio=Q, n_pdb_dom=N_PDB_DOM,
        old_candidate_pool=len(old_pool), now_lit=len(now_lit),
        still_zero=len(now_zero), afdb_only_tested=len(afdb_only),
        fisher_sig=len(fisher_sig), depth_sig=len(depth_sig), both=len(both),
        both_cells=[[int(sk), int(ty), int(A[sk, ty]),
                     float(q[idx[(sk, ty)]]), float(qd[idx[(sk, ty)]])]
                    for (sk, ty) in both],
        contamination=contam, old_five=tracked,
    ), open(G / "significance_expfull.json", "w"), indent=1)
    print(f"\nwrote {G/'significance_expfull.npz'} + significance_expfull.json")


if __name__ == "__main__":
    main()
