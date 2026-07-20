#!/usr/bin/env python3
"""The n=3/4/5 curve: why S5, made empirical.

Three quantities, each as a function of motif size n, computed on the SAME record
subset at every n (afdb_sub_records / pdb_sub_records) so the trend is about n and
not about which records were searched:

  1. saturation  -- fraction of enumerated cells that are lit. Is topology space
                    filling up as n grows?
  2. sharing     -- mean distinct ECOD T-groups per lit cell. This is the load-
                    bearing one: the largest n at which a motif is still SHARED
                    across unrelated folds is where the abstraction still
                    abstracts. Beyond that a cell is a per-domain fingerprint.
  3. agreement   -- Jaccard of the AFDB-lit and PDB-lit cell sets. Flat across n
                    => the "same alphabet" claim is scale-robust; narrowing =>
                    divergence begins, which would locate the vocabulary boundary
                    and justify the ~24x cost of n=6.

n=3,4 come from the dedicated sweeps (hitparts); n=5 AFDB from the full sweep's
per-cell rep dump and n=5 PDB from its hitparts -- both restricted to the subset.
"""
import subprocess, re, glob, pickle
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
import numpy as np

R = Path("/home/rschaeff/work/prosmos_2026")
G = R / "s5_grid"
ZSTD = "/sw/apps/Anaconda3-2023.09-0/bin/zstd"
QP = re.compile(r'^s([345])-(\d{4})-(\d{4})$')
TOTAL = {3: 5 * 8, 4: 42 * 16, 5: 198 * 32}

afdb_sub = pickle.load(open(G / "afdb_sub_records.pkl", "rb"))
pdb_sub = pickle.load(open(G / "pdb_sub_records.pkl", "rb"))
print(f"sub-records: AFDB {len(afdb_sub):,}  PDB {len(pdb_sub):,}", flush=True)

# ---- T-group maps
afdb_tg = {}
for line in open(G / "afdb_name_tgroup.tsv"):
    a, b = line.rstrip("\n").split("\t")
    afdb_tg[a] = b            # keyed by full record name (accession_Dn)
pdb_tg = {}
for line in open(R / "pdb_exp_build" / "uid_class.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 3 and f[0].isdigit() and f[2]:
        pdb_tg[f[0]] = f[2]
print(f"tgroup maps: AFDB {len(afdb_tg):,}  PDB {len(pdb_tg):,}", flush=True)


def scan_hitparts(pattern, sub):
    """cell(sk,ty) -> set(record), restricted to sub."""
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
    """cell -> set(record) for n=5 from the full sweep's rep dump, restricted."""
    out = defaultdict(set)
    for line in open(R / "design_scaffold" / "afdb_cell_reps.tsv"):
        sk, ty, acc = line.rstrip("\n").split("\t")
        if acc in sub:
            out[(int(sk), int(ty))].add(acc)
    return out


def metrics(afdb_cells, pdb_cells, n):
    def sharing(cells, tg, norm=lambda r: r):
        vals = []
        for c, recs in cells.items():
            folds = {tg[k] for r in recs if (k := norm(r)) in tg}
            if folds:
                vals.append(len(folds))
        return vals
    a_share = sharing(afdb_cells, afdb_tg)
    # PDB hitpart uids are zero-padded ('000432069'); the T-group map is keyed
    # unpadded ('432069'), so normalise via int before lookup.
    p_share = sharing(pdb_cells, pdb_tg, norm=lambda r: str(int(r)))
    a_lit = {c for c, s in afdb_cells.items() if s}
    p_lit = {c for c, s in pdb_cells.items() if s}
    union = a_lit | p_lit
    inter = a_lit & p_lit
    jac = len(inter) / len(union) if union else float("nan")
    return {
        "n": n, "total_cells": TOTAL[n],
        "afdb_lit": len(a_lit), "pdb_lit": len(p_lit),
        "afdb_sat": len(a_lit) / TOTAL[n], "pdb_sat": len(p_lit) / TOTAL[n],
        "afdb_share_mean": float(np.mean(a_share)) if a_share else 0.0,
        "afdb_share_med": float(np.median(a_share)) if a_share else 0.0,
        "pdb_share_mean": float(np.mean(p_share)) if p_share else 0.0,
        "pdb_share_med": float(np.median(p_share)) if p_share else 0.0,
        "jaccard": jac, "shared": len(inter), "afdb_only": len(a_lit - p_lit),
        "pdb_only": len(p_lit - a_lit),
    }


rows = []
for n in (3, 4, 5):
    print(f"\n--- n={n} ---", flush=True)
    if n < 5:
        ac = scan_hitparts(str(R / f"n{n}_afdb" / "hitparts" / "*.tsv.zst"), afdb_sub)
        pc = scan_hitparts(str(R / f"n{n}_pdb" / "hitparts" / "*.tsv.zst"), pdb_sub)
    else:
        ac = afdb_n5(afdb_sub)
        pc = scan_hitparts(str(R / "s5_pdb_inv" / "hitparts" / "*.tsv.zst"), pdb_sub)
    na = len(set().union(*ac.values())) if ac else 0
    npd = len(set().union(*pc.values())) if pc else 0
    print(f"  records hit (in subset): AFDB {na:,}  PDB {npd:,}", flush=True)
    m = metrics(ac, pc, n)
    rows.append(m)
    print(f"  saturation : AFDB {m['afdb_sat']:.1%} ({m['afdb_lit']}/{m['total_cells']})   "
          f"PDB {m['pdb_sat']:.1%} ({m['pdb_lit']})")
    print(f"  sharing    : AFDB mean {m['afdb_share_mean']:.1f} T-groups/cell (med {m['afdb_share_med']:.0f})   "
          f"PDB mean {m['pdb_share_mean']:.1f} (med {m['pdb_share_med']:.0f})")
    print(f"  agreement  : Jaccard {m['jaccard']:.3f}  shared {m['shared']}  "
          f"AFDB-only {m['afdb_only']}  PDB-only {m['pdb_only']}")

print("\n================ CURVE ================")
print(f"{'n':>2} {'cells':>6} {'AFDBsat':>8} {'PDBsat':>7} {'AFDBshare':>9} {'PDBshare':>8} {'Jaccard':>8}")
for m in rows:
    print(f"{m['n']:>2} {m['total_cells']:>6} {m['afdb_sat']:>7.1%} {m['pdb_sat']:>6.1%} "
          f"{m['afdb_share_mean']:>9.1f} {m['pdb_share_mean']:>8.1f} {m['jaccard']:>8.3f}")

import json
json.dump(rows, open(G / "curve_n345.json", "w"), indent=1)
print(f"\nwrote {G/'curve_n345.json'}")
