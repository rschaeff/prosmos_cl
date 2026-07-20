#!/usr/bin/env python3
"""Instrument agreement: does PALSSE read a predicted model the same way it reads
the crystal structure of the SAME molecule?

Every cross-database claim in this project assumes it does. If predicted models --
which have no missing density, no static disorder and idealised termini -- yield a
systematically different SSE inventory, then "this domain lights cell X" is partly
a statement about provenance, and the AFDB-only cells are the first place that
contamination would show up.

Four measurements, in increasing order of how much they matter:

  1. SSE inventory     -- does PALSSE call the same number of helices/strands?
  2. Cell agreement    -- do the two halves of a pair light the same cells?
  3. Directional bias  -- when they disagree, does the PREDICTED side light more?
                          (idealisation would inflate the predicted side)
  4. The one that matters: are the 619 AFDB-only cells enriched for
     predicted-side-only lighting? A cell that is "AFDB-only" in the big sweep AND
     preferentially lit by the predicted half of a same-molecule pair is a
     detector artifact, not biology.

(4) is the falsification test. (1)-(3) are how we interpret it.
"""
import json, re, subprocess, sys
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np

R = Path("/home/rschaeff/work/prosmos_2026/instr_agree")
G = Path("/home/rschaeff/work/prosmos_2026/s5_grid")
ZSTD = "/sw/apps/Anaconda3-2023.09-0/bin/zstd"

# The fixed-width metamatrix format collides fields (HA1011, 14.176-119.171), so a
# \s+ regex silently mis-parses; this is the corrected form (chain at group 2,
# length at group 5) validated against declared SSE counts.
SSE = re.compile(
    r'([HE])([A-Za-z0-9])\s*(-?\d+[A-Za-z]?)\s*--\s*(-?\d+[A-Za-z]?)\s+(\d+)\s+-?\d+\.\d')
REC = re.compile(r'^(\S+)\.ssd')
QP = re.compile(r'^s5-(\d{4})-(\d{4})$')

pairs = json.load(open(R / "pairs.json"))
print(f"pairs: {len(pairs):,}")

# record id == basename of the domain .pdb (the ECOD store's numeric uid)
exp_rec = {Path(p["exp_pdb"]).stem: p for p in pairs}
pred_rec = {Path(p["pred_pdb"]).stem: p for p in pairs}
for p in pairs:
    p["exp_rec"] = Path(p["exp_pdb"]).stem
    p["pred_rec"] = Path(p["pred_pdb"]).stem


def sse_counts(side):
    """record -> (n_helix, n_strand), concatenated over all shard DBs."""
    out, cur = {}, None
    for dbpath in sorted((R / "dbs").glob(f"{side}_*.db")):
        with open(dbpath, errors="replace") as fh:
            for line in fh:
                m = REC.match(line)
                if m:
                    cur = m.group(1)
                    out.setdefault(cur, [0, 0])
                    continue
                if cur:
                    s = SSE.match(line.strip())
                    if s:
                        out[cur][0 if s.group(1) == "H" else 1] += 1
    return {k: tuple(v) for k, v in out.items()}


def cells(side):
    """record -> set of (skeleton, typing) cells lit, concatenated over all shards."""
    d = defaultdict(set)
    parts = sorted((R / "hitparts").glob(f"{side}_*.tsv.zst"))
    if not parts:
        raise SystemExit(f"no hitparts for {side} -- did the array finish?")
    for tsv in parts:
        p = subprocess.Popen([ZSTD, "-dc", str(tsv)], stdout=subprocess.PIPE, text=True)
        for line in p.stdout:
            t = line.rstrip("\n").split("\t")
            if len(t) < 2:
                continue
            m = QP.match(t[0])
            if m:
                d[t[1]].add((int(m.group(1)), int(m.group(2))))
        p.wait()
    return d


print("\n=== 1. SSE INVENTORY ===")
se, sp = sse_counts("exp"), sse_counts("pred")
print(f"records in exp.db {len(se):,}   pred.db {len(sp):,}")
both_sse = [(p, se[p["exp_rec"]], sp[p["pred_rec"]])
            for p in pairs if p["exp_rec"] in se and p["pred_rec"] in sp]
print(f"pairs with both sides parsed: {len(both_sse):,}")
dh = np.array([b[2][0] - b[1][0] for b in both_sse])
ds = np.array([b[2][1] - b[1][1] for b in both_sse])
dt = dh + ds
print(f"  helices : pred - exp  mean {dh.mean():+.3f}  median {np.median(dh):+.0f}  "
      f"identical {100*(dh==0).mean():.1f}%")
print(f"  strands : pred - exp  mean {ds.mean():+.3f}  median {np.median(ds):+.0f}  "
      f"identical {100*(ds==0).mean():.1f}%")
print(f"  total   : pred - exp  mean {dt.mean():+.3f}  identical inventory {100*(dt==0).mean():.1f}%")
print(f"  pred richer {100*(dt>0).mean():.1f}%   exp richer {100*(dt<0).mean():.1f}%")

print("\n=== 2. CELL AGREEMENT ===")
ce, cp = cells("exp"), cells("pred")
print(f"records with >=1 hit: exp {len(ce):,}  pred {len(cp):,}")
rows = []
for p in pairs:
    a, b = ce.get(p["exp_rec"], set()), cp.get(p["pred_rec"], set())
    if not a and not b:
        continue                       # neither side lights anything: uninformative
    rows.append((p, a, b))
print(f"pairs where at least one side lights a cell: {len(rows):,}")
exact = sum(1 for _, a, b in rows if a == b)
jac = np.array([len(a & b) / len(a | b) for _, a, b in rows])
both_lit = [(p, a, b) for p, a, b in rows if a and b]
jac_bl = np.array([len(a & b) / len(a | b) for _, a, b in both_lit])
print(f"  identical cell set        : {exact:,}/{len(rows):,} ({100*exact/len(rows):.1f}%)")
print(f"  mean Jaccard (all)        : {jac.mean():.3f}   median {np.median(jac):.3f}")
print(f"  mean Jaccard (both lit)   : {jac_bl.mean():.3f}   median {np.median(jac_bl):.3f}  n={len(both_lit):,}")
onlyexp = sum(1 for _, a, b in rows if a and not b)
onlypred = sum(1 for _, a, b in rows if b and not a)
print(f"  exp lights, pred dark     : {onlyexp:,}")
print(f"  pred lights, exp dark     : {onlypred:,}")

print("\n=== 3. DIRECTIONAL BIAS ===")
ne = np.array([len(a) for _, a, b in rows])
np_ = np.array([len(b) for _, a, b in rows])
print(f"  cells lit per domain: exp {ne.mean():.2f}   pred {np_.mean():.2f}   "
      f"ratio {np_.mean()/max(ne.mean(),1e-9):.3f}")
d = np_ - ne
print(f"  pred lights MORE cells in {100*(d>0).mean():.1f}% of pairs, "
      f"FEWER in {100*(d<0).mean():.1f}%, equal in {100*(d==0).mean():.1f}%")
# sign test on the discordant pairs -- a fair detector gives a 50/50 split
nplus, nminus = int((d > 0).sum()), int((d < 0).sum())
nd = nplus + nminus
if nd:
    z = (nplus - nd / 2) / np.sqrt(nd / 4)
    print(f"  sign test on {nd:,} discordant pairs: {nplus:,} pred-richer vs "
          f"{nminus:,} exp-richer, z = {z:+.2f}")

print("\n=== 4. ARE AFDB-ONLY CELLS DETECTOR ARTIFACTS? ===")
aonly = np.load(G / "afdb_only_mask.npy")
imp = np.load(G / "impossible_mask.npy")
sat = ~imp
print(f"AFDB-only cells from the full sweep: {int(aonly.sum()):,}")

# per cell, across all pairs: how often lit by pred-only vs exp-only
pred_only_c, exp_only_c, both_c = Counter(), Counter(), Counter()
for _, a, b in rows:
    for c in b - a:
        pred_only_c[c] += 1
    for c in a - b:
        exp_only_c[c] += 1
    for c in a & b:
        both_c[c] += 1

touched = set(pred_only_c) | set(exp_only_c) | set(both_c)
in_ao = [c for c in touched if aonly[c]]
not_ao = [c for c in touched if sat[c] and not aonly[c]]
print(f"cells this paired set touches at all: {len(touched):,}  "
      f"(of which AFDB-only: {len(in_ao):,})")


def bias(cs):
    """Fraction of provenance-discordant observations that are predicted-side-only."""
    po = sum(pred_only_c[c] for c in cs)
    eo = sum(exp_only_c[c] for c in cs)
    return po, eo, (po / (po + eo) if (po + eo) else float("nan"))


po_a, eo_a, f_a = bias(in_ao)
po_n, eo_n, f_n = bias(not_ao)
print(f"  AFDB-only cells   : pred-only {po_a:,}  exp-only {eo_a:,}  "
      f"predicted-side share {f_a:.3f}")
print(f"  all other cells   : pred-only {po_n:,}  exp-only {eo_n:,}  "
      f"predicted-side share {f_n:.3f}")
if not (np.isnan(f_a) or np.isnan(f_n)):
    # two-proportion z-test; 0.5 is the fair-detector null for either group
    p1, n1 = f_a, po_a + eo_a
    p2, n2 = f_n, po_n + eo_n
    pp = (po_a + po_n) / (n1 + n2)
    zz = (p1 - p2) / np.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
    print(f"  difference between them: z = {zz:+.2f}")
    print(f"  AFDB-only vs fair-detector null (0.5): "
          f"z = {(p1-0.5)/np.sqrt(0.25/n1):+.2f}")

json.dump({
    "n_pairs": len(pairs),
    "n_informative": len(rows),
    "exact_match": exact,
    "exact_match_frac": exact / len(rows),
    "jaccard_mean": float(jac.mean()),
    "jaccard_mean_bothlit": float(jac_bl.mean()),
    "sse_identical_frac": float((dt == 0).mean()),
    "sse_delta_mean": float(dt.mean()),
    "cells_per_domain": {"exp": float(ne.mean()), "pred": float(np_.mean())},
    "afdb_only_pred_share": None if np.isnan(f_a) else float(f_a),
    "other_pred_share": None if np.isnan(f_n) else float(f_n),
    "afdb_only_counts": {"pred_only": po_a, "exp_only": eo_a},
    "other_counts": {"pred_only": po_n, "exp_only": eo_n},
}, open(R / "agreement.json", "w"), indent=1)
print(f"\nwrote {R/'agreement.json'}")
