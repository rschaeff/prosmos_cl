#!/usr/bin/env python3
"""Per-cell lineage composition of the AFDB S5 grid.

Source: afdb_200m.protein (214.7M rows, uniprot_acc -> organism_tax_id), joined to
afdb_200m.taxonomy for superkingdom. 99% coverage and contemporaneous with the
models we searched -- NOT current UniProt, which resolves only 34% of our
accessions and whose survivors are 51% eukaryotic against a true 27% (a
survivorship bias that would manufacture the expected eukaryote signal).

Question: RS expects the all-helix (HHHHH) cells to be eukaryote-enriched because
that is where repeat proteins (ARM/HEAT/ankyrin) live. Is that the ONLY signal?

Two passes:
 (1) raw   -- eukaryote fraction per cell vs the global fraction, binomial, BH.
 (2) fold-controlled -- the raw signal is mostly mediated by WHICH FOLDS occupy a
     cell (ARM repeats are eukaryotic and light HHHHH). So we also ask, within each
     cell, whether its lineage mix differs from what its own fold composition
     predicts: expected euk = sum over folds of (fold's global euk rate x n in cell).
     A cell that is still skewed after that is skewed for a reason OTHER than
     "eukaryotic folds live here".
"""
import numpy as np, glob, pickle, json
from collections import defaultdict, Counter
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

G = "/home/rschaeff/work/prosmos_2026/s5_grid"

acc_sk = {}
for line in open(f"{G}/acc_tax.tsv"):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 3 and p[2]:
        acc_sk[p[0]] = p[2]
print(f"accessions with a superkingdom: {len(acc_sk):,}")

tg = {}
for line in open(f"{G}/afdb_name_tgroup.tsv"):
    n, t = line.rstrip("\n").split("\t"); tg[n] = t

pairs = []
for f in glob.glob(f"{G}/afdb_rare/*.pkl"):
    pairs.extend(pickle.load(open(f, "rb")))

def acc_of(rec):
    return rec.split("_D")[0].replace("dpam_", "")

# global composition + per-fold composition (over DISTINCT records, not cell-hits,
# so a promiscuous domain doesn't vote many times)
rec_sk, rec_t = {}, {}
for _, r in pairs:
    if r in rec_sk:
        continue
    s = acc_sk.get(acc_of(r))
    if s:
        rec_sk[r] = s; rec_t[r] = tg.get(r)
glob_c = Counter(rec_sk.values())
tot = sum(glob_c.values())
print("global lineage mix over distinct lit records:")
for k, n in glob_c.most_common():
    print(f"  {k:12s} {n:>8,}  {100*n/tot:5.1f}%")
P_EUK = glob_c.get("Eukaryota", 0) / tot

fold_euk = defaultdict(lambda: [0, 0])          # fold -> [euk, total]
for r, s in rec_sk.items():
    t = rec_t.get(r)
    if t:
        fold_euk[t][1] += 1
        if s == "Eukaryota":
            fold_euk[t][0] += 1

cell_recs = defaultdict(set)
for cid, r in pairs:
    if r in rec_sk:
        cell_recs[cid].add(r)

A = np.load(f"{G}/grid_afdb_nT_rebuilt.npy")
IMP = np.load(f"{G}/impossible_mask.npy")

rows = []
for cid, recs in cell_recs.items():
    sk, ty = cid // 32, cid % 32
    if IMP[sk, ty] or len(recs) < 30:
        continue
    n = len(recs)
    e = sum(1 for r in recs if rec_sk[r] == "Eukaryota")
    # fold-predicted euk count for this cell
    exp = 0.0
    for r in recs:
        t = rec_t.get(r)
        if t and fold_euk[t][1] >= 5:
            exp += fold_euk[t][0] / fold_euk[t][1]
        else:
            exp += P_EUK
    p_raw = binomtest(e, n, P_EUK).pvalue
    p_fold = binomtest(e, n, min(max(exp / n, 1e-9), 1 - 1e-9)).pvalue
    rows.append((sk, ty, n, e, e / n, exp / n, p_raw, p_fold))

q_raw = multipletests([r[6] for r in rows], method="fdr_bh")[1]
q_fold = multipletests([r[7] for r in rows], method="fdr_bh")[1]
print(f"\ncells tested (>=30 lineage-typed records, satisfiable): {len(rows):,}")
print(f"  euk-skewed vs GLOBAL   q<0.05: {(q_raw < 0.05).sum():,}")
print(f"  euk-skewed vs OWN FOLDS q<0.05: {(q_fold < 0.05).sum():,}   <- not explained by fold content")


def strands(ty): return bin(ty).count("1")


print("\nTop 15 eukaryote-enriched cells vs global (raw):")
print(f"{'sk':>4}{'ty':>3} {'E':>2} {'n':>6} {'euk%':>6} {'foldpred%':>9} {'q_raw':>9} {'q_fold':>9}")
for i in np.argsort([r[4] for r in rows])[::-1][:15]:
    r = rows[i]
    print(f"{r[0]:>4}{r[1]:>3} {strands(r[1]):>2} {r[2]:>6,} {100*r[4]:>5.1f}% {100*r[5]:>8.1f}% "
          f"{q_raw[i]:>9.1e} {q_fold[i]:>9.1e}")

sig_fold = [i for i in range(len(rows)) if q_fold[i] < 0.05]
print(f"\nCells still skewed AFTER controlling for fold content: {len(sig_fold)}")
print(f"{'sk':>4}{'ty':>3} {'E':>2} {'n':>6} {'euk%':>6} {'foldpred%':>9} {'q_fold':>9}")
for i in sorted(sig_fold, key=lambda i: q_fold[i])[:15]:
    r = rows[i]
    print(f"{r[0]:>4}{r[1]:>3} {strands(r[1]):>2} {r[2]:>6,} {100*r[4]:>5.1f}% {100*r[5]:>8.1f}% {q_fold[i]:>9.1e}")

# is the raw signal just the all-helix cells?
raw_sig = [i for i in range(len(rows)) if q_raw[i] < 0.05 and rows[i][4] > P_EUK]
by_str = Counter(strands(rows[i][1]) for i in raw_sig)
allc = Counter(strands(r[1]) for r in rows)
print("\nEuk-enriched cells by strand count (is it only the all-helix ones?):")
for s in sorted(allc):
    k = by_str.get(s, 0)
    print(f"  {s} strands: {k:>4}/{allc[s]:<5} cells enriched ({100*k/allc[s]:4.0f}%)")
json.dump({"rows": [[int(a) for a in r[:4]] + [float(x) for x in r[4:]] for r in rows],
           "q_raw": q_raw.tolist(), "q_fold": q_fold.tolist()},
          open(f"{G}/tax_cells.json", "w"))
print("\nwrote tax_cells.json")
