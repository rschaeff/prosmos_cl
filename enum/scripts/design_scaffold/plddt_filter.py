#!/usr/bin/env python3
"""pLDDT-gate the under-templated cells' AFDB scaffolds.

A design audience will not template off low-confidence models, so the 53
under-templated cells only count if their AFDB scaffolds are actually
well-predicted. Recount each cell keeping only scaffolds at pLDDT >= 70 and >= 80
and ask whether the PDB<=5 / AFDB>=50 asymmetry survives.

pLDDT source: afdb_200m.protein_plddt (per-protein mean). This is a proxy for the
per-domain value; it is exact for single-domain proteins and an average otherwise,
so we also record how many of these scaffolds are single-domain (where the proxy
is exact) and validate a sample against the DPAM-chopped B-factor separately.
"""
import json, io
from pathlib import Path
from collections import defaultdict
import psycopg2

D = Path("/home/rschaeff/work/prosmos_2026/design_scaffold")

scaffolds = [s.strip() for s in open(D / "under_scaffolds.txt") if s.strip()]
# scaffold names come in three shapes: ACC_D<n>, ACC_nD<n>, and dpam_ACC[...].
# Per-protein pLDDT keys on the accession only, so extract that robustly.
def parse(s):
    s2 = s[5:] if s.startswith("dpam_") else s
    for sep in ("_nD", "_D"):
        if sep in s2:
            acc, dn = s2.split(sep, 1)
            try:
                return acc, int(dn)
            except ValueError:
                return acc, None
    return s2, None
acc_of = {s: parse(s) for s in scaffolds}
accs = sorted({a for a, _ in acc_of.values()})
print(f"scaffolds {len(scaffolds):,}  distinct accessions {len(accs):,}")

c = psycopg2.connect(host='lotta', port=45000, user='ecod', dbname='ecod_protein')
cur = c.cursor()
cur.execute("create temp table qa(acc text primary key)")
cur.copy_from(io.StringIO("\n".join(accs) + "\n"), 'qa', columns=('acc',))
# accession -> protein id, per-protein mean pLDDT
cur.execute("""select p.uniprot_acc, p.id, pp.avg_plddt
               from afdb_200m.protein p
               join qa on qa.acc = p.uniprot_acc
               join afdb_200m.protein_plddt pp on pp.protein_id = p.id""")
acc_pl, acc_pid = {}, {}
for acc, pid, pl in cur.fetchall():
    acc_pl[acc] = pl
    acc_pid[acc] = pid
print(f"accessions with per-protein pLDDT: {len(acc_pl):,}/{len(accs):,}")

# how many DPAM domains does each protein have? (single-domain -> proxy is exact)
pids = sorted(set(acc_pid.values()))
cur.execute("create temp table qp(pid bigint primary key)")
cur.copy_from(io.StringIO("\n".join(map(str, pids)) + "\n"), 'qp', columns=('pid',))
cur.execute("""select r.protein_id, count(*) from afdb_200m.ecod_domain_range r
               join qp on qp.pid = r.protein_id group by r.protein_id""")
ndom = {pid: n for pid, n in cur.fetchall()}
c.close()

single = sum(1 for a in acc_pl if ndom.get(acc_pid[a], 1) == 1)
print(f"single-domain proteins (proxy exact): {single:,}/{len(acc_pl):,} "
      f"({100*single/max(len(acc_pl),1):.0f}%)")

# scaffold -> pLDDT
sc_pl = {s: acc_pl.get(acc_of[s][0]) for s in scaffolds}
have = [v for v in sc_pl.values() if v is not None]
import numpy as np
have = np.array(have)
print(f"\nscaffold pLDDT: n={len(have):,}  median {np.median(have):.1f}  "
      f">=70 {100*(have>=70).mean():.0f}%  >=80 {100*(have>=80).mean():.0f}%  "
      f">=90 {100*(have>=90).mean():.0f}%")

# ---- recount the 53 cells at pLDDT gates
m = json.load(open(D / "matched.json"))
under = {(u["sk"], u["ty"]): u for u in m["under_templated_matched"]}
cell_sc = defaultdict(set)
for line in open(D / "afdb_cell_reps.tsv"):
    sk, ty, acc = line.rstrip("\n").split("\t")
    k = (int(sk), int(ty))
    if k in under:
        cell_sc[k].add(acc)

print(f"\n{'cell':>12} {'PDB':>4} {'AFDBraw':>7} {'AFDB>=70':>8} {'AFDB>=80':>8} {'mult80':>6}")
surv70, surv80 = 0, 0
rows = []
for k, u in sorted(under.items(), key=lambda kv: -kv[1]["afdb"]):
    scs = cell_sc[k]
    n70 = sum(1 for s in scs if (sc_pl.get(s) or 0) >= 70)
    n80 = sum(1 for s in scs if (sc_pl.get(s) or 0) >= 80)
    pdb = u["pdb"]
    if pdb <= 5 and n70 >= 50:
        surv70 += 1
    if pdb <= 5 and n80 >= 50:
        surv80 += 1
    rows.append((k, pdb, u["afdb"], n70, n80))
for k, pdb, raw, n70, n80 in rows[:20]:
    print(f"  sk{k[0]:03d} ty{k[1]:02d} {pdb:>4} {raw:>7} {n70:>8} {n80:>8} {n80/max(pdb,1):>5.0f}x")

print(f"\nunder-templated cells surviving pLDDT>=70 (AFDB>=50): {surv70}/53")
print(f"under-templated cells surviving pLDDT>=80 (AFDB>=50): {surv80}/53")

json.dump({
    "n_scaffolds": len(scaffolds),
    "pl_coverage": len(acc_pl) / len(accs),
    "single_domain_frac": single / max(len(acc_pl), 1),
    "scaffold_plddt_median": float(np.median(have)),
    "frac_ge70": float((have >= 70).mean()), "frac_ge80": float((have >= 80).mean()),
    "survive70": surv70, "survive80": surv80,
    "cells": [{"sk": k[0], "ty": k[1], "pdb": pdb, "afdb_raw": raw,
               "afdb_ge70": n70, "afdb_ge80": n80} for k, pdb, raw, n70, n80 in rows],
}, open(D / "plddt_filter.json", "w"), indent=1)
print(f"\nwrote {D/'plddt_filter.json'}")
