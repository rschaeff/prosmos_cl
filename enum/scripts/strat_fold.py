#!/usr/bin/env python3
"""Are any S5 topologies serviced ONLY by a particular structural context?
Two cuts over the experimental PDB: experimental METHOD (is anything EM-only?)
and SOURCE SUPERKINGDOM (is anything viral-only?).

Counted at FOLD level (distinct ECOD T-groups), not domain level. At domain level
the PDB's redundancy dominates: e.g. sk178/ty16 has 168 viral domains and looks
overwhelming (p~1e-150 against a 7.7% baseline), but they are ONE T-group
(4967.1.1, hepatitis C) solved 168 times for drug design. Fold level is the only
unit in which "this cell is serviced only by X" means anything.

Taxonomy: ecod_commons.domain_taxonomy on dione (PDB-native, keyed by ecod_uid,
3.1M rows, INCLUDES Viruses). NOT afdb_200m.taxonomy -- that is AlphaFold-derived
and contains no viruses at all, so routing through it returns 0 viral by
construction.
"""
import psycopg2, glob, subprocess, re, io, json
import numpy as np
from collections import defaultdict, Counter
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

G = "/home/rschaeff/work/prosmos_2026/s5_grid"
P = "/home/rschaeff/work/prosmos_2026/s5_pdb_inv"

# uid -> T-group, uid -> pdb
uid_t, uid_pdb = {}, {}
for line in open("/home/rschaeff/work/prosmos_2026/pdb_exp_build/uid_class.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 5:
        if f[2]:
            uid_t[int(f[0])] = f[2]
        if f[4]:
            uid_pdb[int(f[0])] = f[4].lower()

QP = re.compile(r'^s5-(\d{4})-(\d{4})$')
cell_uids = defaultdict(set)
for f in sorted(glob.glob(P + "/hitparts/*.tsv.zst")):
    p = subprocess.Popen(["/sw/apps/Anaconda3-2023.09-0/bin/zstd", "-dc", f],
                         stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        t = line.rstrip("\n").split("\t")
        if len(t) < 2:
            continue
        m = QP.match(t[0])
        if m:
            cell_uids[(int(m.group(1)), int(m.group(2)))].add(int(t[1]))
    p.wait()
uids = sorted({u for cs in cell_uids.values() for u in cs})
print(f"cells {len(cell_uids):,}  distinct lit PDB domains {len(uids):,}", flush=True)

# method
c = psycopg2.connect(host='sangala', port=45000, user='ecod', dbname='ecod_af2_pdb')
cur = c.cursor(); cur.execute("select lower(pdb), method from pdb_info")
pdb_m = {a: b for a, b in cur.fetchall()}; c.close()
# superkingdom (PDB-native)
c = psycopg2.connect(host='dione', port=45000, user='ecod', dbname='ecod_protein')
cur = c.cursor()
cur.execute("create temp table qu(uid bigint primary key)")
cur.copy_from(io.StringIO("\n".join(map(str, uids)) + "\n"), 'qu', columns=('uid',))
cur.execute("select d.ecod_uid, d.superkingdom from ecod_commons.domain_taxonomy d join qu on qu.uid=d.ecod_uid")
uid_sk = {u: s for u, s in cur.fetchall()}; c.close()
print(f"  method-typed {sum(1 for u in uids if pdb_m.get(uid_pdb.get(u))):,} | "
      f"superkingdom-typed {sum(1 for u in uids if uid_sk.get(u)):,}", flush=True)


def fold_label(getter, name):
    """For each cell: the set of T-GROUPS, and how they split by the attribute.
    A T-group counts for a category if ANY of its domains in that cell has it."""
    rows, only = [], []
    # global fold-level base rate: fraction of T-groups that ever carry the category
    fold_cat, fold_all = defaultdict(set), set()
    for cs in cell_uids.values():
        for u in cs:
            t = uid_t.get(u)
            if not t:
                continue
            fold_all.add(t)
            v = getter(u)
            if v:
                fold_cat[v].add(t)
    return fold_cat, fold_all


def analyse(getter, cat, label):
    # per cell: distinct T-groups, and distinct T-groups carrying `cat`
    rows, only = [], []
    tot_f, cat_f = set(), set()
    for cs in cell_uids.values():
        for u in cs:
            t = uid_t.get(u)
            if not t:
                continue
            tot_f.add(t)
            if getter(u) == cat:
                cat_f.add(t)
    base = len(cat_f) / len(tot_f)
    for (a, b), us in cell_uids.items():
        tg_all, tg_cat = set(), set()
        for u in us:
            t = uid_t.get(u)
            if not t:
                continue
            tg_all.add(t)
            if getter(u) == cat:
                tg_cat.add(t)
        n = len(tg_all)
        if n == 0:
            continue
        k = len(tg_cat)
        if k == n:
            only.append((a, b, n))
        if n >= 5:
            rows.append((a, b, n, k, k / n, binomtest(k, n, base, alternative='greater').pvalue))
    q = multipletests([r[5] for r in rows], method='fdr_bh')[1] if rows else np.array([])
    print(f"\n=== {label}  (FOLD level)")
    print(f"  global base rate: {100*base:.1f}% of T-groups")
    print(f"  cells tested (>=5 T-groups): {len(rows):,}   enriched q<0.05: {int((q<0.05).sum()):,}")
    multi = [o for o in only if o[2] >= 3]
    print(f"  {label}-ONLY cells: {len(only)}  |  with >=3 distinct T-groups: {len(multi)}")
    for a, b, n in sorted(only, key=lambda x: -x[2])[:8]:
        print(f"      sk{a} ty{b}: {n} T-group(s) — {'MULTI-FOLD' if n>=3 else 'single/near-single fold'}")
    if len(rows):
        print(f"  top enriched: ", end="")
        for i in np.argsort(q)[:5]:
            r = rows[i]
            print(f"sk{r[0]}/ty{r[1]} {r[3]}/{r[2]}f ", end="")
        print()
    return dict(base=base, tested=len(rows), enriched=int((q < 0.05).sum()),
                only=len(only), only_multifold=len(multi))


res = {}
res["EM"] = analyse(lambda u: pdb_m.get(uid_pdb.get(u)), "ELECTRON MICROSCOPY", "EM")
res["viral"] = analyse(lambda u: uid_sk.get(u), "Viruses", "viral")
json.dump(res, open(f"{G}/strat_fold.json", "w"), indent=1)
print("\nwrote strat_fold.json")
