#!/usr/bin/env python3
"""Rewind the PDB and watch the "AFDB-only" cells fill in.

Every cross-database delta we compute is confounded by sampling depth, and the
depth argument is unfalsifiable on its own. This makes it falsifiable: rebuild the
PDB side using only structures RELEASED BEFORE year T, recount the cells that are
AFDB-lit but PDB-dark, and ask what happened to them as the PDB grew.

If "AFDB-only" means "real but not yet crystallised", those cells should fill in
at a measurable rate as T advances -- and the ones that resist decades of growth
are the interesting residue. If instead the count is flat, the cells are telling
us about the grid's vocabulary rather than about experimental coverage.

Fold level throughout (a cell is PDB-lit at vintage T if >=1 distinct ECOD T-group
has a domain in a structure released <= T). Domain level would let one
heavily-redeposited drug target masquerade as broad coverage.
"""
import psycopg2, glob, subprocess, re, io, json
import numpy as np
from collections import defaultdict

G = "/home/rschaeff/work/prosmos_2026/s5_grid"
P = "/home/rschaeff/work/prosmos_2026/s5_pdb_inv"
VINTAGES = [1995, 2000, 2005, 2010, 2013, 2016, 2019, 2022, 2026]

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
print(f"cells with PDB hits: {len(cell_uids):,}", flush=True)

pdbs = sorted({v for v in uid_pdb.values()})
c = psycopg2.connect(host='dione', port=45000, user='ecod', dbname='ecod_protein')
cur = c.cursor()
cur.execute("create temp table qp(pdb text primary key)")
cur.copy_from(io.StringIO("\n".join(pdbs) + "\n"), 'qp', columns=('pdb',))
cur.execute("""select lower(e.pdb_id), extract(year from coalesce(e.release_date, e.deposition_date))
               from pdb_analysis.pdb_entries e join qp on qp.pdb = lower(e.pdb_id)""")
pdb_year = {a: int(b) for a, b in cur.fetchall() if b}
c.close()
print(f"pdb -> year: {len(pdb_year):,}", flush=True)

A = np.load(f"{G}/grid_afdb_nT_rebuilt.npy")
IMP = np.load(f"{G}/impossible_mask.npy")
SAT = ~IMP
afdb_lit = (A > 0) & SAT

# per cell: fold -> earliest year that fold appears in this cell
cell_fold_year = defaultdict(dict)
for (sk, ty), us in cell_uids.items():
    for u in us:
        t, pd = uid_t.get(u), uid_pdb.get(u)
        if not t or not pd:
            continue
        y = pdb_year.get(pd)
        if y is None:
            continue
        d = cell_fold_year[(sk, ty)]
        if t not in d or y < d[t]:
            d[t] = y

series, masks = [], {}
for T in VINTAGES:
    pdb_lit = np.zeros_like(afdb_lit)
    for (sk, ty), fy in cell_fold_year.items():
        if any(y <= T for y in fy.values()):
            pdb_lit[sk, ty] = True
    pdb_lit &= SAT
    only = afdb_lit & ~pdb_lit
    masks[T] = only.copy()
    series.append((T, int(pdb_lit.sum()), int(only.sum()), int((pdb_lit & ~afdb_lit).sum())))
    print(f"  {T}: PDB-lit {pdb_lit.sum():>5,}   AFDB-only {only.sum():>5,}   PDB-only {int((pdb_lit&~afdb_lit).sum()):>4,}", flush=True)

print("\nFILL RATE — of the cells AFDB-only at vintage T, how many are PDB-lit by 2026?")
final = masks[VINTAGES[-1]]
for T in VINTAGES[:-1]:
    m = masks[T]
    n = int(m.sum())
    filled = int((m & ~final).sum())
    print(f"  AFDB-only in {T}: {n:>5,}  ->  filled by 2026: {filled:>5,} ({100*filled/max(n,1):4.1f}%)  "
          f"still dark: {n-filled:,}")

# how long has each still-dark cell resisted? (first vintage at which it was AFDB-only)
first_seen = {}
for T in VINTAGES:
    for sk, ty in zip(*np.where(masks[T])):
        first_seen.setdefault((int(sk), int(ty)), T)
resist = [(c, y) for c, y in first_seen.items() if final[c]]
print(f"\nstill AFDB-only in 2026: {len(resist):,}")
from collections import Counter
cnt = Counter(y for _, y in resist)
for y in VINTAGES:
    if cnt.get(y):
        print(f"  AFDB-only continuously since {y}: {cnt[y]:>4,}")
json.dump({"series": series,
           "fill": [[T, int(masks[T].sum()), int((masks[T] & ~final).sum())] for T in VINTAGES[:-1]]},
          open(f"{G}/pdb_vintage.json", "w"), indent=1)
print("\nwrote pdb_vintage.json")
