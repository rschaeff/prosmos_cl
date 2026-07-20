"""Size-matched rarefaction: are the 619 AFDB-only cells a depth artifact?

AFDB searched 4.92M records; PDB 496K -- 10x depth. Subsample AFDB's searched
records down toward PDB depth, recompute distinct-fold-per-cell, and count cells
that are AFDB-lit but PDB-dark (the AFDB-only set) against the FULL 100% PDB grid.
If AFDB-only collapses toward the PDB-only count (168) at matched depth, the excess
is depth; if it holds, it is real AFDB-specific topology sampling.

Fold assignment replicates nt_shard.py exactly (tg keyed on full afdb_name_tgroup
names, looked up by the 15-char record key) so the full-population rebuild must
reproduce grid_afdb_nT_rebuilt.npy -- a sanity gate before trusting the rarefaction.
"""
import numpy as np, glob, pickle, random
from collections import defaultdict

G = "/home/rschaeff/work/prosmos_2026/s5_grid"

# fold map, exactly as nt_shard.py built it
tg = {}
for line in open(f"{G}/afdb_name_tgroup.tsv"):
    a, b = line.rstrip("\n").split("\t")
    p = b.split(".")
    tg[a] = ".".join(p[:2]) if len(p) >= 2 else b   # H-group... but grid used T? check below

# grid_afdb_nT counts distinct T-GROUPS. nt_shard used the FULL t_group id.
tgT = {}
for line in open(f"{G}/afdb_name_tgroup.tsv"):
    a, b = line.rstrip("\n").split("\t")
    tgT[a] = b

# load all (cell_id, record15) pairs
pairs = []
for f in glob.glob(f"{G}/afdb_rare/*.pkl"):
    pairs.extend(pickle.load(open(f, "rb")))
print(f"loaded {len(pairs):,} (cell,record) pairs")

# record -> its fold (T-group); None if unannotated
recs = sorted({r for _, r in pairs})
recfold = {r: tgT.get(r) for r in recs}
print(f"distinct lit records: {len(recs):,}  annotated: {sum(1 for r in recs if recfold[r]):,}")

# cell -> set of records (only need annotated ones for fold counting)
cell_recs = defaultdict(set)
for cid, r in pairs:
    if recfold[r]:
        cell_recs[cid].add(r)


def grid_from(keep):
    """distinct folds per cell over the kept record set -> 198x32 int grid"""
    g = np.zeros((198, 32), int)
    for cid, rs in cell_recs.items():
        folds = {recfold[r] for r in rs if r in keep}
        if folds:
            g[cid // 32, cid % 32] = len(folds)
    return g


# --- sanity gate: full rebuild must match grid_afdb_nT_rebuilt
allrec = set(recs)
full = grid_from(allrec)
ref = np.load(f"{G}/grid_afdb_nT_rebuilt.npy")
# rebuilt array counts folds; compare where both nonzero pattern
match = (full > 0) == (ref > 0)
print(f"\nSANITY: full rebuild vs grid_afdb_nT_rebuilt -- lit-pattern agreement "
      f"{100*match.mean():.1f}%  (cells differing: {(~match).sum()})")
corr = np.corrcoef(full.ravel(), ref.ravel())[0, 1]
print(f"        fold-count correlation: {corr:.4f}")

# --- rarefaction
P = np.load(f"{G}/grid_pdb_nT.npy").astype(float)
imp = np.load(f"{G}/impossible_mask.npy")
sat = ~imp
pl = (P > 0) & sat
NSEARCH = 4_921_931
PDB_SEARCH = 496_359
# rarefy the SEARCHED population; lit records survive iff their (searched) draw survives.
# we only hold lit records, so approximate by subsampling lit records at the same rate
# (uniform subsampling of searched == uniform subsampling of lit, since lit is a subset).
rng = random.Random(0)
reclist = list(allrec)
print(f"\n{'AFDB depth':>11} {'kept recs':>10} {'AFDB-lit':>9} {'AFDB-only':>10} {'PDB-only(ref 168)':>17}")
for frac in [1.0, 0.75, 0.5, 0.25, PDB_SEARCH / NSEARCH]:
    aos = []
    for rep in range(3 if frac < 1.0 else 1):
        k = max(1, int(len(reclist) * frac))
        keep = set(rng.sample(reclist, k))
        g = grid_from(keep)
        al = (g > 0) & sat
        ao = (al & ~pl).sum()
        po = (pl & ~al).sum()
        aos.append((al.sum(), ao, po))
    al_m = np.mean([x[0] for x in aos]); ao_m = np.mean([x[1] for x in aos]); po_m = np.mean([x[2] for x in aos])
    print(f"{frac:>10.0%} {int(len(reclist)*frac):>10,} {al_m:>9.0f} {ao_m:>10.0f} {po_m:>17.0f}")

print("\nread: if AFDB-only at ~10% depth (PDB-matched) approaches PDB-only(168),")
print("the 619 excess is depth; if it stays well above, it is real.")
