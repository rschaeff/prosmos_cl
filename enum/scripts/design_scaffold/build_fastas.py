#!/usr/bin/env python3
"""Build matched FASTAs for the 50%/90% clustering of both databases' hitters.

AFDB: the distinct grey non-singleton reps that lit >=1 cell (distinct_hitters),
      sequences from grey's cluster-rep parquet.
PDB : the distinct ECOD domains that lit >=1 cell, sequences via
      ecod_uid -> ecod_commons.domains.id -> domain_sequences.

Also dumps per-cell PDB uid sets (pickle) so cells can be recounted in
50%-cluster units after mmseqs.
"""
import subprocess, re, glob, io, pickle
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
import psycopg2
import pyarrow.parquet as pq

R = Path("/home/rschaeff/work/prosmos_2026")
D = R / "design_scaffold"
ZSTD = "/sw/apps/Anaconda3-2023.09-0/bin/zstd"
QP = re.compile(r'^s5-(\d{4})-(\d{4})$')

# ---------- PDB: per-cell uid sets + distinct uids ----------
def scan(fn):
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

parts = sorted(glob.glob(str(R / "s5_pdb_inv" / "hitparts" / "*.tsv.zst")))
pdb_cell_uids = defaultdict(set)
with Pool(16) as pool:
    for d in pool.imap_unordered(scan, parts, chunksize=4):
        for c, us in d.items():
            pdb_cell_uids[c] |= us
pickle.dump(dict(pdb_cell_uids), open(D / "pdb_cell_uids.pkl", "wb"))
all_uids = sorted(set().union(*pdb_cell_uids.values()))
print(f"PDB: {len(pdb_cell_uids):,} cells, {len(all_uids):,} distinct lit domains", flush=True)

# sequences for those uids
c = psycopg2.connect(host='dione', port=45000, user='ecod', dbname='ecod_protein')
cur = c.cursor()
cur.execute("create temp table qu(ecod_uid int primary key)")
cur.copy_from(io.StringIO("\n".join(map(str, all_uids)) + "\n"), 'qu', columns=('ecod_uid',))
cur.execute("""select d.ecod_uid, ds.sequence
               from ecod_commons.domains d
               join qu on qu.ecod_uid = d.ecod_uid
               join ecod_commons.domain_sequences ds on ds.domain_id = d.id
               where not d.is_obsolete""")
n = 0
with open(D / "pdb_hitters.fasta", "w") as fh:
    for uid, seq in cur:
        if seq:
            fh.write(f">{uid}\n{seq}\n"); n += 1
c.close()
print(f"PDB: wrote {n:,} sequences to pdb_hitters.fasta", flush=True)

# ---------- AFDB: distinct hitting reps + sequences from grey's parquet ----------
want = set()
for line in open(R / "s5_inv" / "distinct_hitters.txt"):
    a = line.strip()
    if a.startswith("pdb"):
        a = a[3:]
    if a:
        want.add(a)
print(f"AFDB: {len(want):,} distinct hitting reps", flush=True)

pf = pq.ParquetFile("/home/grey/afdb.200m/summary/non_singleton_4p9M_cluster_reps_seqs.parquet")
n = 0
with open(D / "afdb_hitters.fasta", "w") as fh:
    for batch in pf.iter_batches(columns=["cluster_name", "sequence"], batch_size=200000):
        names = batch.column("cluster_name").to_pylist()
        seqs = batch.column("sequence").to_pylist()
        for name, seq in zip(names, seqs):
            if name in want and seq:
                fh.write(f">{name}\n{seq}\n"); n += 1
print(f"AFDB: wrote {n:,} sequences to afdb_hitters.fasta  "
      f"(missing {len(want)-n:,} reps not found in parquet)", flush=True)
