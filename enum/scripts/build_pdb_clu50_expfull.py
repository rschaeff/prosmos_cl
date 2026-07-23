#!/usr/bin/env python3
"""Build a MATCHED 50%-id sequence-cluster map for the new 1.4M-sweep PDB domains.

The cluster-unit significance correction requires PDB clustered at the SAME
threshold as AFDB (afdb_clu50 = 50%-id / 90%-cov mmseqs). The old pdb_clu50 map
covers only the 496k old-sweep domains; using ECOD F-group instead is too coarse
(10k families vs 1.18M AFDB clusters -> 100x scale asymmetry that INFLATES rather
than deflates, exactly what the correction is meant to remove).

This pulls sequences for the distinct domains that lit >=1 S5 cell in
s5_exp_full_g198, clusters them at 50%/90% with mmseqs (matching afdb_clu50), and
writes pdb_clu50_expfull_cluster.tsv (rep<TAB>member, ecod_uid keys) -- the same
format cluster_significance_expfull.py consumes.

Sequence pull mirrors design_scaffold/build_fastas.py.
"""
import io
import pickle
import subprocess
from pathlib import Path

import psycopg2

R = Path("/home/rschaeff/work/prosmos_2026")
G = R / "s5_grid"
OUT = R / "design_scaffold"
MMSEQS = "/sw/apps/mmseqs/bin/mmseqs"
FASTA = OUT / "pdb_hitters_expfull.fasta"
TMP = R / "mmseqs_tmp_expfull"
PREFIX = OUT / "pdb_clu50_expfull"

lit = pickle.load(open(G / "exp_full_uids.pkl", "rb"))
print(f"distinct lit PDB domains (s5_exp_full_g198): {len(lit):,}", flush=True)

# ---- sequences via ecod_commons.domain_sequences (build_fastas.py idiom)
c = psycopg2.connect(host="dione", port=45000, user="ecod", dbname="ecod_protein")
cur = c.cursor()
cur.execute("create temp table qu(ecod_uid int primary key)")
cur.copy_from(io.StringIO("\n".join(map(str, sorted(lit))) + "\n"), "qu",
              columns=("ecod_uid",))
cur.execute("""select d.ecod_uid, ds.sequence
               from ecod_commons.domains d
               join qu on qu.ecod_uid = d.ecod_uid
               join ecod_commons.domain_sequences ds on ds.domain_id = d.id
               where not d.is_obsolete""")
n = 0
with open(FASTA, "w") as fh:
    for uid, seq in cur:
        if seq:
            fh.write(f">{uid}\n{seq}\n")
            n += 1
c.close()
print(f"wrote {n:,} sequences -> {FASTA}", flush=True)

# ---- mmseqs easy-cluster at 50% id / 90% cov (matched to afdb_clu50)
TMP.mkdir(exist_ok=True)
subprocess.run([MMSEQS, "easy-cluster", str(FASTA), str(PREFIX), str(TMP),
                "--min-seq-id", "0.5", "-c", "0.9", "--cov-mode", "1",
                "-v", "1"], check=True)
# easy-cluster writes <prefix>_cluster.tsv as rep<TAB>member
clu = PREFIX.with_name(PREFIX.name + "_cluster.tsv")
final = OUT / "pdb_clu50_expfull_cluster.tsv"
clu.rename(final)
nmem = sum(1 for _ in open(final))
nrep = len({line.split("\t")[0] for line in open(final)})
print(f"clustered: {nmem:,} members -> {nrep:,} clusters "
      f"({100*nrep/nmem:.1f}% retained) -> {final}", flush=True)
