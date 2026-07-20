#!/usr/bin/env python3
"""Populate a local store of AFDB DPAM domain structures for afdb_assigned exemplars.

Source (grey's, world-readable): /home/grey/afdb.200m/non_singleton_4p9m_structures
  new/<dir>/<tar>.tar.gz   -> '<acc>_D<n>.pdb'      (98.6% of our sweep)
  dpam/<dir>/<tar>.tar.gz  -> 'dpam_<acc>_nD<n>.pdb'(1.2%)
  ecod/ecod_<uid>.pdb.gz   -> (0.2%)
Index: /home/grey/resources/afdb200M/356M_pdb_folder_maps/<dir>.list
       columns: dir, tar, <n>, accession    and tar == accession[-4:-2]

These are ALREADY chopped to the DPAM domain (verified: A0A1J4NAJ4_D1 spans 4-115,
matching afdb_200m.ecod_domain_range), so no server-side range chop is needed —
we serve the domain directly, keyed by did.

We COPY what we need into our own store rather than reading grey's paths at
runtime, so the app cannot break if that tree is reorganised.
Writes: ~/work/prosmos_2026/afdb_domain_struct/<did>.pdb
"""
import json, glob, os, sys, tarfile, re
from pathlib import Path
from collections import defaultdict

SRC = Path("/home/grey/afdb.200m/non_singleton_4p9m_structures")
MAPS = Path("/home/grey/resources/afdb200M/356M_pdb_folder_maps")
CELLS = Path.home() / "dev/prosmos_inspect/data/afdb_assigned/cell"
OUT = Path.home() / "work/prosmos_2026/afdb_domain_struct"
OUT.mkdir(parents=True, exist_ok=True)

# 1) which domains do the exemplars need?
need = {}          # did -> accession
for cf in sorted(CELLS.glob("*.json")):
    d = json.loads(cf.read_text())
    for e in d["exemplars"]:
        did = e["did"]
        acc = (e.get("unp") or did.split("_D")[0]).replace("dpam_", "")
        need[did] = acc
have = {p.stem for p in OUT.glob("*.pdb")}
todo = {k: v for k, v in need.items() if k not in have}
print(f"exemplar domains needed {len(need):,}  already have {len(have):,}  to fetch {len(todo):,}", flush=True)
if not todo:
    sys.exit(0)

# 2) accession -> shard dir, from the folder maps (tar name is computable)
want_acc = set(todo.values())
acc_dir = {}
for lf in sorted(MAPS.glob("*.list")):
    with open(lf, errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[3] in want_acc:
                acc_dir[p[3]] = p[0]
    if len(acc_dir) == len(want_acc):
        break
print(f"accessions located in the folder maps: {len(acc_dir):,}/{len(want_acc):,}", flush=True)

# 3) group the work by archive so each tar is opened exactly once
by_tar = defaultdict(list)      # (family, dir, tar) -> [(did, member_name)]
for did, acc in todo.items():
    d = acc_dir.get(acc)
    if not d:
        continue
    tar = acc[-4:-2]
    if did.startswith("dpam_"):
        by_tar[("dpam", d, tar)].append((did, f"{did}.pdb"))
    else:
        by_tar[("new", d, tar)].append((did, f"{did}.pdb"))
print(f"archives to open: {len(by_tar):,}", flush=True)

n = miss = 0
for i, ((fam, d, tar), members) in enumerate(sorted(by_tar.items()), 1):
    path = SRC / fam / d / f"{tar}.tar.gz"
    if not path.exists():
        miss += len(members); continue
    wanted = {m: did for did, m in members}
    try:
        with tarfile.open(path, "r:gz") as tf:
            for ti in tf:
                did = wanted.get(os.path.basename(ti.name))
                if not did:
                    continue
                f = tf.extractfile(ti)
                if f:
                    (OUT / f"{did}.pdb").write_bytes(f.read())
                    n += 1
    except (tarfile.TarError, OSError) as ex:
        print(f"  archive fail {path}: {ex}", file=sys.stderr)
        miss += len(members)
    if i % 500 == 0:
        print(f"  {i}/{len(by_tar)} archives, {n:,} written", flush=True)
print(f"wrote {n:,} domain structures ({miss:,} unavailable) -> {OUT}")
