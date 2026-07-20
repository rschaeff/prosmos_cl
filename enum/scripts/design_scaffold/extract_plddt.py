#!/usr/bin/env python3
"""Gold-standard per-domain pLDDT for the surviving-cell scaffolds.

The grey DPAM structures are ALREADY chopped to the domain, and the AF model's
per-residue pLDDT lives in the B-factor column, so the mean CA B-factor of the
chopped domain IS the exact per-domain pLDDT -- no boundary approximation (unlike
ted_domain.plddt) and no whole-protein averaging (unlike protein_plddt, which was
exact for only 37% of these).

Locate each scaffold via grey's folder maps (dir, tar, n, accession; tar =
accession[-4:-2]), open each archive once, read the member's CA B-factors in
memory, discard the bytes. Writes did -> mean pLDDT.
"""
import os, sys, tarfile, io
from pathlib import Path
from collections import defaultdict

SRC = Path("/home/grey/afdb.200m/non_singleton_4p9m_structures")
MAPS = Path("/home/grey/resources/afdb200M/356M_pdb_folder_maps")
D = Path("/home/rschaeff/work/prosmos_2026/design_scaffold")
OUT = D / "surv_plddt.tsv"

dids = [s.strip() for s in open(D / "surv_scaffolds.txt") if s.strip()]

def acc_of(did):
    s = did[5:] if did.startswith("dpam_") else did
    for sep in ("_nD", "_D"):
        if sep in s:
            return s.split(sep, 1)[0]
    return s

want_acc = {acc_of(d): d for d in dids}      # accession -> did (1:1 here)
print(f"scaffolds {len(dids):,}  distinct accessions {len(want_acc):,}", flush=True)

# accession -> shard dir, from the folder maps (break early once all found)
acc_dir = {}
need = set(want_acc)
for lf in sorted(MAPS.glob("*.list")):
    with open(lf, errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[3] in need:
                acc_dir[p[3]] = p[0]
                need.discard(p[3])
    if not need:
        break
print(f"located {len(acc_dir):,}/{len(want_acc):,} accessions  (unlocated {len(need):,})", flush=True)

# group by archive so each tar opens once
by_tar = defaultdict(list)   # (fam, dir, tar) -> [(did, member)]
for acc, did in want_acc.items():
    d = acc_dir.get(acc)
    if not d:
        continue
    tar = acc[-4:-2]
    fam = "dpam" if did.startswith("dpam_") else "new"
    by_tar[(fam, d, tar)].append((did, f"{did}.pdb"))
print(f"archives to open: {len(by_tar):,}", flush=True)


def mean_ca_plddt(raw):
    vals = []
    for line in io.TextIOWrapper(io.BytesIO(raw), errors="replace"):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                pass
    return sum(vals) / len(vals) if vals else None


n = miss = 0
with open(OUT, "w") as out:
    for i, ((fam, d, tar), members) in enumerate(sorted(by_tar.items()), 1):
        path = SRC / fam / d / f"{tar}.tar.gz"
        if not path.exists():
            miss += len(members); continue
        wanted = {m: did for did, m in members}
        remaining = len(wanted)
        try:
            with tarfile.open(path, "r:gz") as tf:
                for ti in tf:
                    did = wanted.get(os.path.basename(ti.name))
                    if not did:
                        continue
                    f = tf.extractfile(ti)
                    if f:
                        pl = mean_ca_plddt(f.read())
                        if pl is not None:
                            out.write(f"{did}\t{pl:.2f}\n"); n += 1
                    remaining -= 1
                    if remaining == 0:      # stop scanning this archive early
                        break
        except (tarfile.TarError, OSError) as ex:
            print(f"  archive fail {path}: {ex}", file=sys.stderr); miss += len(members)
        if i % 200 == 0:
            print(f"  {i}/{len(by_tar)} archives, {n:,} pLDDT read", flush=True)
print(f"wrote {n:,} per-domain pLDDT ({miss:,} unavailable) -> {OUT}", flush=True)
