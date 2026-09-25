#!/usr/bin/env python3
"""Extract the AFDB domain structures a census needs, from both source trees.

The AFDB store keeps domains as tarballs under two trees that name their members
DIFFERENTLY, plus a flat gzip directory:

  new/{shard}/{KEY}.tar.gz    members `{ACC}_D{n}.pdb`        91,619 tarballs, ~99% of domains
  dpam/{shard}/{KEY}.tar.gz   members `dpam_{ACC}_nD{n}.pdb`   6,977 tarballs, ~1%
  ecod/ecod_{UID:09d}.pdb.gz  flat, gzipped, one per file      7,746 files

Globbing only `dpam/` returns ~5% of the structures and looks like missing data rather
than a glob bug -- that cost a cancelled 200-task array on 2026-08-07. This script globs
both trees and indexes both spellings.

KEY is `acc[-4:-2]`, but the SHARD IS NOT DERIVABLE from the accession, and with ~80
tarballs sharing each KEY a per-record search would scan ~80 archives. So this does the
transpose: one streaming pass per tarball, extracting every needed member found in it.
Members absent from a tarball are simply not extracted -- no error, unlike `tar -T`.

Output is a flat directory of `<record>.pdb`, which is what `acceptance.py --structures`
expects. Copies out; never symlinks (the source is another user's home and will move).

Usage:
  extract_afdb_structures.py --records recs.txt --out structures/ [--shard K/N]

`--records` holds one record name per line, without `.ssd` -- e.g. the record column of a
sieve output table. Shard over tarballs for a SLURM array; every shard reads the same
record list and writes into the same output directory.
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import sys
import tarfile

ROOT = "/home/grey/afdb.200m/non_singleton_4p9m_structures"


def load_records(path):
    recs = set()
    with open(path) as fh:
        for line in fh:
            r = line.strip()
            if not r or r.startswith("#"):
                continue
            if r.endswith(".ssd"):
                r = r[:-4]
            recs.add(r)
    return recs


def key_of(record):
    """KEY = acc[-4:-2]. Handles both member spellings."""
    acc = record
    if acc.startswith("dpam_"):
        acc = acc[5:]
        acc = acc.rsplit("_nD", 1)[0]
    else:
        acc = acc.rsplit("_D", 1)[0]
    return acc[-4:-2] if len(acc) >= 4 else None


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--records", required=True, help="one record name per line")
    p.add_argument("--out", required=True, help="output directory of <record>.pdb")
    p.add_argument("--tarballs", default=None,
                   help="precomputed tarball list (default: walk both trees)")
    p.add_argument("--root", default=ROOT)
    p.add_argument("--shard", default=None, metavar="K/N", help="process tarballs K, K+N, ...")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()

    recs = load_records(a.records)
    want = {r + ".pdb": r for r in recs}
    os.makedirs(a.out, exist_ok=True)

    # the flat ecod/ tree: experimental reference domains, gzipped one per file
    ecod = sorted(r for r in recs if r.startswith("ecod_"))
    n_ecod = 0
    if ecod and (a.shard is None or a.shard.split("/")[0] == "0"):
        for r in ecod:
            src = os.path.join(a.root, "ecod", r + ".pdb.gz")
            dst = os.path.join(a.out, r + ".pdb")
            if os.path.exists(dst) or not os.path.exists(src):
                continue
            if not a.dry_run:
                with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
                    shutil.copyfileobj(fi, fo)
            n_ecod += 1

    if a.tarballs:
        tars = [l.strip() for l in open(a.tarballs) if l.strip()]
    else:
        tars = []
        for tree in ("new", "dpam"):
            for dirpath, _, files in os.walk(os.path.join(a.root, tree)):
                tars.extend(os.path.join(dirpath, f) for f in files if f.endswith(".tar.gz"))
        tars.sort()

    # only tarballs whose KEY some wanted record could live in
    keys = {key_of(r) for r in recs if not r.startswith("ecod_")} - {None}
    tars = [t for t in tars if os.path.basename(t)[:-len(".tar.gz")] in keys]

    if a.shard:
        k, n = (int(x) for x in a.shard.split("/"))
        tars = tars[k::n]

    found = 0
    scanned = 0
    for t in tars:
        scanned += 1
        try:
            with tarfile.open(t, "r:gz") as tf:
                for m in tf:
                    if not m.isfile():
                        continue
                    base = os.path.basename(m.name)
                    rec = want.get(base)
                    if rec is None:
                        continue
                    dst = os.path.join(a.out, rec + ".pdb")
                    if os.path.exists(dst):
                        continue
                    found += 1
                    if a.dry_run:
                        continue
                    src = tf.extractfile(m)
                    if src is None:
                        continue
                    tmp = dst + ".part%d" % os.getpid()
                    with open(tmp, "wb") as fo:
                        shutil.copyfileobj(src, fo)
                    os.replace(tmp, dst)
        except (tarfile.TarError, OSError) as e:
            print(f"warning: {t}: {e}", file=sys.stderr)
    print(f"records wanted {len(recs):,}; tarballs scanned {scanned:,} of {len(keys):,} keys; "
          f"extracted {found:,} from tarballs, {n_ecod:,} from ecod/")


if __name__ == "__main__":
    main()
