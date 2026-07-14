#!/usr/bin/env python3
"""Fold a sweep's hits/ tree into one compressed table, then verify it.

searchmatrix emits one file per (query, hit record): hits/<query>/pdb<record>.txt.
A full sweep therefore lands tens of millions of ~350-byte files on NFS -- the
space is trivial but the inode count and the write storm are not. This walks the
tree once and emits a single zstd TSV holding everything the files contain, so
the tree can be deleted without losing information.

Each hit file holds one or more MOTIF blocks -- the distinct ways the query
matched that record -- so the row unit is (query, record, motif_index), NOT one
row per file. A file with 90 MOTIF blocks becomes 90 rows.

  query          record          motif  segments
  s5-0003-0027   A0A058YZM2_D11  0      E:1:794-801,E:2:802-813,H:5:825-837,...

Worker mode (one SLURM array task): archive a slice of the query dirs.
    archive_hits.py worker  <hits_root> <out_dir> <shard_idx> <n_shards>
Finalize mode: concatenate shards, emit derived summaries, verify.
    archive_hits.py finalize <hits_root> <out_dir> <n_shards>
"""
import os
import sys
import re
import glob
import subprocess

SEG = re.compile(
    r"segment-Type:\s*(\S+)\s+Position:\s*(\d+)\s+Range:\s*(\d+)\s*--\s*(\d+)"
)


def parse_hit_file(path):
    """Yield one compact segment-string per MOTIF block in the file."""
    blocks, cur = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("MOTIF:"):
                if cur:
                    blocks.append(cur)
                cur = []
            else:
                m = SEG.match(line.strip())
                if m:
                    t, pos, lo, hi = m.groups()
                    cur.append(f"{t}:{int(pos)}:{int(lo)}-{int(hi)}")
    if cur:
        blocks.append(cur)
    return [",".join(b) for b in blocks if b]


def worker(hits_root, out_dir, shard, n_shards):
    qdirs = sorted(os.listdir(hits_root))
    mine = qdirs[shard::n_shards]
    os.makedirs(out_dir, exist_ok=True)
    tsv = os.path.join(out_dir, f"part_{shard:04d}.tsv.zst")
    cnt = os.path.join(out_dir, f"count_{shard:04d}.tsv")

    proc = subprocess.Popen(
        ["zstd", "-3", "-q", "-o", tsv, "-f", "-"],
        stdin=subprocess.PIPE, text=True,
    )
    n_files = n_rows = 0
    with open(cnt, "w") as cfh:
        for q in mine:
            qpath = os.path.join(hits_root, q)
            try:
                files = os.listdir(qpath)
            except OSError:
                continue
            nrec = 0
            for fn in files:
                if not fn.endswith(".txt"):
                    continue
                rec = fn[:-4]
                if rec.startswith("pdb"):
                    rec = rec[3:]
                blocks = parse_hit_file(os.path.join(qpath, fn))
                n_files += 1
                nrec += 1
                for i, seg in enumerate(blocks):
                    proc.stdin.write(f"{q}\t{rec}\t{i}\t{seg}\n")
                    n_rows += 1
            # per-query record count: the unit the dark-fraction analysis uses
            cfh.write(f"{q}\t{nrec}\n")
    proc.stdin.close()
    rc = proc.wait()
    with open(os.path.join(out_dir, f"stat_{shard:04d}.tsv"), "w") as sfh:
        sfh.write(f"{shard}\t{len(mine)}\t{n_files}\t{n_rows}\t{rc}\n")
    print(f"shard {shard}: {len(mine)} qdirs, {n_files} files, {n_rows} rows, zstd rc={rc}")
    return rc


def finalize(hits_root, out_dir, n_shards):
    stats = sorted(glob.glob(os.path.join(out_dir, "stat_*.tsv")))
    if len(stats) != n_shards:
        sys.exit(f"FAIL: {len(stats)} shard stats, expected {n_shards}")
    tot_q = tot_f = tot_r = 0
    for s in stats:
        _, nq, nf, nr, rc = open(s).read().strip().split("\t")
        if rc != "0":
            sys.exit(f"FAIL: {s} reports zstd rc={rc}")
        tot_q += int(nq); tot_f += int(nf); tot_r += int(nr)

    n_qdirs = len(os.listdir(hits_root))
    if tot_q != n_qdirs:
        sys.exit(f"FAIL: shards covered {tot_q} query dirs, tree has {n_qdirs}")

    # single archive
    parts = sorted(glob.glob(os.path.join(out_dir, "part_*.tsv.zst")))
    archive = os.path.join(out_dir, "hits.tsv.zst")
    with open(archive, "wb") as out:
        for p in parts:
            with open(p, "rb") as fh:
                out.write(fh.read())          # zstd frames concatenate cleanly

    # derived summaries, straight from the shard counts (no second tree walk)
    counts = {}
    for c in sorted(glob.glob(os.path.join(out_dir, "count_*.tsv"))):
        for line in open(c):
            q, n = line.rstrip("\n").split("\t")
            counts[q] = int(n)
    with open(os.path.join(out_dir, "query_counts.tsv"), "w") as fh:
        fh.write("query\thit_records\n")
        for q in sorted(counts):
            fh.write(f"{q}\t{counts[q]}\n")

    # verify the archive reads back with exactly the row count we wrote
    n_read = 0
    recs = set()
    p = subprocess.Popen(["zstd", "-dc", archive], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) == 4:
            n_read += 1
            recs.add(f[1])
    p.wait()
    if n_read != tot_r:
        sys.exit(f"FAIL: archive reads back {n_read} rows, shards wrote {tot_r}")

    with open(os.path.join(out_dir, "distinct_hitters.txt"), "w") as fh:
        for r in sorted(recs):
            fh.write(r + "\n")

    size = os.path.getsize(archive)
    print(f"VERIFIED")
    print(f"  query dirs      : {n_qdirs}")
    print(f"  hit files       : {tot_f:,}")
    print(f"  motif rows      : {tot_r:,}")
    print(f"  distinct records: {len(recs):,}")
    print(f"  archive         : {archive} ({size/1e9:.2f} GB)")
    print(f"  queries w/ hits : {sum(1 for v in counts.values() if v)} / {len(counts)}")
    print(f"SAFE TO DELETE {hits_root} -- {tot_f:,} files fold into 1")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "worker":
        sys.exit(worker(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])))
    elif mode == "finalize":
        finalize(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    else:
        sys.exit(__doc__)
