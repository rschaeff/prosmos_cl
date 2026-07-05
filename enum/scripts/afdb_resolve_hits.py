#!/usr/bin/env python3
"""Resolve protein-level AFDB ProSMoS hits to their DPAM domain.

AFDB metamatrix entries are whole AlphaFold models (hit files `pdbdpam_<acc>.txt`,
where <acc> is either a bare UniProt for a protein-level entry or `<UniProt>_nD<k>`
for a DPAM-domain entry). A hit's matched 5-SSE ranges fall within one DPAM domain
of that protein. Each hit is assigned to a domain uid so that its ECOD T-group,
pre-cut structure file, and contacts all resolve through the SAME pipeline as the
pdb_exp dataset:

  - domain-level entry (`_nD<k>`): look up (unp_acc, k) -> uid directly.
  - protein-level entry (bare UniProt): pick the domain whose ECOD residue range
    overlaps the matched segments most.

Classification comes from ecod_af2_pdb.domain (uid, unp_acc, id, tid, range).

Outputs (--out-dir):
  afdb_hits_resolved.txt  lines 's5-SK-TY/pdb<uid09d>.txt'  (feeds build_promiscuity)
  afdb_segments.json      {uid: [segment,...]} matched motif per resolved domain
  resolve_stats.txt       coverage / ambiguity report

Usage:
  python3 afdb_resolve_hits.py \
      --hitroots ~/work/prosmos_2026/s5_full_afdb/hits \
                 ~/work/prosmos_2026/s5_full_afdb_missing/hits \
      --out-dir  ~/work/prosmos_2026/afdb_build
"""
import argparse, json, os, re, sys
from collections import defaultdict
import psycopg2

seg_re = re.compile(r"segment-Type:\s*(\S+)\s+Position:\s*(\d+)\s+Range:\s*(\d+)\s*--\s*(\d+)\s+(\S+)\s+Length:\s*(\d+)")
cell_re = re.compile(r"s5-(\d{4})-(\d{4})$")
ND_RE = re.compile(r"^(.*)_nD(\d+)$")   # hit acc 'K7KP71_nD2' -> base 'K7KP71', idx 2
ID_ND_RE = re.compile(r"_nD(\d+)$")      # DB id 'K7KP71_F1_nD2' -> idx 2


def parse_range(s):
    """'161-320,356-370' -> [(161,320),(356,370)]"""
    out = []
    for part in s.split(","):
        m = re.match(r"(-?\d+)-(-?\d+)", part.strip())
        if m:
            out.append((int(m.group(1)), int(m.group(2))))
    return out


def overlap(a, b):
    """total residue overlap between two interval lists"""
    t = 0
    for lo1, hi1 in a:
        for lo2, hi2 in b:
            t += max(0, min(hi1, hi2) - max(lo1, lo2) + 1)
    return t


def first_motif(path):
    segs = []
    try:
        for ln in open(path):
            m = seg_re.match(ln)
            if m:
                segs.append({"type": m.group(1), "position": int(m.group(2)),
                             "start": int(m.group(3)), "end": int(m.group(4)),
                             "chain": m.group(5), "length": int(m.group(6))})
            elif ln.strip() == "END" and segs:
                break
    except OSError:
        return []
    return segs[:5]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hitroots", nargs="+", required=True, help="AFDB hits/ dirs (one or more)")
    ap.add_argument("--out-dir", required=True, help="output directory")
    ap.add_argument("--db-host", default="sangala")
    ap.add_argument("--db-port", type=int, default=45000)
    ap.add_argument("--db-name", default="ecod_af2_pdb")
    ap.add_argument("--db-user", default="ecod")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # pass 1: read every hit file -> (cell, base_acc, nd_idx|None, segments)
    hits, bases, n_files = [], set(), 0
    for root in args.hitroots:
        for cell in os.listdir(root):
            cm = cell_re.match(cell)
            if not cm:
                continue
            sk, ty = int(cm.group(1)), int(cm.group(2))
            cdir = os.path.join(root, cell)
            for fn in os.listdir(cdir):
                if not fn.startswith("pdbdpam_") or not fn.endswith(".txt"):
                    continue
                acc = fn[len("pdbdpam_"):-4]
                nd = ND_RE.match(acc)
                base, idx = (nd.group(1), int(nd.group(2))) if nd else (acc, None)
                segs = first_motif(os.path.join(cdir, fn))
                if not segs:
                    continue
                hits.append((sk, ty, base, idx, segs))
                bases.add(base)
                n_files += 1
                if n_files % 20000 == 0:
                    print(f"  read {n_files} hits...", file=sys.stderr)
    print(f"read {n_files} hits, {len(bases)} distinct UniProt", file=sys.stderr)

    # pass 2: all domains for those UniProts -> by protein (range) and by (acc, nD idx)
    conn = psycopg2.connect(host=args.db_host, port=args.db_port, dbname=args.db_name, user=args.db_user)
    cur = conn.cursor()
    by_prot, by_idx = defaultdict(list), {}
    base_list = list(bases)
    for i in range(0, len(base_list), 5000):
        cur.execute("select uid, unp_acc, id, tid, range from domain "
                    "where unp_acc = any(%s) and tid is not null", (base_list[i:i + 5000],))
        for uid, acc, did, tid, rng in cur.fetchall():
            by_prot[acc].append((uid, tid, parse_range(rng or "")))
            m = ID_ND_RE.search(did or "")
            if m:
                by_idx[(acc, int(m.group(1)))] = uid
    print(f"got domains for {len(by_prot)} UniProts", file=sys.stderr)

    # pass 3: resolve. domain-level -> direct; protein-level -> range overlap.
    resolved, segmap = [], {}
    n_ok = n_noprot = n_noov = 0
    for sk, ty, base, idx, segs in hits:
        uid = by_idx.get((base, idx)) if idx is not None else None
        if uid is None:
            doms = by_prot.get(base)
            if not doms:
                n_noprot += 1
                continue
            span = [(s["start"], s["end"]) for s in segs]
            best_ov = 0
            for u, tid, ivs in doms:
                ov = overlap(span, ivs)
                if ov > best_ov:
                    best_ov, uid = ov, u
            if uid is None:
                n_noov += 1
                continue
        resolved.append((sk, ty, uid))
        segmap[str(uid)] = segs
        n_ok += 1

    with open(os.path.join(args.out_dir, "afdb_hits_resolved.txt"), "w") as f:
        for sk, ty, uid in resolved:
            f.write(f"s5-{sk:04d}-{ty:04d}/pdb{uid:09d}.txt\n")
    json.dump(segmap, open(os.path.join(args.out_dir, "afdb_segments.json"), "w"))
    rep = (f"hits read: {n_files}\nresolved: {n_ok}\n"
           f"unresolved (protein not in {args.db_name}): {n_noprot}\n"
           f"unresolved (no domain overlaps match): {n_noov}\n"
           f"distinct resolved domains: {len(segmap)}\n")
    open(os.path.join(args.out_dir, "resolve_stats.txt"), "w").write(rep)
    print(rep)


if __name__ == "__main__":
    main()
