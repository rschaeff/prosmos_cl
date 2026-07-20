#!/usr/bin/env python3
"""Rebuild afdb_assigned from the CORRECTED AFDB sweep (s5_inv), accession-keyed.

Why a rebuild and not a re-run: the old dataset resolved protein-level hits
(pdbdpam_<UniProt>.txt) to numeric ecod_af2_pdb uids. The corrected sweep is
DOMAIN-level (pdb<acc>_D<n>.txt) over the afdb_200m universe, and only ~1.5% of
its accessions exist in ecod_af2_pdb — so there is no uid to key on. We therefore
key exemplars by UniProt accession + DPAM domain, exactly like afdb_unassigned:
did=<acc>_D<n>, unp=<acc>, ecodUid=0 (viewer loads /api/structure/<unp>).

Unlike afdb_unassigned these domains DO have ECOD T-groups (afdb_name_tgroup.tsv),
so nT / nH / group are real.

Inputs
  s5_grid/afdb_rare/*.pkl      (cell_id, record15) pairs from the corrected tree
  s5_grid/afdb_name_tgroup.tsv record -> T-group
  s5_inv/hits/<cell>/pdb<rec>.txt   read ONLY for chosen exemplars (segments)
Writes data/afdb_assigned/{cells.json, cell/*.json}
"""
import json, os, re, sys, glob, pickle
from pathlib import Path
from collections import defaultdict, Counter

REPO = Path.home() / "dev/prosmos_cl"
G = Path.home() / "work/prosmos_2026/s5_grid"
HITS = Path.home() / "work/prosmos_2026/s5_inv/hits"
STRUCT = Path.home() / "work/prosmos_2026/afdb_struct"
OUT = Path.home() / "dev/prosmos_inspect/data/afdb_assigned"
K, M_PER_GROUP = 15, 3
sys.path.insert(0, str(REPO / "enum" / "src"))
sys.path.insert(0, str(REPO / "enum" / "scripts"))
from ssp_enum.enumerate import enumerate_skeletons          # noqa: E402
from plot_skeleton_schematic import compute_matrix          # noqa: E402
SKELS = enumerate_skeletons(5)
seg_re = re.compile(r"segment-Type:\s*(\S+)\s+Position:\s*(\d+)\s+Range:\s*(\d+)\s*--\s*(\d+)\s+(\S+)\s+Length:\s*(\d+)")


def tystr(ty):
    return "".join("E" if (ty >> b) & 1 else "H" for b in range(4, -1, -1))


def skeleton_geom(sk):
    s = SKELS[sk]; adjm = s.adjacency_matrix()
    return {"sk": sk, "nodes": [{"q": p.q, "r": p.r} for p in s.points],
            "adj": [[bool(adjm[i][j]) for j in range(5)] for i in range(5)],
            "orientations": list(s.orientations)}


def hgroup(t):
    p = t.split("."); return ".".join(p[:2]) if len(p) >= 2 else t


def first_motif(sk, ty, rec):
    p = HITS / f"s5-{sk:04d}-{ty:04d}" / f"pdb{rec}.txt"
    segs = []
    try:
        for ln in p.open():
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
    tg = {}
    for line in open(G / "afdb_name_tgroup.tsv"):
        a, b = line.rstrip("\n").split("\t")
        tg[a] = b
    pairs = []
    for f in glob.glob(str(G / "afdb_rare" / "*.pkl")):
        pairs.extend(pickle.load(open(f, "rb")))
    print(f"tgroups {len(tg):,}  (cell,record) pairs {len(pairs):,}", flush=True)

    cell_recs = defaultdict(list)
    for cid, r in pairs:
        cell_recs[cid].append(r)
    have = {p.stem for p in STRUCT.glob("*.pdb")}
    print(f"cells with hits {len(cell_recs):,}  imported models on disk {len(have):,}", flush=True)

    (OUT / "cell").mkdir(parents=True, exist_ok=True)
    index = []
    n_ex_total = n_with_struct = 0
    for cid, recs in sorted(cell_recs.items()):
        sk, ty = cid // 32, cid % 32
        tcnt, hcnt = Counter(), Counter()
        by_t = defaultdict(list)
        for r in recs:
            t = tg.get(r)
            if not t:
                continue
            tcnt[t] += 1; hcnt[hgroup(t)] += 1
            by_t[t].append(r)
        nHits = len(set(recs))
        nT, nH = len(tcnt), len(hcnt)
        typing = tystr(ty)
        exs = []
        for t, cnt in tcnt.most_common(K):
            # prefer domains whose AF model is already on disk so the viewer works
            cand = sorted(set(by_t[t]),
                          key=lambda r: (r.split("_D")[0] not in have, r))
            for rec in cand[:M_PER_GROUP]:
                segs = first_motif(sk, ty, rec)
                if not segs:
                    continue
                acc = rec.split("_D")[0].replace("dpam_", "")
                exs.append({"did": rec, "ecodUid": 0, "pdb": None, "unp": acc,
                            "chain": segs[0]["chain"], "group": t,
                            "count": cnt, "segments": segs})
                if acc in have:
                    n_with_struct += 1
        n_ex_total += len(exs)
        view = "unitypical" if nT == 1 else "promiscuous"
        detail = {"sk": sk, "ty": ty, "typing": typing, "view": view,
                  "nT": nT, "nH": nH, "nHits": nHits,
                  "matrix": {"typing": typing, "codes": compute_matrix(SKELS[sk], tuple(typing))},
                  "skeleton": skeleton_geom(sk), "exemplars": exs, "stats": None}
        (OUT / "cell" / f"{sk:04d}_{ty:02d}.json").write_text(json.dumps(detail))
        index.append({"sk": sk, "ty": ty, "typing": typing, "view": view,
                      "nT": nT, "nH": nH, "nHits": nHits})
    index.sort(key=lambda r: (-r["nHits"], r["sk"], r["ty"]))
    (OUT / "cells.json").write_text(json.dumps(index))
    print(f"wrote {len(index):,} cells  exemplars {n_ex_total:,}  "
          f"({n_with_struct:,} with a model already on disk) -> {OUT}")


if __name__ == "__main__":
    main()
