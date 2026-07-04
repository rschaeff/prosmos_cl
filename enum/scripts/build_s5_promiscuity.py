#!/usr/bin/env python3
"""Per-cell ECOD-group promiscuity for the S5 matrix.

For each (skeleton, typing) cell, how many distinct ECOD H-groups and T-groups
do the hit domains span?  nT==1 -> cell is UNITYPICAL (that geometry uniquely
tags one topology); large nT -> PROMISCUOUS (a common motif shared across folds).

Also records, per cell, the group->exemplar-domains breakdown so we can show an
exemplar (unitypical) or a gallery (promiscuous) downstream.

Inputs:  manual_hits.txt, ecod_uid_class.tsv (ecod_uid, did, t_id, f_id),
         s5_timeseries.json (for the shared fixed row order)
Output:  s5_promiscuity.json
"""
import json, re
from pathlib import Path
from collections import defaultdict, Counter

SP = Path("/tmp/claude-1219/-home-rschaeff-dev-prosmos-cl/02f20625-920b-46f1-a2bf-bc06d84727af/scratchpad")
NSK, NTY = 198, 32
hit_re = re.compile(r"s5-(\d{4})-(\d{4})/pdb0*(\d+)\.txt$")
# For the exemplar gallery: keep up to M_STORE experimental domains for each of
# the top TOP_GROUPS T-groups per cell (breadth = groups, depth = domains/group).
TOP_GROUPS, M_STORE = 20, 6


def hgroup(tid):  # "11.1.1" -> "11.1"
    p = tid.split(".")
    return ".".join(p[:2]) if len(p) >= 2 else tid


def main():
    # euid -> (did, H, T)
    cls = {}
    for line in (SP / "ecod_uid_class.tsv").open():
        euid, did, tid, fid = line.rstrip("\n").split("\t")
        if not tid:
            continue
        cls[int(euid)] = (did, hgroup(tid), tid)

    # per cell: euid set
    cell_euids = defaultdict(set)
    n = 0
    for line in (SP / "manual_hits.txt").open():
        m = hit_re.search(line)
        if not m:
            continue
        sk, ty, euid = int(m.group(1)), int(m.group(2)), int(m.group(3))
        cell_euids[(sk, ty)].add(euid)
        n += 1

    # shared row order from the time-series build
    ts = json.load((SP / "s5_timeseries.json").open())
    row_order = ts["rowSkeleton"]
    row_pos = {sk: i for i, sk in enumerate(row_order)}

    # per-cell metrics + group breakdown
    nH = [[0] * NTY for _ in range(NSK)]
    nT = [[0] * NTY for _ in range(NSK)]
    nHits = [[0] * NTY for _ in range(NSK)]
    cell_groups = {}   # "r,ty" -> {"H":[[hg,cnt,exemplarDid,exemplarEuid],...sorted], "T":[...], "nhit":..}
    for (sk, ty), euids in cell_euids.items():
        r = row_pos[sk]
        hcnt = Counter()
        tcnt = Counter()
        h_ex = {}
        t_ex = {}
        t_doms = defaultdict(list)   # T-group -> [(did, euid), ...] experimental only
        for e in euids:
            if e not in cls:
                continue
            did, H, T = cls[e]
            hcnt[H] += 1
            tcnt[T] += 1
            # exemplar = first PDB-style domain we see (prefer e-prefixed experimental)
            if H not in h_ex or (did.startswith("e") and not h_ex[H][0].startswith("e")):
                h_ex[H] = (did, e)
            if T not in t_ex or (did.startswith("e") and not t_ex[T][0].startswith("e")):
                t_ex[T] = (did, e)
            if did.startswith("e"):
                t_doms[T].append((did, e))
        nH[r][ty] = len(hcnt)
        nT[r][ty] = len(tcnt)
        nHits[r][ty] = len(euids)
        Hlist = [[g, c, h_ex[g][0], h_ex[g][1]] for g, c in hcnt.most_common()]
        Tlist = [[g, c, t_ex[g][0], t_ex[g][1]] for g, c in tcnt.most_common()]
        # depth: up to M_STORE domains for each of the top TOP_GROUPS T-groups
        Tdoms = {g: [[d, u] for d, u in sorted(set(t_doms[g]))[:M_STORE]]
                 for g, _c in tcnt.most_common(TOP_GROUPS) if t_doms[g]}
        cell_groups[f"{r},{ty}"] = {"nhit": len(euids), "H": Hlist, "T": Tlist,
                                    "Tdoms": Tdoms}

    # distributions
    flatT = [nT[r][ty] for r in range(NSK) for ty in range(NTY) if nHits[r][ty] > 0]
    flatH = [nH[r][ty] for r in range(NSK) for ty in range(NTY) if nHits[r][ty] > 0]
    occ = len(flatT)
    uni_T = sum(1 for v in flatT if v == 1)
    uni_H = sum(1 for v in flatH if v == 1)

    out = {
        "nsk": NSK, "nty": NTY,
        "rowSkeleton": row_order,
        "nH": nH, "nT": nT, "nHits": nHits,
        "cellGroups": cell_groups,
        "stats": {
            "occupied": occ,
            "unitypicalT": uni_T, "unitypicalH": uni_H,
            "maxT": max(flatT), "maxH": max(flatH),
            "medianT": sorted(flatT)[len(flatT) // 2],
        },
    }
    (SP / "s5_promiscuity.json").write_text(json.dumps(out, separators=(",", ":")))
    print(f"occupied cells: {occ}")
    print(f"unitypical (nT==1): {uni_T} ({100*uni_T/occ:.0f}%)   "
          f"unitypical (nH==1): {uni_H} ({100*uni_H/occ:.0f}%)")
    print(f"nT: median {out['stats']['medianT']}, max {out['stats']['maxT']}")
    # histogram of nT
    from collections import Counter as C
    hist = C(flatT)
    print("nT distribution (bucketed):")
    buckets = [(1, 1), (2, 2), (3, 5), (6, 10), (11, 20), (21, 50), (51, 10**9)]
    for lo, hi in buckets:
        c = sum(v for k, v in hist.items() if lo <= k <= hi)
        lab = f"{lo}" if lo == hi else (f"{lo}+" if hi > 10**8 else f"{lo}-{hi}")
        print(f"  nT {lab:>6}: {c:5d} cells")


if __name__ == "__main__":
    main()
