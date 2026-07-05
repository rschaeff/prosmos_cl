#!/usr/bin/env python3
"""Geometric determinism of S5 matrix cells vs SSE composition.

For each occupied cell, the within-cell spread of the matched 5-SSE motifs
(median pairwise Kabsch RMSD of fitted SSE-axis endpoints) measures how
determinate the cell's geometry is. Swept across typings within a fixed
skeleton (skeleton-controlled), it tests whether beta-sheet H-bond registration
makes cells geometrically tight: a single unpaired strand is the loosest
composition, then determinism increases monotonically toward EEEEE.

Motif representation: each matched SSE -> two PCA-axis endpoints (denoised,
ordered N->C); 5 SSEs -> 10 points, ordered by query position.

Usage:
  python3 s5_geom_determinism.py \
      --hits-root   ~/work/prosmos_2026/s5_full_pdb_exp/hits \
      --promiscuity ~/work/prosmos_2026/pdb_exp_build/s5_promiscuity.json \
      --out         enum/docs/figures/s5_geom_determinism.png
"""
import argparse, os, re, json
import numpy as np

SEG = re.compile(r"segment-Type:\s*(\S+)\s+Position:\s*(\d+)\s+Range:\s*(\d+)\s*--\s*(\d+)\s+(\S+)\s+Length:\s*(\d+)")
_cache = {}


def dom_path(base, uid):
    p = str(uid).zfill(9)
    return f"{base}/{p[2:7]}/{p}/{p}.pdb"


def load_ca(base, uid):
    if uid in _cache:
        return _cache[uid]
    ca = {}
    try:
        for l in open(dom_path(base, uid)):
            if l[:4] == "ATOM" and l[12:16].strip() == "CA":
                try:
                    ca[(l[21], int(l[22:26]))] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
                except ValueError:
                    pass
    except OSError:
        ca = None
    _cache[uid] = ca
    return ca


def fit_axis(P):
    """CA coords along one SSE -> two PCA-axis endpoints, ordered N->C."""
    P = np.asarray(P, float)
    if len(P) < 3:
        a, b = P[0], P[-1]
    else:
        c = P.mean(0)
        u = np.linalg.svd(P - c, full_matrices=False)[2][0]
        t = (P - c) @ u
        a, b = c + t.min() * u, c + t.max() * u
    if np.linalg.norm(a - P[0]) > np.linalg.norm(b - P[0]):
        a, b = b, a
    return a, b


def motif_axes(ca, segs):
    segs = sorted(segs, key=lambda s: s["position"])[:5]
    if len(segs) != 5:
        return None
    pts = []
    for s in segs:
        cas = [ca[(s["chain"], r)] for r in range(s["start"], s["end"] + 1) if (s["chain"], r) in ca]
        if len(cas) < 3:
            cas = [v for (c, rr), v in ca.items() if s["start"] <= rr <= s["end"]]
        if len(cas) < 2:
            return None
        a, b = fit_axis(cas)
        pts += [a, b]
    return np.array(pts, float)


def parse_hitfile(path):
    segs = []
    try:
        for l in open(path):
            m = SEG.match(l)
            if m:
                segs.append({"position": int(m.group(2)), "start": int(m.group(3)),
                             "end": int(m.group(4)), "chain": m.group(5)})
            elif l.strip() == "END" and segs:
                break
    except OSError:
        return None
    return segs[:5]


def kabsch_rmsd(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1, 1, d]) @ Wt
    return np.sqrt(((Pc @ R - Qc) ** 2).sum(1).mean())


def cell_spread(hits_root, dom_base, sk, ty, n):
    cd = os.path.join(hits_root, f"s5-{sk:04d}-{ty:04d}")
    if not os.path.isdir(cd):
        return None
    motifs = []
    for fn in sorted(os.listdir(cd)):
        if not (fn.startswith("pdb") and fn.endswith(".txt")):
            continue
        uid = int(re.sub(r"\D", "", fn))
        ca = load_ca(dom_base, uid)
        if ca is None:
            continue
        segs = parse_hitfile(os.path.join(cd, fn))
        mp = motif_axes(ca, segs) if segs else None
        if mp is not None:
            motifs.append(mp)
        if len(motifs) >= n:
            break
    if len(motifs) < 5:
        return None
    r = [kabsch_rmsd(motifs[i], motifs[j]) for i in range(len(motifs)) for j in range(i + 1, len(motifs))]
    return float(np.median(r))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits-root", required=True, help="dataset hits/ dir (per-cell hit files)")
    ap.add_argument("--promiscuity", required=True, help="dataset s5_promiscuity.json (for occupancy)")
    ap.add_argument("--domain-dir", default="/data/ecod/af2_pdb_domain_data", help="pre-cut domain PDB store")
    ap.add_argument("--out", required=True, help="output figure PNG")
    ap.add_argument("--n-skels", type=int, default=16)
    ap.add_argument("--n-motifs", type=int, default=50)
    ap.add_argument("--min-hits", type=int, default=100)
    args = ap.parse_args()

    d = json.load(open(args.promiscuity)); row = d["rowSkeleton"]; nH = d["nHits"]
    occ = {}
    for r in range(len(nH)):
        tys = [ty for ty in range(32) if nH[r][ty] >= args.min_hits]
        if tys:
            occ[row[r]] = tys
    skels = sorted(occ, key=lambda s: (-len({bin(t).count('1') for t in occ[s]}), -len(occ[s])))[:args.n_skels]

    data = {}
    for sk in skels:
        data[sk] = {}
        for ty in occ[sk]:
            sp = cell_spread(args.hits_root, args.domain_dir, sk, ty, args.n_motifs)
            if sp is not None:
                data[sk][ty] = sp
    ncells = sum(len(v) for v in data.values())
    print(f"skeletons: {len(skels)}  cells measured: {ncells}  distinct domains read: {len(_cache)}")

    cent = {s: [] for s in range(6)}
    slopes = []
    for sk in skels:
        vals = data[sk]
        if len(vals) < 4:
            continue
        mu = np.mean(list(vals.values()))
        ss, vv = [], []
        for ty, v in vals.items():
            s = bin(ty).count("1"); cent[s].append(v - mu); ss.append(s); vv.append(v)
        mask = [i for i, s in enumerate(ss) if s >= 1]
        if len(set(ss[i] for i in mask)) >= 2:
            slopes.append(np.polyfit([ss[i] for i in mask], [vv[i] for i in mask], 1)[0])

    print("\n=== composition effect, controlling for skeleton (all occupied typings) ===")
    print(f"{'#strands':>8} {'mean_centered(A)':>17} {'n_cells':>8}")
    xs, ys = [], []
    for s in range(6):
        if cent[s]:
            m = np.mean(cent[s]); xs.append(s); ys.append(m)
            print(f"{s:8d} {m:+17.2f} {len(cent[s]):8d}")
    if slopes:
        print(f"\nmean within-skeleton slope (strands 1->5): {np.mean(slopes):+.2f} A/strand (n={len(slopes)})")

    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    base = np.mean([v for sk in data for v in data[sk].values()])
    for sk in skels:
        by = {}
        for t, v in data[sk].items():
            by.setdefault(bin(t).count("1"), []).append(v)
        by = {s: np.mean(v) for s, v in sorted(by.items())}
        ax.plot(list(by), list(by.values()), "-", color="0.75", lw=0.8, alpha=0.7, zorder=1)
    for s in range(6):
        yss = [base + c for c in cent[s]]
        if yss:
            ax.scatter(np.full(len(yss), s) + np.linspace(-0.12, 0.12, len(yss)), yss, s=12,
                       color="#d62728", alpha=0.35, zorder=2)
    ax.plot(xs, np.array(ys) + base, "k-o", lw=2.5, ms=6, zorder=3, label="skeleton-controlled mean")
    ax.set_xlabel("# strands in typing (0 = HHHHH -> 5 = EEEEE)")
    ax.set_ylabel("within-cell motif spread, fitted axes (median pairwise RMSD, A)")
    ax.set_title(f"Geometric determinism vs SSE composition\n"
                 f"{len(skels)} skeletons, {ncells} cells, fitted SSE axes, within-skeleton control")
    ax.set_xticks(range(6)); ax.grid(alpha=0.25); ax.legend()
    plt.tight_layout(); plt.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
