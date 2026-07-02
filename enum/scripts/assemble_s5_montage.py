#!/usr/bin/env python3
"""Assemble curated exemplar renders into labeled figures:
   - s5_exemplars_unitypical.png : 2x2, one fold-specific geometry each
   - s5_montage_skXXXX_tyNN.png   : promiscuous cell, 6 folds sharing one geometry
"""
import json, math, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

SP = Path("/tmp/claude-1219/-home-rschaeff-dev-prosmos-cl/02f20625-920b-46f1-a2bf-bc06d84727af/scratchpad")
R = SP / "renders"
SPEC = json.load((SP / "render_spec.json").open())
OUT = SP

# skeleton geometry + ProSMoS matrix (same source as plot_skeleton_schematic.py)
REPO = Path.home() / "dev/prosmos_cl"
sys.path.insert(0, str(REPO / "enum" / "src"))
sys.path.insert(0, str(REPO / "enum" / "scripts"))
from ssp_enum.enumerate import enumerate_skeletons  # noqa: E402
from plot_skeleton_schematic import compute_matrix, draw_matrix  # noqa: E402
SKELS = enumerate_skeletons(5)
# rainbow N->C, matching the SSE highlight colors in the PyMOL renders
NODE_COLORS = ["#0a6cd4", "#2ca02c", "#e6c700", "#ff8c00", "#e10600"]


def hex_to_xy(q, r):
    return q + 0.5 * r, -r * (math.sqrt(3) / 2)


def draw_schematic(ax, sk, typing):
    """Draw skeleton `sk` on the hex lattice; nodes rainbow N->C, ○=H □=E."""
    skel = SKELS[sk]
    pts = [hex_to_xy(p.q, p.r) for p in skel.points]
    adj = skel.adjacency_matrix()
    n = skel.dim
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i][j]:
                ax.plot([pts[i][0], pts[j][0]], [pts[i][1], pts[j][1]],
                        "-", color="#888", lw=2.2, zorder=1)
    for idx, (x, y) in enumerate(pts):
        marker = "o" if typing[idx] == "H" else "s"
        ax.scatter([x], [y], s=560, c=NODE_COLORS[idx], edgecolors="#222",
                   linewidths=1.6, marker=marker, zorder=3)
        ax.text(x, y, str(idx + 1), ha="center", va="center",
                fontsize=11, fontweight="bold", color="white", zorder=4)
    xs = [p[0] for p in pts]
    ax.text(min(xs) - 1.2, 0, "shared\ngeometry →", ha="right", va="center",
            fontsize=9.5, color="#333")
    ax.text(max(xs) + 1.1, 0, "○ helix\n□ strand", ha="left", va="center",
            fontsize=9, color="#555")
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(min(xs) - 3.0, max(xs) + 2.6)
    ax.set_ylim(-1.2, 1.2)
    for s in ax.spines.values():
        s.set_visible(False)


def panel(ax, png, title, sub=None):
    ax.imshow(mpimg.imread(str(R / png)))
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, fontsize=10, pad=2)
    if sub:
        ax.set_xlabel(sub, fontsize=8.5, color="#555")


# ---- unitypical 2x2 --------------------------------------------------------
uni = SPEC["unitypical"]
fig, axes = plt.subplots(2, 2, figsize=(8.5, 8.8))
for ax, j in zip(axes.flat, uni):
    e = j["exemplar"]
    panel(ax, f"uni_sk{j['sk']:04d}_ty{j['ty']:02d}_{e['did']}.png",
          f"{e['did']}  ·  T-group {e['group']}",
          f"skeleton {j['sk']:04d} · typing {j['typing']} · all {j['nhit']} hits in ONE topology")
for ax in axes.flat[len(uni):]:
    ax.axis("off")
fig.suptitle("Unitypical S5 cells — the 5-SSE geometry (rainbow N→C) uniquely tags one ECOD fold",
             fontsize=12, y=0.99)
plt.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.05, hspace=0.28, wspace=0.06)
plt.savefig(OUT / "s5_exemplars_unitypical.png", dpi=150)
plt.close()
print("wrote s5_exemplars_unitypical.png")

# ---- promiscuous montages --------------------------------------------------
for j in SPEC["promiscuous"]:
    ex = j["exemplars"]
    n = len(ex)
    fig = plt.figure(figsize=(11.5, 9.6))
    # top banner: skeleton schematic + ProSMoS query matrix; below: 2x3 structures
    outer = fig.add_gridspec(2, 1, height_ratios=[0.72, 4], hspace=0.16)
    banner = outer[0].subgridspec(1, 2, width_ratios=[1.5, 1], wspace=0.12)
    draw_schematic(fig.add_subplot(banner[0, 0]), j["sk"], j["typing"])
    types = tuple(j["typing"])
    mat = compute_matrix(SKELS[j["sk"]], types)
    draw_matrix(fig.add_subplot(banner[0, 1]), mat, types,
                title="ProSMoS query matrix")
    grid = outer[1].subgridspec(2, 3, hspace=0.3, wspace=0.05)
    for k, e in enumerate(ex):
        ax = fig.add_subplot(grid[k // 3, k % 3])
        png = f"pro_sk{j['sk']:04d}_ty{j['ty']:02d}_{e['group']}_{e['did']}.png"
        panel(ax, png, f"{e['did']}  ·  T-group {e['group']}",
              f"{e['count']} rep domains in this group")
    fig.suptitle(f"Promiscuous S5 cell — skeleton {j['sk']:04d}, typing {j['typing']} "
                 f"(one geometry shared by {j['nT']} topologies; 6 shown)\n"
                 f"same 5-SSE arrangement (rainbow N→C), different folds",
                 fontsize=12, y=0.995)
    plt.subplots_adjust(left=0.03, right=0.97, top=0.9, bottom=0.04)
    name = f"s5_montage_sk{j['sk']:04d}_ty{j['ty']:02d}.png"
    plt.savefig(OUT / name, dpi=150)
    plt.close()
    print("wrote", name)
