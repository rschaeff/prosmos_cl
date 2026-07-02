#!/usr/bin/env python3
"""Phase C: compose one composite figure per occupied cell.
  banner (skeleton schematic + ProSMoS query matrix) + exemplar render(s).
Array-capable: NTASKS / SLURM_ARRAY_TASK_ID slice the cell list.
Output: s5_gallery/cells/{sk:04d}_{ty:02d}.png
"""
import json, math, os, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

GAL = Path("/home/rschaeff/work/prosmos_2026/s5_gallery")
REPO = Path.home() / "dev/prosmos_cl"
RND = GAL / "renders"
OUT = GAL / "cells"; OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(REPO / "enum" / "src"))
sys.path.insert(0, str(REPO / "enum" / "scripts"))
from ssp_enum.enumerate import enumerate_skeletons          # noqa: E402
from plot_skeleton_schematic import compute_matrix, draw_matrix  # noqa: E402
SKELS = enumerate_skeletons(5)
NODE_COLORS = ["#0a6cd4", "#2ca02c", "#e6c700", "#ff8c00", "#e10600"]


def hex_to_xy(q, r):
    return q + 0.5 * r, -r * (math.sqrt(3) / 2)


def draw_schematic(ax, sk, typing):
    skel = SKELS[sk]
    pts = [hex_to_xy(p.q, p.r) for p in skel.points]
    adj = skel.adjacency_matrix(); n = skel.dim
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i][j]:
                ax.plot([pts[i][0], pts[j][0]], [pts[i][1], pts[j][1]], "-",
                        color="#888", lw=2.0, zorder=1)
    for idx, (x, y) in enumerate(pts):
        mk = "o" if typing[idx] == "H" else "s"
        ax.scatter([x], [y], s=430, c=NODE_COLORS[idx], edgecolors="#222",
                   linewidths=1.4, marker=mk, zorder=3)
        ax.text(x, y, str(idx + 1), ha="center", va="center", fontsize=9.5,
                fontweight="bold", color="white", zorder=4)
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    cy = (min(ys) + max(ys)) / 2
    ax.text(max(xs) + 1.0, cy, "○ H\n□ E", ha="left", va="center", fontsize=8, color="#555")
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(min(xs) - 1.0, max(xs) + 2.4)
    ax.set_ylim(min(ys) - 0.8, max(ys) + 0.8)
    for s in ax.spines.values():
        s.set_visible(False)


def struct_panel(ax, png, title, sub):
    p = RND / png
    if p.exists():
        ax.imshow(mpimg.imread(str(p)))
    else:
        ax.text(0.5, 0.5, "(render missing)", ha="center", va="center", fontsize=8, color="#999")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, fontsize=9, pad=2)
    ax.set_xlabel(sub, fontsize=8, color="#555")


def compose(c):
    sk, ty, typing, view = c["sk"], c["ty"], c["typing"], c["view"]
    ex = c["exemplars"]
    ncol = max(1, len(ex))
    fig = plt.figure(figsize=(3.4 + 2.6 * ncol, 5.2))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.7], hspace=0.28)
    banner = outer[0].subgridspec(1, 2, width_ratios=[1.4, 1], wspace=0.15)
    draw_schematic(fig.add_subplot(banner[0, 0]), sk, typing)
    draw_matrix(fig.add_subplot(banner[0, 1]), compute_matrix(SKELS[sk], tuple(typing)),
                tuple(typing), title="ProSMoS query")
    grid = outer[1].subgridspec(1, ncol, wspace=0.05)
    for k, e in enumerate(ex):
        png = f"{sk:04d}_{ty:02d}_{e['did']}.png"
        struct_panel(fig.add_subplot(grid[0, k]), png,
                     f"{e['did']} · T {e['group']}", f"{e['count']} dom")
    tag = ("UNITYPICAL — 1 topology" if view == "unitypical"
           else f"PROMISCUOUS — {c['nT']} topologies (top {len(ex)})")
    fig.suptitle(f"skeleton {sk:04d} · typing {typing} · {c['nhit']} hits · {tag}",
                 fontsize=11, y=0.99)
    plt.subplots_adjust(left=0.04, right=0.97, top=0.9, bottom=0.06)
    fig.savefig(OUT / f"{sk:04d}_{ty:02d}.png", dpi=100)
    plt.close(fig)


def main():
    N = int(os.environ.get("NTASKS", "1"))
    T = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    cells = json.load((GAL / "full_spec.json").open())["cells"]
    sz = (len(cells) + N - 1) // N
    mine = cells[T * sz:(T + 1) * sz]
    done = 0
    for c in mine:
        out = OUT / f"{c['sk']:04d}_{c['ty']:02d}.png"
        if out.exists():
            done += 1; continue
        compose(c); done += 1
    print(f"task {T}/{N}: composed {done}/{len(mine)}")


main()
