"""Schematic 2x3 figure showing how an S5 skeleton + an H/E typing combine
into a ProSMoS interaction matrix. Top row: the RLM-canonical hit-rich
skeleton s5-0098. Bottom row: a zero-hit negspace skeleton s5-0001.

  col 1: hex-lattice topology (5 nodes, edges from adjacency, numbered 1-5)
  col 2: same with EHEHE typing applied (H=blue, E=red)
  col 3: resulting ProSMoS interaction matrix (per-pair code letters)

Code legend for col 3 (from assignment.py:_interaction_code):
  -  diagonal (self) or two paired-strand cells inside a sheet
  X  non-adjacent SSEs (no interaction)
  u  H-H or E-E (out-of-sheet) packed parallel
  v  H-H or E-E (out-of-sheet) packed antiparallel
  c  paired strands inside a sheet, parallel
  t  paired strands inside a sheet, antiparallel
  C  mixed H-E adjacent, same direction
  T  mixed H-E adjacent, opposite direction
  X  unconstrained wildcard — matches any relationship (NOT a non-contact requirement)
"""

import math
import sys
from itertools import product
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "enum" / "src"))

from ssp_enum.assignment import _find_sheets, _interaction_code  # noqa: E402
from ssp_enum.enumerate import enumerate_skeletons  # noqa: E402

OUT = REPO / "enum" / "docs" / "figures" / "s5_skeleton_schematic.png"

H_COLOR = "#3a73c4"  # blue for helix
E_COLOR = "#c44e3a"  # red-orange for strand
NODE_OUTLINE = "#222222"
EDGE_COLOR = "#666666"


def hex_to_xy(q: int, r: int) -> tuple[float, float]:
    """Axial hex coords -> Cartesian for plotting (pointy-top)."""
    x = q + 0.5 * r
    y = -r * (math.sqrt(3) / 2)  # flip y for natural orientation
    return x, y


def draw_skeleton(ax, skel, types=None, title="", node_size=900):
    pts = [hex_to_xy(p.q, p.r) for p in skel.points]
    adj = skel.adjacency_matrix()
    n = skel.dim

    # Edges from upper-triangular adjacency
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i][j]:
                xi, yi = pts[i]
                xj, yj = pts[j]
                ax.plot([xi, xj], [yi, yj], "-", color=EDGE_COLOR, lw=2.5, zorder=1)

    # Nodes
    for idx, (x, y) in enumerate(pts):
        if types is None:
            color = "#dddddd"
        else:
            color = H_COLOR if types[idx] == "H" else E_COLOR
        ax.scatter([x], [y], s=node_size, c=color, edgecolors=NODE_OUTLINE,
                   linewidths=2, zorder=3)
        ax.text(x, y, str(idx + 1), ha="center", va="center",
                fontsize=14, fontweight="bold", zorder=4,
                color="white" if types is not None else "black")

    # Background hex-lattice hint (faint dots at nearby lattice points)
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    xpad, ypad = 1.0, 1.0
    for q in range(-4, 5):
        for r in range(-4, 5):
            x, y = hex_to_xy(q, r)
            if (min(xs) - xpad <= x <= max(xs) + xpad) and \
               (min(ys) - ypad <= y <= max(ys) + ypad):
                ax.scatter([x], [y], s=8, c="#cccccc", zorder=0)

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=11)


def compute_matrix(skel, types):
    n = skel.dim
    adj = [[bool(skel.adjacency_matrix()[i][j] or skel.adjacency_matrix()[j][i])
            for j in range(n)] for i in range(n)]
    sheets = _find_sheets(types, adj)
    matrix = [["" for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                matrix[i][j] = "-"
                continue
            same_dir = (skel.orientations[i] == skel.orientations[j])
            same_sheet = any((i + 1 in sh) and (j + 1 in sh) for sh in sheets)
            matrix[i][j] = _interaction_code(
                types[i], types[j],
                adjacent=adj[i][j], same_dir=same_dir, same_sheet=same_sheet,
            )
    return matrix


CODE_COLORS = {
    "-": "#ffffff",  # diagonal / paired-sheet
    "X": "#f4f4f4",  # non-adjacent (faint grey)
    "u": "#a8d8a8",  # H-H/E-E parallel
    "v": "#88c088",  # H-H/E-E antiparallel
    "c": "#f4b942",  # sheet-paired parallel
    "t": "#e88a2a",  # sheet-paired antiparallel
    "C": "#c089d6",  # mixed H-E parallel
    "T": "#9a5fc4",  # mixed H-E antiparallel
}


def draw_matrix(ax, matrix, types, title=""):
    n = len(matrix)
    for i in range(n):
        for j in range(n):
            code = matrix[i][j]
            color = CODE_COLORS.get(code, "#dddddd")
            ax.add_patch(plt.Rectangle((j, n - 1 - i), 1, 1, facecolor=color,
                                       edgecolor="#333333", lw=0.5))
            ax.text(j + 0.5, n - 1 - i + 0.5, code, ha="center", va="center",
                    fontsize=12, fontweight="bold")
    # Row/col labels showing type
    for i in range(n):
        t = types[i]
        c = H_COLOR if t == "H" else E_COLOR
        ax.text(-0.4, n - 1 - i + 0.5, f"{i+1}{t}", ha="right", va="center",
                fontsize=11, fontweight="bold", color=c)
        ax.text(i + 0.5, n + 0.15, f"{i+1}{t}", ha="center", va="bottom",
                fontsize=11, fontweight="bold", color=c)

    ax.set_xlim(-0.6, n + 0.2)
    ax.set_ylim(-0.2, n + 0.7)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=11)


def main() -> None:
    skels = enumerate_skeletons(5)

    # Use EHEHE typing (typing index 21 — Rossmann-canonical) for both rows
    # so the type-vector is held constant and ONLY the skeleton differs.
    ehehe = ("E", "H", "E", "H", "E")

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))

    # --- Top row: hit-rich RLM-canonical skeleton s5-0098 ---
    s = skels[98]
    draw_skeleton(axes[0, 0], s, types=None,
                  title=f"s5-0098 (hit-rich, RLM-canonical)\nhex-lattice skeleton: 5 SSE positions, edges = adjacency")
    draw_skeleton(axes[0, 1], s, types=ehehe,
                  title=f"+ EHEHE typing (H = blue, E = red)")
    m = compute_matrix(s, ehehe)
    draw_matrix(axes[0, 2], m, ehehe,
                title=f"ProSMoS query matrix (one of 2^5 = 32 typings)\n"
                      f"hits in v4 manual-rep DB: {2910}")

    # --- Bottom row: zero-hit negspace skeleton s5-0001 ---
    s2 = skels[1]
    draw_skeleton(axes[1, 0], s2, types=None,
                  title=f"s5-0001 (one of 25 zero-hit skeletons)\nsame lattice model, different topology")
    draw_skeleton(axes[1, 1], s2, types=ehehe,
                  title=f"+ EHEHE typing (same as above)")
    m2 = compute_matrix(s2, ehehe)
    draw_matrix(axes[1, 2], m2, ehehe,
                title=f"Different ProSMoS query → different result\n"
                      f"hits in v4 manual-rep DB: 0 (and zero for all 32 typings)")

    # Shared legend for matrix code colors
    legend_handles = []
    legend_labels = {
        "-": "diagonal / sheet-pair (-)",
        "X": "unconstrained / any relationship (X)",
        "u": "α-α or β-β (out-of-sheet) parallel (u)",
        "v": "α-α or β-β antiparallel (v)",
        "c": "β-β paired, parallel sheet (c)",
        "t": "β-β paired, antiparallel sheet (t)",
        "C": "α-β adjacent, parallel (C)",
        "T": "α-β adjacent, antiparallel (T)",
    }
    for code, lbl in legend_labels.items():
        legend_handles.append(mpatches.Patch(facecolor=CODE_COLORS[code],
                                             edgecolor="#333", label=lbl))
    fig.legend(handles=legend_handles, loc="lower center", ncol=4,
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(
        "From a lattice skeleton to a ProSMoS query — and why some are zero-hit",
        fontsize=13, y=0.99,
    )
    plt.tight_layout(rect=(0, 0.04, 1, 0.97))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUT, dpi=180, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
