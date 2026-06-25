"""Plot a 198 × 32 heatmap of S5 (skeleton × typing) hit counts from the
v4 manual-reps sweep, with the 25 fully-zero-hit skeletons grouped at the
bottom to visually demonstrate the negspace = 25 rows × 32 typings = 800
queries.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


REPO = Path(__file__).resolve().parents[2]
SUMMARY = REPO / "results" / "v4_ecod_manual_reps_summary.tsv"
OUT = REPO / "enum" / "docs" / "figures" / "s5_hit_grid.png"


def load_s5_grid(path: Path) -> np.ndarray:
    """Return a (198, 32) array of hit counts indexed by (skeleton_id, typing_id)."""
    grid = np.zeros((198, 32), dtype=int)
    with path.open() as f:
        next(f)  # header
        for line in f:
            name, _runtime, _rc, hits = line.rstrip().split("\t")
            if not name.startswith("s5-"):
                continue
            _, sk, ty = name.split("-")
            grid[int(sk), int(ty)] = int(hits)
    return grid


def main() -> None:
    grid = load_s5_grid(SUMMARY)

    # Sort skeletons by total hit count (descending). This puts the 25 fully-
    # zero-hit skeletons together at the bottom, where they form a visually
    # obvious block.
    totals = grid.sum(axis=1)
    order = np.argsort(-totals)  # descending
    sorted_grid = grid[order]
    sorted_ids = order  # original skeleton ids in new row order

    n_zero = int((totals == 0).sum())  # = 25
    assert n_zero == 25, f"expected 25 fully-zero-hit S5 skeletons, got {n_zero}"

    # Log-scale color (0 stays white; lit-up cells ramp through viridis).
    # Use SymLogNorm so 0 is exact white and we can still see 1 vs 10 vs 100.
    vmax = sorted_grid.max()
    norm = mcolors.SymLogNorm(linthresh=1, vmin=0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(9, 12))
    im = ax.imshow(
        sorted_grid,
        aspect="auto",
        cmap="viridis",
        norm=norm,
        interpolation="nearest",
    )

    # Outline the negspace block: bottom 25 rows × all 32 columns.
    neg_start = 198 - n_zero - 0.5  # row boundary just above the zero block
    ax.axhline(neg_start, color="red", linewidth=1.5, linestyle="-")
    ax.annotate(
        f"25 zero-hit skeletons × 32 typings\n= 800 negspace queries",
        xy=(31.5, neg_start + (n_zero / 2)),
        xytext=(34, neg_start + (n_zero / 2)),
        fontsize=11,
        color="red",
        va="center",
        ha="left",
        annotation_clip=False,
        arrowprops=dict(arrowstyle="->", color="red", lw=1.2),
    )

    ax.set_xlabel("H/E typing index (0 = HHHHH, 31 = EEEEE)")
    ax.set_ylabel("S5 skeleton (sorted by total hit count, descending)")
    ax.set_title(
        "S5 sweep against v4 manual-rep DB (19,015 experimental PDB domains)\n"
        "Each cell = one of 198 × 32 = 6,336 queries; color = hit count (symlog)"
    )

    # Mark composition tick groups (H-count) on the x-axis as a hint.
    # typing_idx bit count == strand count; group by composition.
    e_counts = [bin(t).count("1") for t in range(32)]
    # show transitions between compositions
    transitions = [i for i in range(1, 32) if e_counts[i] != e_counts[i - 1]]
    for t in transitions:
        ax.axvline(t - 0.5, color="white", linewidth=0.3, alpha=0.4)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("ECOD manual-rep hits (symlog scale; 0 = absent)")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(OUT, dpi=180, bbox_inches="tight")
    print(f"wrote {OUT}")
    print(f"  total queries:    {grid.size}")
    print(f"  zero-hit queries: {(grid == 0).sum()}")
    print(f"  fully-zero skeletons: {n_zero} -> {n_zero * 32} queries forming the negspace block")


if __name__ == "__main__":
    main()
