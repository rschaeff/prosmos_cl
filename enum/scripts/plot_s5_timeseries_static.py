#!/usr/bin/env python3
"""Static verification / paper figures from s5_timeseries.json:
  fig1  the answer: occupancy (new cells) vs cumulative hits over deposition year
  fig2  first-appearance-year heatmap + 4 cumulative snapshot frames
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

SP = Path("/tmp/claude-1219/-home-rschaeff-dev-prosmos-cl/02f20625-920b-46f1-a2bf-bc06d84727af/scratchpad")
D = json.load((SP / "s5_timeseries.json").open())
YEARS = D["years"]; NY = len(YEARS); NSK = D["nsk"]; NTY = D["nty"]

# reconstruct cumulative cube from sparse series ------------------------------
inc = np.zeros((NY, NSK, NTY))
for r, ty, yi, c in D["series"]:
    inc[yi, r, ty] += c
cum = np.cumsum(inc, axis=0)              # cum[yi] = counts <= YEARS[yi], rows in display order
first = np.array(D["firstGrid"]).reshape(NSK, NTY).astype(float)
first[first < 0] = np.nan                 # never-occupied -> grey

occ = np.array(D["occByYear"]); tot = np.array(D["totByYear"])

# ---- fig1: the answer curve -------------------------------------------------
# marginal discovery rate: new cells per 1000 new hits, in each year
inc_hits = np.diff(tot, prepend=0)
inc_cells = np.diff(occ, prepend=0)
with np.errstate(divide="ignore", invalid="ignore"):
    marg = np.where(inc_hits > 0, inc_cells / inc_hits * 1000, np.nan)

fig, (ax, axm) = plt.subplots(2, 1, figsize=(9, 7), height_ratios=[2.4, 1],
                              sharex=True, gridspec_kw=dict(hspace=0.12))
ax2 = ax.twinx()
l1, = ax.plot(YEARS, occ, color="#2ca02c", lw=2.6, label="distinct occupied cells (breadth)")
l2, = ax2.plot(YEARS, tot / 1000, color="#555", lw=2.2, ls="--",
               label="cumulative hits, thousands (depth)")
ax.set_ylabel("occupied skeleton×typing cells", color="#2ca02c")
ax2.set_ylabel("cumulative hits (thousands)", color="#555")
ax.tick_params(axis="y", colors="#2ca02c")
ax.set_ylim(0, NSK * NTY * 0.42)
ax.axhline(occ[-1], color="#2ca02c", lw=0.6, ls=":", alpha=0.6)
sat = YEARS[int(np.argmax(occ >= 0.9 * occ[-1]))]
ax.axvline(sat, color="#999", lw=0.8, ls=":")
axm.axvline(sat, color="#999", lw=0.8, ls=":")
ax.annotate(f"90% of ever-occupied\ncells reached by {sat}",
            xy=(sat, 0.9 * occ[-1]), xytext=(sat - 20, occ[-1] * 0.62),
            fontsize=9, color="#2ca02c",
            arrowprops=dict(arrowstyle="->", color="#2ca02c", lw=1))
ax.legend(handles=[l1, l2], loc="lower right", fontsize=9, framealpha=0.9)
ax.set_title("Did the PDB change the sampled S5 matrix, or fill it in?\n"
             "The reachable region is mapped by the mid-2000s; later structures densify it → FILLED IN",
             fontsize=11)
# marginal panel
axm.plot(YEARS, marg, color="#8b4fc9", lw=2)
axm.fill_between(YEARS, marg, color="#8b4fc9", alpha=0.15)
axm.set_ylabel("new cells per\n1000 new hits", fontsize=9)
axm.set_xlabel("PDB deposition year")
axm.set_ylim(0, np.nanmax(marg[3:]) * 1.05)
axm.annotate("each new structure is\nless and less likely to\nopen a new cell",
             xy=(2010, marg[YEARS.index(2010)]), xytext=(1996, np.nanmax(marg[3:]) * 0.6),
             fontsize=8.5, color="#6a3ba0",
             arrowprops=dict(arrowstyle="->", color="#6a3ba0", lw=0.9))
plt.savefig(SP / "fig1_fill_vs_change.png", dpi=170, bbox_inches="tight")
plt.close()

# ---- fig2: first-appearance + snapshots ------------------------------------
frames = [1990, 2000, 2010, YEARS[-1]]
fig = plt.figure(figsize=(15, 8))
gs = fig.add_gridspec(1, 5, width_ratios=[1.25, 1, 1, 1, 1], wspace=0.25)

vmax = cum[-1].max()
norm = mcolors.SymLogNorm(linthresh=1, vmin=0, vmax=vmax)

axf = fig.add_subplot(gs[0, 0])
im = axf.imshow(first, aspect="auto", cmap="plasma",
                vmin=0, vmax=NY - 1, interpolation="nearest")
axf.set_title("First-appearance year\n(when each cell first got a hit)", fontsize=10)
axf.set_xlabel("typing 0→31"); axf.set_ylabel("skeleton rank (fixed order)")
axf.set_facecolor("#cccccc")
cb = fig.colorbar(im, ax=axf, shrink=0.5, ticks=[0, NY // 2, NY - 1])
cb.ax.set_yticklabels([YEARS[0], YEARS[NY // 2], YEARS[-1]])

for j, yr in enumerate(frames):
    yi = YEARS.index(yr)
    axs = fig.add_subplot(gs[0, j + 1])
    axs.imshow(cum[yi], aspect="auto", cmap="viridis", norm=norm,
               interpolation="nearest")
    n_occ = int((cum[yi] > 0).sum())
    axs.set_title(f"≤ {yr}\n{n_occ} cells, {int(cum[yi].sum()):,} hits", fontsize=10)
    axs.set_xlabel("typing 0→31")
    axs.set_yticks([])
fig.suptitle("S5 matrix (198 skeletons × 32 typings) sampled by experimental PDB over time — "
             "cells light up early, then just deepen", fontsize=12, y=0.99)
plt.savefig(SP / "fig2_frames.png", dpi=150, bbox_inches="tight")
plt.close()
print("wrote fig1_fill_vs_change.png, fig2_frames.png")
