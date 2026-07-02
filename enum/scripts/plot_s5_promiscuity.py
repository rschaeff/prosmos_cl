#!/usr/bin/env python3
"""Static promiscuity heatmap: distinct ECOD T-groups and H-groups per S5 cell."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

SP = Path("/tmp/claude-1219/-home-rschaeff-dev-prosmos-cl/02f20625-920b-46f1-a2bf-bc06d84727af/scratchpad")
D = json.load((SP / "s5_promiscuity.json").open())
NSK, NTY = D["nsk"], D["nty"]
nT = np.array(D["nT"], float)
nH = np.array(D["nH"], float)
nHits = np.array(D["nHits"], float)

for A in (nT, nH):
    A[nHits == 0] = np.nan  # unoccupied -> grey

fig = plt.figure(figsize=(13, 9))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[3, 1], hspace=0.28, wspace=0.22)

cmap = plt.cm.turbo.copy()
cmap.set_bad("#d9d9d9")
vmax = np.nanmax(nT)
norm = mcolors.LogNorm(vmin=1, vmax=vmax)

for ax, A, lab in [(fig.add_subplot(gs[0, 0]), nT, "distinct T-groups (topology)"),
                   (fig.add_subplot(gs[0, 1]), nH, "distinct H-groups (homology)")]:
    im = ax.imshow(A, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
    ax.set_title(f"# {lab} per cell", fontsize=11)
    ax.set_xlabel("H/E typing 0→31")
    ax.set_ylabel("skeleton rank (shared order)")
    cb = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cb.set_label("distinct ECOD groups (log)")

# distribution panel
axd = fig.add_subplot(gs[1, :])
flatT = nT[np.isfinite(nT)]
flatH = nH[np.isfinite(nH)]
bins = np.logspace(0, np.log10(vmax + 1), 26)
axd.hist(flatT, bins=bins, color="#d1495b", alpha=0.7, label="T-groups")
axd.hist(flatH, bins=bins, color="#4c72b0", alpha=0.55, label="H-groups")
axd.set_xscale("log")
axd.set_xlabel("distinct ECOD groups spanned by a cell's hits")
axd.set_ylabel("# cells")
uni = D["stats"]["unitypicalT"]; occ = D["stats"]["occupied"]
axd.axvline(1, color="#2ca02c", lw=1.2, ls=":")
axd.annotate(f"unitypical (nT=1): {uni} cells\n= {100*uni/occ:.0f}% of {occ} occupied",
             xy=(1, axd.get_ylim()[1] * 0.7), xytext=(1.4, axd.get_ylim()[1] * 0.7),
             fontsize=9, color="#2ca02c", va="center")
axd.legend(fontsize=9)
axd.set_title("How fold-specific is each geometric cell?  left = unitypical (one fold),  right = promiscuous (shared across folds)",
             fontsize=10)

fig.suptitle("ECOD-group promiscuity of the S5 skeleton×typing matrix "
             f"(ECOD manual reps · {occ} occupied cells)", fontsize=12.5, y=0.98)
plt.savefig(SP / "s5_promiscuity.png", dpi=155, bbox_inches="tight")
print("wrote s5_promiscuity.png")
