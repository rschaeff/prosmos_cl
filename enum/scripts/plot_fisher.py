#!/usr/bin/env python3
"""4-panel: the two occupancy heatmaps + Fisher significance + depth-corrected."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.colors as mcolors

G = "/home/rschaeff/work/prosmos_2026/s5_grid"
A = np.load(f"{G}/grid_afdb_nT_rebuilt.npy").astype(float)
P = np.load(f"{G}/grid_pdb_nT.npy").astype(float)
IMP = np.load(f"{G}/impossible_mask.npy")
Z = np.load(f"{G}/fisher_cells.npz")
NA, NP = 3656, 3935
cells = Z["cells"]

sigF = np.full(A.shape, np.nan)
sigD = np.full(A.shape, np.nan)
for (sk, ty), qf, qd in zip(cells, Z["q"], Z["q_depth"]):
    # signed -log10(q): + = commoner in AFDB, - = commoner in PDB (per-fold rate)
    rateA, rateP = A[sk, ty] / NA, P[sk, ty] / NP
    s = 1.0 if rateA >= rateP else -1.0
    sigF[sk, ty] = s * -np.log10(max(qf, 1e-12))
    sigD[sk, ty] = s * -np.log10(max(qd, 1e-12))

order = np.argsort(-(A / NA).sum(1))
def srt(M): return M[order]
imp = IMP[order]

fig, ax = plt.subplots(1, 4, figsize=(21, 11))
norm = mcolors.SymLogNorm(linthresh=1 / NA, vmin=0, vmax=max((A/NA).max(), (P/NP).max()))
for a, M, t in [(ax[0], srt(A / NA), f"AFDB — fold rate\n{NA:,} folds searched"),
                (ax[1], srt(P / NP), f"PDB — fold rate\n{NP:,} folds searched")]:
    a.imshow(M, aspect="auto", cmap="viridis", norm=norm, interpolation="nearest")
    a.set_title(t, fontsize=11)
    a.set_xlabel("H/E typing (0=HHHHH -> 31=EEEEE)", fontsize=9)

lim = 6
cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("#7a7a7a")
for a, M, t in [(ax[2], np.ma.masked_invalid(srt(sigF)), "FISHER on fold counts\nsigned -log10(q), BH"),
                (ax[3], np.ma.masked_invalid(srt(sigD)), "DEPTH-CORRECTED\nAFDB thinned to PDB depth")]:
    M = np.ma.masked_where(imp, M)
    im = a.imshow(-M, aspect="auto", cmap=cmap, vmin=-lim, vmax=lim, interpolation="nearest")
    a.set_title(t, fontsize=11)
    a.set_xlabel("H/E typing", fontsize=9)
for a in ax[1:]:
    a.set_yticklabels([])
ax[0].set_ylabel("S5 skeleton (sorted by AFDB fold rate)", fontsize=10)
cb = fig.colorbar(im, ax=ax[2:], orientation="horizontal", fraction=0.035, pad=0.07)
cb.set_label("signed -log10(q):  blue = commoner in AFDB   ·   red = commoner in PDB   ·   grey = no test", fontsize=9)

nF = int((Z["q"] < 0.05).sum()); nD = int((Z["q_depth"] < 0.05).sum())
fig.suptitle("Per-topology PDB vs AFDB significance — Fisher, and what survives a depth correction", fontsize=14, y=0.965)
fig.text(0.5, 0.055,
         f"Fisher calls {nF:,} of {len(cells):,} tested cells significant (q<0.05); after thinning AFDB to PDB's structure depth only {nD:,} remain.\n"
         "In the zero-cells: PDB=0 & AFDB>0 -> Fisher 42 significant, depth-corrected 5. AFDB=0 & PDB>0 -> 0 either way (the test cannot\n"
         "reach significance in that direction at these counts). Fold-level counts throughout: structure-level would treat ~30x within-fold redundancy as evidence.",
         ha="center", fontsize=10, color="#933")
out = f"{G}/fisher_panels.png"
plt.savefig(out, dpi=130, bbox_inches="tight"); plt.close()
print("wrote", out)
