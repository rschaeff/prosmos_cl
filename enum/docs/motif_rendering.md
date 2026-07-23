# Motif topology rendering — status 2026-07-23

Generated SVG topology cartoons for the enumerated motifs, so the **full**
enumeration can be depicted (not just the 48 S4 panels Ruczinski traced in 2002).
Two renderers, both driven from the encoding — nothing is traced from the PNGs.

## What exists

### `enum/scripts/render_motif_svg.py` — 4-strand Ruczinski cartoons (WORKS)
`render(strand_order, ori_drawn)` → SVG. Outlined up/down arrows in spatial
order, N→C backbone loops (over-top / under-bottom, with a rounded S for
parallel-neighbour crossovers), N/C labels. Validated against the originals:
m15/p23 (the pretzel) and m56/p1 reproduce their panels.
- Full set: `enum/docs/figures/ruczinski_all96_topology.{svg,png}` — all 96
  four-strand motifs, the 48 Fig-1 panels highlighted.
- Faithful and clean. This one is done.

### `enum/scripts/render_s5_iso.py` — isometric S*n* with helices (USABLE, not polished)
`render_cell(skeleton, typing)` → SVG. The 4-strand side-view does NOT generalise:
S5 skeletons are 2-D hex arrangements (only 12/198 are planar; most span 2–3
rows), so a side view collapses the depth and a top-down view collapses up/down.
Uses an **orthographic isometric** projection instead — hex (q,r) as a tilted
floor, each SSE a vertical rod:
- strand = arrow (arrowhead = C-terminus, up/down = N→C direction),
- helix = cylinder with a directional C-terminal cone + direction-slanted coils,
- backbone = arcs in the linker's actual plane (top when the SSE points up,
  bottom when down — orientations alternate, so each linker sits cleanly in one
  plane), painter's-ordered back-to-front,
- N/C at the true chain termini.
- H/E palette `#5a8fd6` / `#d2703c` passes the dataviz validator in BOTH light and
  dark; shape (cylinder vs arrow) carries H/E independent of colour. Strokes use
  currentColor.
- Examples: `enum/docs/figures/s5_iso_examples.png` (sk171 ty07 HHEEE — the
  AFDB-only survivor — plus EHEHH / HHHHH / EEHEE).

## The limitation (why this is parked)

Flat SVG primitives hit a ceiling on **depth**: linkers and rods behind a front
cylinder get occluded, so densely-stacked skeletons (e.g. HHHHH, sk171) read less
cleanly than the linear ones (sk98 EHEHH is crisp). The current renderer draws all
backbone behind all rods; a per-segment depth interleave (draw each arc at its
rearmost endpoint's depth, between rod layers) would help but does not remove the
fundamental problem — two rods projecting to nearly the same screen-x still crowd.

Getting arbitrary S5/S6 packing to read well from 2-D primitives is genuinely
hard. If we revisit, the options worth exploring are NOT more primitive-tuning:
- a real 2.5-D/3-D scene (three.js / a ribbon library) with lighting + true
  occlusion, rendered once to static SVG/PNG per cell;
- a molecular cartoon (PyMOL/Mol\* from an idealised backbone built from the
  lattice) — heavier, but photoreal and unambiguous;
- a TOPS-style abstraction that deliberately drops the 3-D and shows adjacency as
  a graph (loses the spatial intuition the isometric was meant to give).

## Reproduce
```
python enum/scripts/render_motif_svg.py <outdir>     # 4-strand smoke set
python enum/scripts/render_s5_iso.py   <outdir>      # S5 iso smoke set
```
Both emit currentColor SVGs (~1–2 KB); ImageMagick can't resolve currentColor, so
for a raster preview substitute a colour first: `sed 's/currentColor/black/g'`.
