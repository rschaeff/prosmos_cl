#!/usr/bin/env python3
"""Compose the full 96 four-strand topology cartoons into one grid figure.

Reconstructed 2026-07-24 (the original all-96 driver was a lost one-off, like the
pdb_exp mkseg step). Reads the motif annotation table and renders each motif with
render_motif_svg.render(), tiling 96 cells into a single SVG. The 48 Ruczinski
Fig-1 "forbidden" motifs are highlighted.

ENCODING PER MOTIF. render() draws (strand_order, ori_drawn) in the AS-DRAWN
convention. The 48 forbidden motifs carry the as-drawn panel values
(panel_strand_order, panel_ori_drawn) transcribed from Fig 1 — used so those cells
reproduce the 2002 panels exactly. The other 48 (never panels) have no as-drawn
transcription, so they fall back to the canonical (spatial_seq, orientation); each
is a valid depiction of the same topology.

Annotations: work/prosmos_2026/ruczinski/queries/motifs.tsv (96 rows).
Output: enum/docs/figures/ruczinski_all96_topology.svg  (currentColor; theme-safe).
PNG (optional raster preview): ImageMagick can't resolve currentColor, so
  sed 's/currentColor/black/g' out.svg | convert - out.png
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from render_motif_svg import render  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
TSV = Path.home() / "work/prosmos_2026/ruczinski/queries/motifs.tsv"
OUT = REPO / "docs/figures/ruczinski_all96_topology.svg"

# grid + cell geometry
COLS, ROWS = 8, 12                       # 8 x 12 = 96
MOTIF_W, MOTIF_H = 146, 122              # render() viewBox for a 4-strand motif
CW = 170                                 # cell width
PAD = 10
IN_W = CW - 2 * PAD                       # nested motif width
IN_H = round(IN_W * MOTIF_H / MOTIF_W)    # keep aspect
CAP = 22                                  # caption strip
CH = IN_H + CAP + PAD                      # cell height
H_TINT = "#ffab21"                        # forbidden-panel highlight (theme-safe accent)

_svg_wrap = re.compile(r'^<svg\b[^>]*?>(.*)</svg>\s*$', re.DOTALL)


def place(motif_svg: str, x: float, y: float) -> str:
    """Embed a render() output as a translated+scaled <g> group. A <g> transform is
    far more portable across SVG rasterizers than a nested <svg> viewport (which
    ImageMagick mis-renders), and the app inlines it the same way."""
    inner = _svg_wrap.match(motif_svg).group(1)
    s = IN_W / MOTIF_W
    return f'<g transform="translate({x:.1f},{y:.1f}) scale({s:.4f})">{inner}</g>'


def main():
    rows = {int(r["motif_id"]): r for r in
            csv.DictReader(TSV.open(), delimiter="\t")}
    assert len(rows) == 96, f"expected 96 motifs, got {len(rows)}"

    W, H = COLS * CW, ROWS * CH + 34
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" '
        f'width="{W}" role="img" aria-label="all 96 four-strand topologies">',
        f'<text x="{PAD}" y="22" font-size="16" font-weight="bold" '
        f'fill="currentColor" font-family="sans-serif">'
        f'All 96 four-strand &#946;-sheet topologies '
        f'(the 48 Ruczinski 2002 &#8220;forbidden&#8221; panels highlighted)</text>',
    ]
    Y0 = 34
    for mid in range(1, 97):
        r = rows[mid]
        col, row = (mid - 1) % COLS, (mid - 1) // COLS
        cx, cy = col * CW, Y0 + row * CH
        forbidden = r["ruczinski_forbidden"] == "yes"
        # as-drawn for the 48 panels; canonical for the rest
        if forbidden and r["panel_ori_drawn"]:
            so, od = r["panel_strand_order"], r["panel_ori_drawn"]
        else:
            so, od = r["spatial_seq"], r["orientation"]
        if forbidden:                    # highlight tint behind the cell
            parts.append(
                f'<rect x="{cx+3:.0f}" y="{cy+3:.0f}" width="{CW-6}" height="{CH-6}" '
                f'rx="8" fill="{H_TINT}" fill-opacity="0.12" '
                f'stroke="{H_TINT}" stroke-opacity="0.55" stroke-width="1.5"/>')
        parts.append(place(render(so, od), cx + PAD, cy + PAD))
        panel = f' &#183; panel {r["fig1_panel"]}' if forbidden else ""
        parts.append(
            f'<text x="{cx + CW/2:.0f}" y="{cy + IN_H + PAD + 14:.0f}" '
            f'font-size="11" text-anchor="middle" fill="currentColor" '
            f'font-family="sans-serif">m{mid}{panel}</text>')
    parts.append("</svg>")
    OUT.write_text("".join(parts))
    n_forb = sum(1 for r in rows.values() if r["ruczinski_forbidden"] == "yes")
    print(f"wrote {OUT} ({len(rows)} motifs, {n_forb} highlighted)")


if __name__ == "__main__":
    main()
