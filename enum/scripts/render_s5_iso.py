#!/usr/bin/env python3
"""Isometric SVG renderer for S5 (and any hex-lattice S*n*) motifs, with helix
and strand glyphs.

WHY ISOMETRIC. In the Chitturi model each SSE is a unit rod PERPENDICULAR to the
hex xy-plane, sandwiched between z=0 and z=1, standing at a hex (q,r) position.
A top-down view collapses the up/down direction; a side view (the Ruczinski
4-strand cartoon) collapses the hex depth. Only an isometric projection shows
all three at once: the hex floor tilted into 2.5-D, each SSE a vertical rod, and
its N->C direction as which end the arrow/coil points.

GLYPHS:
  strand (E) : a flat vertical arrow rod. up -> arrowhead at top (C at z=1),
               down -> arrowhead at bottom (C at z=0).
  helix  (H) : a vertical cylinder with a couple of coil hatches.
Orientation alternates from the skeleton's start_up (paper Appendix 1.1), so it
comes straight off `skeleton.orientations`.

BACKBONE. Consecutive SSEs are joined N->C by a linker that lies in the z=1 plane
(if the C-strand end is up) or the z=0 plane (if down) -- again per the paper.
Drawn as a line on that plane between the two floor positions.

Painter's algorithm orders rods back-to-front by projected depth so nearer rods
occlude farther ones.

Strokes use currentColor; H/E are distinguished by SHAPE first (cylinder vs
arrow) and tinted with two theme-safe hues as reinforcement, so the encoding
never rests on color alone.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "enum" / "src"))

from ssp_enum.enumerate import enumerate_skeletons  # noqa: E402

# ---- isometric projection parameters
XSCALE = 34.0        # hex q,r horizontal spread on screen
DEPTH = 30.0         # how much hex depth (r) pushes up-and-back (the squash)
SSE_H = 40.0         # rod height (z=0 -> z=1)
ROD_W = 9.0          # half-width of a strand rod / helix radius
HEAD = 13.0          # arrowhead height

H_TINT = "#5a8fd6"   # helix (blue-ish); shape already disambiguates
E_TINT = "#d2703c"   # strand (orange-ish)


def hex_flat(q: int, r: int) -> tuple[float, float]:
    """Axial hex -> flat plane coords (pointy-top), before iso tilt."""
    return (q + 0.5 * r, r * (math.sqrt(3) / 2))


def floor_xy(hx: float, hy: float) -> tuple[float, float]:
    """Flat plane point -> screen coords of its z=0 (floor) projection.
    hy (depth) pushes the point UP the screen (further back) and is squashed."""
    return (hx * XSCALE, -hy * DEPTH)


def _strand(sx: float, z0y: float, up: bool, tint: str) -> str:
    top, bot = z0y - SSE_H, z0y
    w = ROD_W
    if up:   # arrowhead at top
        body = f'<rect x="{sx-w:.1f}" y="{top+HEAD:.1f}" width="{2*w:.1f}" height="{SSE_H-HEAD:.1f}"'
        head = (f'<path d="M{sx-w-3:.1f},{top+HEAD:.1f} L{sx:.1f},{top:.1f} '
                f'L{sx+w+3:.1f},{top+HEAD:.1f} Z" fill="{tint}" stroke="currentColor" stroke-width="1.5"/>')
    else:    # arrowhead at bottom
        body = f'<rect x="{sx-w:.1f}" y="{top:.1f}" width="{2*w:.1f}" height="{SSE_H-HEAD:.1f}"'
        head = (f'<path d="M{sx-w-3:.1f},{bot-HEAD:.1f} L{sx:.1f},{bot:.1f} '
                f'L{sx+w+3:.1f},{bot-HEAD:.1f} Z" fill="{tint}" stroke="currentColor" stroke-width="1.5"/>')
    return (f'{body} fill="{tint}" stroke="currentColor" stroke-width="1.5"/>' + head)


def _helix(sx: float, z0y: float, up: bool, tint: str) -> str:
    """Cylinder with a DIRECTIONAL C-terminal cap. Coil hatches slant with the
    N->C sense (up: rising left->right, down: falling), so a helix reads its
    direction the way a strand's arrowhead does -- the fix for 'cylinders look
    the same either way'."""
    top, bot = z0y - SSE_H, z0y
    w, ry = ROD_W, 4.0
    cap = 9.0                                  # directional cone at the C end
    parts = [
        f'<rect x="{sx-w:.1f}" y="{top:.1f}" width="{2*w:.1f}" height="{SSE_H:.1f}" '
        f'fill="{tint}" stroke="currentColor" stroke-width="1.5"/>',
        f'<ellipse cx="{sx:.1f}" cy="{bot:.1f}" rx="{w:.1f}" ry="{ry:.1f}" '
        f'fill="{tint}" stroke="currentColor" stroke-width="1.5"/>',
        f'<ellipse cx="{sx:.1f}" cy="{top:.1f}" rx="{w:.1f}" ry="{ry:.1f}" '
        f'fill="{tint}" stroke="currentColor" stroke-width="1.5"/>',
    ]
    # coil hatches, slanted by direction
    for f in (0.28, 0.5, 0.72):
        y = top + SSE_H * f
        dy = -3 if up else 3                   # rise (up) vs fall (down)
        parts.append(f'<path d="M{sx-w:.1f},{y-dy:.1f} Q{sx:.1f},{y+dy:.1f} '
                     f'{sx+w:.1f},{y-dy:.1f}" fill="none" stroke="currentColor" '
                     f'stroke-width="1"/>')
    # directional cone at the C terminus (top if up, bottom if down)
    cy = top if up else bot
    tipy = cy - cap if up else cy + cap
    parts.append(f'<path d="M{sx-w-3:.1f},{cy:.1f} L{sx:.1f},{tipy:.1f} '
                 f'L{sx+w+3:.1f},{cy:.1f} Z" fill="{tint}" stroke="currentColor" '
                 f'stroke-width="1.5" stroke-linejoin="round"/>')
    return "".join(parts)


def render_cell(skel, typing: str, *, label: str | None = None, size: int = 150) -> str:
    n = skel.dim
    ori = skel.orientations                      # tuple[bool] up/down, N->C
    floors = []
    for p in skel.points:
        hx, hy = hex_flat(p.q, p.r)
        floors.append(floor_xy(hx, hy))

    # terminus screen points per SSE: (N_pt, C_pt) in its up/down plane
    def termini(i):
        sx, z0y = floors[i]
        z1y = z0y - SSE_H
        return ((sx, z0y), (sx, z1y)) if ori[i] else ((sx, z1y), (sx, z0y))

    body = []
    # backbone: C_i -> N_{i+1} as an ARC in the linker's plane. Orientations
    # alternate, so C_i and N_{i+1} always share a plane: top (z=1) when SSE i
    # points up, bottom (z=0) when it points down. Top-plane linkers bulge
    # up-and-back, bottom-plane down-and-front -- lifting them off the rods so the
    # N->C order reads and the depth is unambiguous.
    ARC = 16.0
    for i in range(n - 1):
        (_, ci) = termini(i)
        (nj, _) = termini(i + 1)
        mx, my = (ci[0] + nj[0]) / 2, (ci[1] + nj[1]) / 2
        my += -ARC if ori[i] else ARC          # up -> top plane -> bulge up
        body.append(f'<path d="M{ci[0]:.1f},{ci[1]:.1f} Q{mx:.1f},{my:.1f} '
                    f'{nj[0]:.1f},{nj[1]:.1f}" fill="none" stroke="currentColor" '
                    f'stroke-width="1.6" stroke-linecap="round"/>')

    # rods back-to-front: larger -hy (further back) drawn first == smaller screen y-floor
    order = sorted(range(n), key=lambda i: floors[i][1])
    for i in order:
        sx, z0y = floors[i]
        tint = H_TINT if typing[i] == "H" else E_TINT
        body.append(_helix(sx, z0y, ori[i], tint) if typing[i] == "H"
                    else _strand(sx, z0y, ori[i], tint))
        # sequence number
        body.append(f'<text x="{sx:.1f}" y="{z0y - SSE_H/2 + 4:.1f}" font-size="12" '
                    f'text-anchor="middle" fill="currentColor" '
                    f'font-family="sans-serif" font-weight="bold">{i+1}</text>')

    # N / C chain-terminus labels
    (n0, _) = termini(0)
    (_, cN) = termini(n - 1)
    body.append(f'<text x="{n0[0]-11:.1f}" y="{n0[1]+4:.1f}" font-size="11" '
                f'fill="currentColor" font-family="sans-serif">N</text>')
    body.append(f'<text x="{cN[0]+5:.1f}" y="{cN[1]+4:.1f}" font-size="11" '
                f'fill="currentColor" font-family="sans-serif">C</text>')

    # bounding box
    xs = [f[0] for f in floors]; ys = [f[1] for f in floors]
    minx, maxx = min(xs) - 26, max(xs) + 26
    miny, maxy = min(ys) - SSE_H - 20, max(ys) + 22
    w, h = maxx - minx, maxy - miny
    lab = (f'<text x="{minx+4:.1f}" y="{miny+14:.1f}" font-size="12" '
           f'fill="currentColor" font-family="sans-serif">{label}</text>') if label else ""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" '
            f'viewBox="{minx:.1f} {miny:.1f} {w:.1f} {h:.1f}" width="{size}" '
            f'role="img" aria-label="S{n} motif {label or ""}">' + lab
            + "".join(body) + "</svg>")


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    out.mkdir(parents=True, exist_ok=True)
    s5 = enumerate_skeletons(5)
    # sk171 ty07 = HHEEE, the AFDB-only survivor; plus a couple of shapes
    from itertools import product
    def typing(idx):  # ty index -> H/E string (H=0,E=1, MSB=SSE1) matching the grid
        return "".join("HE"[b] for b in
                       [(idx >> (4 - k)) & 1 for k in range(5)])
    for sk, ty in [(171, 7), (98, 20), (0, 0), (158, 27)]:
        svg = render_cell(s5[sk], typing(ty), label=f"sk{sk} ty{ty:02d} {typing(ty)}")
        (out / f"s5_{sk}_{ty}.svg").write_text(svg)
        print(f"wrote s5_{sk}_{ty}.svg  ({typing(ty)})")
