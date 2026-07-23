#!/usr/bin/env python3
"""Render a Ruczinski-style four-strand (or n-strand) topology cartoon as SVG,
from the motif encoding -- so we can draw the FULL enumeration, not only the 48
panels that happened to be traced as PNGs in 2002.

ENCODING (the notation settled during the 48-panel transcription):
  strand_order : the spatial order of SEQUENCE strands, left-to-right as drawn.
                 strand_order[p] = which sequence-strand occupies spatial slot p.
                 "1423" (the pretzel) = slot0 holds seq1, slot1 seq4, slot2 seq2,
                 slot3 seq3.
  ori_drawn    : up/down per SPATIAL slot, left-to-right. 'u' = arrowhead up
                 (N-terminus at the bottom, C at the top); 'd' = arrowhead down
                 (N at top, C at bottom).

DRAWING. Each strand is an outlined block arrow at its spatial slot. The backbone
is then traced N->C: for each sequence bond i -> i+1 we connect the C-terminus of
seq-strand i to the N-terminus of seq-strand i+1 with a loop that bulges away from
the sheet on the side the C-terminus sits. N and C of the whole chain are labelled.

Strokes use currentColor and fills are 'none' / the page background, so the same
SVG reads correctly on light or dark (the app is dark-themed; the 2002 PNGs are
black-on-white).
"""
from __future__ import annotations

# geometry (viewBox units)
PITCH = 34        # horizontal spacing between strands
BODY_W = 12       # arrow body width
HEAD_W = 20       # arrowhead width
HEAD_H = 12       # arrowhead height
Y_TOP = 30        # top of a strand's vertical span
Y_BOT = 92        # bottom
MARGIN_X = 22
LOOP = 15         # how far loops bulge past the strand ends


def _arrow(cx: float, up: bool) -> str:
    """Outlined block arrow centred at x=cx, spanning Y_TOP..Y_BOT.
    up=True -> head at top (C-terminus up)."""
    hw, bw = HEAD_W / 2, BODY_W / 2
    if up:
        # head at Y_TOP, body down to Y_BOT
        pts = [
            (cx - bw, Y_BOT), (cx - bw, Y_TOP + HEAD_H), (cx - hw, Y_TOP + HEAD_H),
            (cx, Y_TOP), (cx + hw, Y_TOP + HEAD_H), (cx + bw, Y_TOP + HEAD_H),
            (cx + bw, Y_BOT),
        ]
    else:
        pts = [
            (cx - bw, Y_TOP), (cx - bw, Y_BOT - HEAD_H), (cx - hw, Y_BOT - HEAD_H),
            (cx, Y_BOT), (cx + hw, Y_BOT - HEAD_H), (cx + bw, Y_BOT - HEAD_H),
            (cx + bw, Y_TOP),
        ]
    d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts) + " Z"
    return f'<path d="{d}" fill="none" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/>'


def _termini(up: bool):
    """(N_y, C_y) for a strand. up: N at bottom, C at top."""
    return (Y_BOT, Y_TOP) if up else (Y_TOP, Y_BOT)


def render(strand_order: str, ori_drawn: str, *, label: str | None = None,
           size: int = 120) -> str:
    n = len(strand_order)
    assert len(ori_drawn) == n
    ups = [c == "u" for c in ori_drawn]
    cx = [MARGIN_X + p * PITCH for p in range(n)]
    # sequence-strand s (1-based) -> spatial slot
    slot = {int(strand_order[p]): p for p in range(n)}

    parts = []
    # backbone loops first (so arrows draw on top of the join)
    for s in range(1, n):
        pi, pj = slot[s], slot[s + 1]
        xi, xj = cx[pi], cx[pj]
        _, cyi = _termini(ups[pi])          # C-terminus of strand s
        nyj, _ = _termini(ups[pj])          # N-terminus of strand s+1
        c_top = cyi == Y_TOP
        n_top = nyj == Y_TOP
        out_i = (Y_TOP - LOOP) if c_top else (Y_BOT + LOOP)
        out_j = (Y_TOP - LOOP) if n_top else (Y_BOT + LOOP)
        if c_top == n_top:
            # same side: single smooth arc bulging out over the midpoint
            mx = (xi + xj) / 2
            my = out_i
            d = f"M{xi:.1f},{cyi:.1f} Q{mx:.1f},{my:.1f} {xj:.1f},{nyj:.1f}"
        else:
            # opposite sides (parallel-neighbour crossover): a rounded S that
            # leaves C straight out past the strand end, sweeps across in the gap,
            # and enters N straight in -- vertical control handles keep it from
            # slashing diagonally through the sheet.
            k = 22
            ci = cyi - k if c_top else cyi + k
            cj = nyj - k if n_top else nyj + k
            d = (f"M{xi:.1f},{cyi:.1f} C{xi:.1f},{ci:.1f} "
                 f"{xj:.1f},{cj:.1f} {xj:.1f},{nyj:.1f}")
        parts.append(f'<path d="{d}" fill="none" stroke="currentColor" '
                     f'stroke-width="1.6"/>')

    for p in range(n):
        parts.append(_arrow(cx[p], ups[p]))

    # N / C labels at the chain termini
    n_slot = slot[1]
    ny, _ = _termini(ups[n_slot])
    parts.append(_lbl(cx[n_slot], ny, "N"))
    c_slot = slot[n]
    _, cy = _termini(ups[c_slot])
    parts.append(_lbl(cx[c_slot], cy, "C"))

    w = MARGIN_X * 2 + (n - 1) * PITCH
    h = 122
    head = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" '
            f'width="{size}" role="img" '
            f'aria-label="four-strand topology {strand_order} {ori_drawn}">')
    lab = (f'<text x="4" y="14" font-size="12" fill="currentColor" '
           f'font-family="sans-serif">{label}</text>') if label else ""
    return head + lab + "".join(parts) + "</svg>"


def _lbl(cx: float, y: float, txt: str) -> str:
    dy = 13 if y == Y_BOT else -5
    return (f'<text x="{cx:.1f}" y="{y + dy:.1f}" font-size="11" '
            f'text-anchor="middle" fill="currentColor" '
            f'font-family="sans-serif">{txt}</text>')


if __name__ == "__main__":
    import sys
    from pathlib import Path
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    out.mkdir(parents=True, exist_ok=True)
    # smoke set: the pretzel + panel 1, for eyeball vs the PNGs
    for so, od, lab in [("1423", "uddu", "23"), ("1234", "uudu", "1"),
                        ("1243", "uuuu", "2"), ("1342", "uuud", "28")]:
        (out / f"motif_{lab}.svg").write_text(render(so, od, label=lab))
        print(f"wrote motif_{lab}.svg  ({so} {od})")
