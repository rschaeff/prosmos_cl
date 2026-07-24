#!/usr/bin/env python3
"""Extend a typed ProSMoS motif query by one SSE, preserving the seed exactly.

WHY THIS EXISTS (and why the lattice route does not work here). The obvious way
to reach S6 near a Ruczinski panel is via its S5 cell. That route is a dead end,
but NOT for the reasons first claimed. Re-derived against the canonical TYPED
graph198 grid (`ruczinski/panel_parent_geom.py`, 2026-07-24):
  * NOT orientation-blindness. graph198 carries c/t (parallel/antiparallel), so
    the earlier "the grid can't see up/down" argument was an artifact of the
    orientation-blind `adjonly` grid -- it does not apply to the canonical grid.
  * NOT "the grid can't hold the sheet." It can: panels 3 and 10 each embed in
    10 typed S5 cells (all EEEEE) -- verified by exhaustive code-matrix match.
  * The real reason is SHEET CARDINALITY. Those 10 embeddings are all LOOSE: the
    panel's four strands appear only as a SUB-sheet of a 5-strand (EEEEE) sheet,
    never as a complete 4-strand sheet with the 5th SSE outside it. STRICT S5
    parents (panel = a complete S5 sheet) number 0. Since ProSMoS sheet-
    completeness only matches a panel as a COMPLETE sheet, no S5 cell is a parent.
    (Panel 23 is the adjacent case: its best embedding is sk171 ty07 = HHEEE, only
    3 of its 4 strands as a complete sheet + 2 helices -- also not a parent.)

So extensions are grown in the RUCZINSKI ENCODING instead: the seed's own `c`/`t`
codes are carried through untouched, and one SSE is added. Output is again a
valid typed query, so applying the tool twice gives S6.

WHAT IS ENUMERATED. The new SSE is specified MINIMALLY -- one declared contact,
everything else `X` (optional). That is deliberate: a minimal query is a superset
query, so a helix declared against strand 2 still matches structures where it
also packs on strands 1 and 3. Over-specifying would silently exclude them.

  strand_edge  a new strand H-bonded to a sheet EDGE strand, parallel (`c`) or
               antiparallel (`t`), joining that sheet. Edges are derived from the
               matrix: a strand with exactly one c/t partner. Interior strands
               have two and cannot accept another neighbour.
  helix_pack   a new helix packing on one strand, parallel (`C`) or
               antiparallel (`T`), in no sheet.

Each geometric option is crossed with every sequence insertion slot (before the
first SSE, between any two, after the last), because ProSMoS indexes by sequence
order and a motif's N->C placement is part of its identity.

READ THE RESULTS CONDITIONALLY. Extensions of one seed form a sibling set sharing
a parent and a sampling depth. Every added SSE multiplies the constraint count,
so hit counts fall mechanically and some extension is always empty. Ask whether
THIS addition kills realizability when its siblings survive -- never quote a bare
"0 hits => dark".
"""
from __future__ import annotations

import argparse
from pathlib import Path

HBOND = {"c", "t"}          # strand-strand, H-bonded (i.e. sheet-adjacent)
PACK = {"C", "T"}           # helix-strand packing
MIN_LEN = {"E": 5, "H": 8}
MAX_LEN = 1000


class Query:
    """A parsed ProSMoS query. Matrix is a full n x n grid; only i<j is used."""

    def __init__(self, text: str):
        lines = [l.rstrip("\n") for l in text.splitlines() if l.strip()]
        self.n = len(lines[0].split())
        self.types = lines[1].split()
        self.matrix = [["X"] * self.n for _ in range(self.n)]
        for i in range(self.n):
            for j, c in enumerate(lines[2 + i].split()):
                self.matrix[i][i + j] = c
        self.sheets: list[list[int]] = []
        self.handedness: list[tuple[int, int, int, str]] = []
        for l in lines[2 + self.n:]:
            f = l.split()
            if f[0] == "sheetS":
                self.sheets.append([int(x) for x in f[1:]])
            elif f[0] == "handedness":
                self.handedness.append((int(f[1]), int(f[2]), int(f[3]), f[4]))

    def partners(self, i: int) -> list[int]:
        """0-based indices H-bonded to strand i (its sheet neighbours)."""
        out = []
        for j in range(self.n):
            if j == i:
                continue
            a, b = (i, j) if i < j else (j, i)
            if self.matrix[a][b] in HBOND:
                out.append(j)
        return out

    def edge_strands(self) -> list[int]:
        """0-based strands with exactly one H-bond partner: the sheet's edges."""
        return [i for i in range(self.n)
                if self.types[i] == "E" and len(self.partners(i)) == 1]

    def render(self) -> str:
        out = [" ".join(str(i + 1) for i in range(self.n)), " ".join(self.types)]
        for i in range(self.n):
            cells = ["*" if j == i else self.matrix[i][j] for j in range(i, self.n)]
            out.append(" " * (2 * i) + " ".join(cells))
        for sh in self.sheets:
            out.append("sheetS " + " ".join(str(x) for x in sh))
        for i, t in enumerate(self.types, start=1):
            out.append(f"length {i} {t} {MIN_LEN.get(t, 5)} {MAX_LEN}")
        for a, b, c, lr in self.handedness:
            out.append(f"handedness {a} {b} {c} {lr}")
        return "\n".join(out) + "\n"


def insert(q: Query, slot: int, new_type: str, anchor: int, code: str,
           join_sheet: bool) -> Query:
    """Insert `new_type` at sequence `slot` (0..n), contacting old SSE `anchor`
    (0-based, pre-insertion) with `code`. Everything else is 'X'."""
    new = Query.__new__(Query)
    n = q.n + 1
    new.n = n
    # old 0-based index -> new 0-based index
    shift = lambda i: i + 1 if i >= slot else i          # noqa: E731
    new.types = q.types[:slot] + [new_type] + q.types[slot:]
    new.matrix = [["X"] * n for _ in range(n)]
    for i in range(q.n):                                  # carry the seed intact
        for j in range(i + 1, q.n):
            a, b = shift(i), shift(j)
            lo, hi = (a, b) if a < b else (b, a)
            new.matrix[lo][hi] = q.matrix[i][j]
    a, b = slot, shift(anchor)
    lo, hi = (a, b) if a < b else (b, a)
    new.matrix[lo][hi] = code
    new.sheets = [[shift(x - 1) + 1 for x in sh] for sh in q.sheets]
    if join_sheet:
        for sh in new.sheets:
            if shift(anchor) + 1 in sh:
                sh.append(slot + 1)
                sh.sort()
                break
    new.handedness = [(shift(x - 1) + 1, shift(y - 1) + 1, shift(z - 1) + 1, lr)
                      for x, y, z, lr in q.handedness]
    return new


def extensions(q: Query):
    """Yield (query, kind, slot, anchor, code) for every minimal extension."""
    for slot in range(q.n + 1):
        for anchor in q.edge_strands():                   # new strand at a sheet edge
            for code in ("c", "t"):
                yield (insert(q, slot, "E", anchor, code, True),
                       "strand_edge", slot, anchor, code)
        for anchor in range(q.n):                         # new helix packing on a strand
            if q.types[anchor] != "E":
                continue
            for code in ("C", "T"):
                yield (insert(q, slot, "H", anchor, code, False),
                       "helix_pack", slot, anchor, code)


def contains_seed(ext: Query, seed: Query, slot: int) -> bool:
    """The seed must survive cell-for-cell, and none of its sheets may split."""
    shift = lambda i: i + 1 if i >= slot else i           # noqa: E731
    for i in range(seed.n):
        for j in range(i + 1, seed.n):
            a, b = shift(i), shift(j)
            lo, hi = (a, b) if a < b else (b, a)
            if ext.matrix[lo][hi] != seed.matrix[i][j]:
                return False
        if ext.types[shift(i)] != seed.types[i]:
            return False
    for sh in seed.sheets:
        want = {shift(x - 1) + 1 for x in sh}
        if not any(want <= set(h) for h in ext.sheets):
            return False
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("seed", type=Path, help="a typed ProSMoS .query file")
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--tag", default=None, help="filename prefix (default: seed stem)")
    a = ap.parse_args()

    seed = Query(a.seed.read_text())
    tag = a.tag or a.seed.stem
    a.out.mkdir(parents=True, exist_ok=True)
    rows, n_emit, rejected = [], 0, 0
    seen = set()
    for ext, kind, slot, anchor, code in extensions(seed):
        if not contains_seed(ext, seed, slot):
            rejected += 1
            continue
        key = ext.render()
        if key in seen:
            continue
        seen.add(key)
        name = f"{tag}-x{n_emit:03d}"
        (a.out / f"{name}.query").write_text(key)
        rows.append((name, tag, "".join(seed.types), "".join(ext.types),
                     kind, slot + 1, anchor + 1, code))
        n_emit += 1
    mf = a.out / "manifest.tsv"
    header = ("query\tseed\tseed_types\text_types\tkind\tnew_sse_pos\t"
              "anchor_sse\tcontact_code\n")
    with open(mf, "a" if mf.exists() else "w") as fh:
        if not mf.exists() or mf.stat().st_size == 0:
            fh.write(header)
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"{a.seed.name}: seed {''.join(seed.types)} "
          f"edges={[i+1 for i in seed.edge_strands()]} -> {n_emit} extensions "
          f"({rejected} rejected)")


if __name__ == "__main__":
    main()
