"""Assign SSE types and turn Skeletons into ProSMoS-ready SSPRecords.

`skeletons_to_records` walks the Cartesian product of `('H','E')` over each
skeleton's nodes and yields one `SSPRecord` per (skeleton, type-assignment)
combo. Each record carries:

- the interaction-matrix derived from skeleton lattice adjacency, the SSE
  types, the sequence-direction alternation (`Skeleton.orientations`), and
  the sheet partition (connected components of E-E lattice adjacency);
- `same_sheet` directives — one per sheet of >= 2 strands;
- `diff_sheet` directives — one per pair of E SSEs in *different* sheets;
- `handedness` directives — one per (i, j, k) triple with a non-coplanar
  sign in `handedness_signature(skel)`.

Codes follow `prosmos.py` Methods (paper §1.1) plus the lower-case
'u'/'v' helix-helix and non-H-bond conventions seen in CG-2012 IA.txt:

    EE adj, same dir, same sheet -> 'c'   (parallel H-bond)
    EE adj, opp dir,  same sheet -> 't'   (antiparallel H-bond)
    EE adj, same dir, diff sheet -> 'u'   (parallel non-H-bond)
    EE adj, opp dir,  diff sheet -> 'v'   (antiparallel non-H-bond)
    HE adj, same dir              -> 'C'   (parallel helix-strand contact)
    HE adj, opp dir               -> 'T'   (antiparallel helix-strand)
    HH adj, same dir              -> 'u'
    HH adj, opp dir               -> 'v'
    EE non-adj, same sheet        -> '-'   (long-range strand pair within sheet)
    everything else non-adj       -> 'X'   (optional)
"""

from __future__ import annotations

from itertools import combinations, product
from typing import Iterator

from .oracle import SSPRecord
from .skeleton import Skeleton
from .combine import handedness_signature


def _interaction_code(
    ti: str, tj: str, *, adjacent: bool, same_dir: bool, same_sheet: bool,
) -> str:
    if not adjacent:
        if ti == "E" and tj == "E" and same_sheet:
            return "-"
        return "X"
    # Lattice-adjacent below.
    if ti == "E" and tj == "E":
        if same_sheet:
            return "c" if same_dir else "t"
        return "u" if same_dir else "v"
    if ti == "H" and tj == "H":
        return "u" if same_dir else "v"
    # Mixed H/E.
    return "C" if same_dir else "T"


def _find_sheets(
    types: tuple[str, ...], adj: list[list[bool]],
) -> list[list[int]]:
    """Connected components of E-E lattice adjacency. 1-based, sorted."""
    n = len(types)
    seen = [False] * n
    sheets: list[list[int]] = []
    for start in range(n):
        if seen[start] or types[start] != "E":
            continue
        comp: list[int] = []
        stack = [start]
        while stack:
            v = stack.pop()
            if seen[v]:
                continue
            seen[v] = True
            comp.append(v + 1)
            for w in range(n):
                if not seen[w] and types[w] == "E" and adj[v][w]:
                    stack.append(w)
        if comp:
            sheets.append(sorted(comp))
    return sheets


def _build_record(
    skel: Skeleton,
    skel_id: int,
    type_idx: int,
    types: tuple[str, ...],
    adj: list[list[bool]],
    orientations: tuple[bool, ...],
    hsig: tuple[int, ...],
    triples: list[tuple[int, int, int]],
) -> SSPRecord:
    n = len(types)
    sheets = _find_sheets(types, adj)
    sheet_of = {idx: k for k, s in enumerate(sheets) for idx in s}

    # Upper-triangular matrix.
    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            if j == i:
                row.append("*")
            elif j < i:
                row.append("")  # lower triangle empty (parser convention)
            else:
                ti, tj = types[i], types[j]
                same_sheet = (
                    ti == "E" and tj == "E"
                    and (i + 1) in sheet_of
                    and sheet_of[i + 1] == sheet_of.get(j + 1, -1)
                )
                row.append(_interaction_code(
                    ti, tj,
                    adjacent=adj[i][j],
                    same_dir=(orientations[i] == orientations[j]),
                    same_sheet=same_sheet,
                ))
        matrix.append(tuple(row))

    same_sheet_dirs = tuple(tuple(s) for s in sheets if len(s) >= 2)

    diff_sheet_dirs: list[tuple[int, int]] = []
    e_indices = [k + 1 for k, t in enumerate(types) if t == "E"]
    for ii in range(len(e_indices)):
        for jj in range(ii + 1, len(e_indices)):
            a = e_indices[ii]; b = e_indices[jj]
            sa = sheet_of.get(a); sb = sheet_of.get(b)
            if sa is not None and sb is not None and sa != sb:
                diff_sheet_dirs.append((a, b))

    hand: list[tuple[int, int, int, str]] = []
    for (i, j, k), sign in zip(triples, hsig):
        if sign > 0:
            hand.append((i, j, k, "R"))
        elif sign < 0:
            hand.append((i, j, k, "L"))

    return SSPRecord(
        dim=n,
        skeleton_id=skel_id,
        third_idx=0,
        sub_first=type_idx,
        sub_second=0,
        sse_types=types,
        matrix=tuple(matrix),
        same_sheet=same_sheet_dirs,
        diff_sheet=tuple(diff_sheet_dirs),
        handedness=tuple(hand),
    )


def skeletons_to_records(
    skeletons: list[Skeleton],
    *,
    sse_alphabet: tuple[str, ...] = ("H", "E"),
) -> Iterator[SSPRecord]:
    """For each skeleton, yield one SSPRecord per type assignment.

    With the default `sse_alphabet=('H','E')` and skeleton dim n, emits
    `2**n` records per skeleton. Records are deterministic in
    (skeleton position in `skeletons`, lexicographic type index).
    """
    for skel_id, skel in enumerate(skeletons):
        n = skel.dim
        adj_tuple = skel.adjacency_matrix()
        adj = [[bool(c) for c in row] for row in adj_tuple]
        orientations = skel.orientations
        hsig = handedness_signature(skel)
        triples = list(combinations(range(1, n + 1), 3))

        for type_idx, types in enumerate(product(sse_alphabet, repeat=n)):
            yield _build_record(
                skel, skel_id, type_idx, types,
                adj, orientations, hsig, triples,
            )
