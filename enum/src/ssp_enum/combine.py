"""Combine-pairs growth (Chitturi 2016 Appendix §1.1 + §1.2).

Paper Appendix § "secondary compactness criterion":

    A skeleton s with more than two SSEs is formed by combining two
    smaller skeletons s1 and s2 where |s| = |s1| + |s2|. Among s1 and
    s2; s1 is larger. ... we either augment m with a single node
    (|s2| = 1) or another skeleton (1 < |s2| ≤ m) with at most m nodes.

This module currently implements the `|s2| = 1` case only — single-node
growth from a larger skeleton. That alone gives Sn = Sn-1 + S1 for any
n, missing the cases where the smaller piece has 2+ nodes (S4 = S2 + S2,
S5 = S3 + S2). Phase B extends to multi-node combine.

Paper Appendix pseudocode (single-node simplified):

    valid(s1) = {p ∉ s1 :
        p is collinear with at least two points of s1, OR
        p is adjacent to at least two points of s1}

    For each p ∈ valid(s1):
        Form skeleton s by adding p at the end of s1 (label |s1| + 1), and
        also by adding p at the front of s1 (label 1, shift s1 labels +1).
        Check PCC, SCC-1, SCC-2 on s.
        If passes, deduplicate against existing candidates via canonical form.

The front/end joining distinction matters: same lattice arrangement with
different sequence labels is a different skeleton (paper Methods: "The
labels of the grid points are changed to reflect their new positions in
the combined skeleton").
"""

from __future__ import annotations

from typing import Iterator

from .compactness import _are_collinear, is_compact
from .lattice import LatticePoint
from .skeleton import Skeleton


def adjacent_points(skel: Skeleton) -> set[LatticePoint]:
    """All lattice points adjacent to ≥1 point of `skel` but not in `skel`."""
    skel_set = set(skel.points)
    result: set[LatticePoint] = set()
    for p in skel.points:
        for n in p.neighbors():
            if n not in skel_set:
                result.add(n)
    return result


def valid_extension_points(skel: Skeleton) -> set[LatticePoint]:
    """Lattice points where a new node can be placed under SCC-1.

    Paper Appendix § "let valid(s1) ⊂ adjacent(s1) where every point
    p ∈ valid(s1) is either collinear to at least two points in s1 or
    is adjacent to at least two points in s1."

    For points adjacent to exactly one skeleton point we require
    collinearity with a pair of skeleton points (the anchor and another
    skeleton point that is adjacent to the anchor — same form as
    `passes_scc_1`).
    """
    valid: set[LatticePoint] = set()
    for candidate in adjacent_points(skel):
        adj_existing = [p for p in skel.points if candidate.is_adjacent(p)]
        if len(adj_existing) >= 2:
            valid.add(candidate)
            continue
        if len(adj_existing) == 1:
            anchor = adj_existing[0]
            for q in skel.points:
                if q is anchor:
                    continue
                if anchor.is_adjacent(q) and _are_collinear(q, anchor, candidate):
                    valid.add(candidate)
                    break
    return valid


def combine_with_single_node(s1: Skeleton) -> Iterator[Skeleton]:
    """Yield candidate skeletons formed by adding a single node to `s1`.

    For each valid extension point, produces two candidates:
    one where the new node is appended (end-join, label |s1|+1), and
    one where it is prepended (front-join, label 1, old labels shifted +1).

    The produced candidates have `chirality=None`; chirality assignment
    is a separate downstream step. The candidates are NOT filtered or
    deduped here — the caller is responsible for `is_compact()` and
    canonical-form deduplication.
    """
    for new_point in valid_extension_points(s1):
        yield Skeleton(points=s1.points + (new_point,))
        yield Skeleton(points=(new_point,) + s1.points)


# Hex axial-coord lattice symmetries: 6 rotations × 2 reflections = 12 ops.
# Applied to (q, r) only; z is independent.

def _rotate60(p: LatticePoint) -> LatticePoint:
    return LatticePoint(-p.r, p.q + p.r, p.z)


def _reflect_q(p: LatticePoint) -> LatticePoint:
    return LatticePoint(p.q + p.r, -p.r, p.z)


def _z_flip(p: LatticePoint) -> LatticePoint:
    return LatticePoint(p.q, p.r, -p.z)


def canonical_key(skel: Skeleton) -> tuple:
    """Canonical-form key for a labeled skeleton.

    Quotients by translation, the 12 hex-XY symmetries (6 rotations × 2
    reflections), and layer-flip (z → -z). Sequence labels are NOT
    relabeled: two skeletons with the same lattice arrangement but
    different sequence numbering produce different keys (CG-2012 keeps
    these distinct).

    Layer-flip is included because at this stage we are enumerating
    *labeled skeletons* without chirality assignment — chirality (which
    would distinguish layer-flipped pairs) is layered on top later. For
    purposes of matching against CG-2012's per-skeleton record count,
    layer-flip-equivalent skeletons should collapse.
    """
    pts = skel.points

    def translate_to_origin(seq: tuple[LatticePoint, ...]) -> tuple[LatticePoint, ...]:
        o = seq[0]
        return tuple(LatticePoint(p.q - o.q, p.r - o.r, p.z) for p in seq)

    def key_of(seq: tuple[LatticePoint, ...]) -> tuple:
        return tuple((p.q, p.r, p.z) for p in seq)

    candidates: list[tuple] = []
    for z_op in (lambda x: x, _z_flip):
        seq0 = tuple(z_op(p) for p in pts)
        cur = translate_to_origin(seq0)
        for _ in range(6):
            candidates.append(key_of(cur))
            reflected = tuple(_reflect_q(p) for p in cur)
            candidates.append(key_of(translate_to_origin(reflected)))
            cur = tuple(_rotate60(p) for p in cur)
            cur = translate_to_origin(cur)

    return min(candidates)
