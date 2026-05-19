"""Combine-pairs growth (Chitturi 2016 Appendix §1.1 + §1.2).

Paper Appendix § "secondary compactness criterion":

    A skeleton s with more than two SSEs is formed by combining two
    smaller skeletons s1 and s2 where |s| = |s1| + |s2|. Among s1 and
    s2; s1 is larger. ... we either augment m with a single node
    (|s2| = 1) or another skeleton (1 < |s2| ≤ m) with at most m nodes.

This module implements both cases:

  * `combine_with_single_node(s1)` — |s2| = 1. For each valid extension
    point of s1, yields both front- and end-join skeletons.
  * `combine_two_skeletons(s1, s2)` — |s2| ≥ 2. For each anchor point of
    s2, each valid extension point of s1, and each 60° rotation of s2
    around the anchor, positions s2 and yields both front-/end-join
    combined skeletons (overlap with s1 rejects).

Paper Appendix pseudocode:

    valid(s1) = {p ∉ s1 :
        p is collinear with at least two points of s1, OR
        p is adjacent to at least two points of s1}

    L1) For every skeleton s2 of size q, each point of s2 is superposed
        on each point of valid(s1).
    L2) For each such superposition, s2 is rotated around the chosen
        point in xy by 60k degrees where 0 ≤ k ≤ 5.

The front/end joining distinction matters: same lattice arrangement
with different sequence labels is a different skeleton (paper Methods:
"The labels of the grid points are changed to reflect their new
positions in the combined skeleton").

Deferred to Phase C:
  - SSE orientation tracking. Paper Appendix §: "if the orientations of
    the start node of s2 and the end node of s1 are the same then s2 is
    rotated 180 degrees around X-axis." We don't track up/down node
    orientation yet, so the conflict-resolution rotation isn't applied
    here. Skeletons with mismatched orientation are still emitted; they
    will be filtered (or have s2 z-flipped) once orientation lands.
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


def _position_s2(
    s2: Skeleton,
    anchor: LatticePoint,
    target: LatticePoint,
    rotation_k: int,
) -> tuple[LatticePoint, ...]:
    """Place s2 so its `anchor` lands at `target`, rotated by 60°·k around the anchor.

    Anchor is one of s2's own points (the "selected point" in the
    paper's wording). The transformation: translate s2 so anchor is at
    origin, rotate by k·60° about origin, translate so anchor sits at
    target. Z coordinates pass through unchanged (hex rotation is in
    xy only).
    """
    out: list[LatticePoint] = []
    for p in s2.points:
        dq = p.q - anchor.q
        dr = p.r - anchor.r
        for _ in range(rotation_k % 6):
            dq, dr = -dr, dq + dr
        out.append(LatticePoint(dq + target.q, dr + target.r, p.z))
    return tuple(out)


def combine_two_skeletons(s1: Skeleton, s2: Skeleton) -> Iterator[Skeleton]:
    """Multi-node combine of `s1` (larger) with `s2` (1 < |s2| ≤ |s1|).

    For each anchor in s2 × each ext_pt in valid(s1) × each 60° rotation k,
    position s2 with anchor at ext_pt, rotated by k. Reject if any
    positioned s2 point coincides with an s1 point (paper Appendix
    pseudocode L2: "if none of the points of s1 and s2 intersect ...").
    Yield both front-join (combined = positioned_s2 + s1.points) and
    end-join (combined = s1.points + positioned_s2) candidates.

    The caller is responsible for `is_compact()` and canonical-form
    deduplication. For |s2| = 1, delegate to `combine_with_single_node`
    (which is geometrically the same but avoids the rotation loop, since
    a single-point skeleton is rotation-invariant).
    """
    if s2.dim == 1:
        yield from combine_with_single_node(s1)
        return

    s1_set = set(s1.points)
    extension_pts = valid_extension_points(s1)
    for ext_pt in extension_pts:
        for anchor in s2.points:
            for k in range(6):
                positioned = _position_s2(s2, anchor, ext_pt, k)
                pos_set = set(positioned)
                if len(pos_set) < len(positioned):
                    continue   # rotation collapsed two s2 points (shouldn't on hex, but guard)
                if pos_set & s1_set:
                    continue   # overlap with s1: rejected
                yield Skeleton(points=positioned + s1.points)        # front-join
                yield Skeleton(points=s1.points + positioned)        # end-join


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
