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


_SQRT3_OVER_2 = 0.8660254037844386
_LAYER_EPSILON = 0.01  # matches CG-2012 LatticePoint.ZERO_THRESHOLD


def _axial_to_euclidean(p: LatticePoint) -> tuple[float, float]:
    """Convert axial hex coords to Euclidean xy."""
    return (p.q + p.r * 0.5, p.r * _SQRT3_OVER_2)


def dist_sum(skel: Skeleton) -> float:
    """Pair-wise Euclidean distance sum (RCC primary criterion).

    Port of CG-2012 CGMotif._distSum (read in checkEquivalence5Gr at
    line 0x0047 of the IL). Used as the first RCC tie-break per paper
    Appendix §1.2: "the skeleton with least pair-wise distance sum is
    preferred."
    """
    import math
    pts = [_axial_to_euclidean(p) for p in skel.points]
    n = len(pts)
    total = 0.0
    for i in range(n):
        xi, yi = pts[i]
        for j in range(i + 1, n):
            xj, yj = pts[j]
            total += math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2)
    return total


def layers(skel: Skeleton) -> int:
    """Number of unique y-coordinates in the skeleton's Euclidean projection.

    Port of CG-2012 CGMotif._layers (read in checkEquivalence5Gr at
    line 0x006e). RCC secondary tie-break per paper Appendix §1.2:
    "if it is same for both then the skeleton with minimum number of
    unique y coordinates (layers) is chosen." Uses the same
    `LatticePoint.ZERO_THRESHOLD = 0.01` epsilon as CG-2012.
    """
    ys = sorted(_axial_to_euclidean(p)[1] for p in skel.points)
    if not ys:
        return 0
    n_unique = 1
    for i in range(1, len(ys)):
        if abs(ys[i] - ys[i - 1]) > _LAYER_EPSILON:
            n_unique += 1
    return n_unique


def rcc_dedup(skeletons: list[Skeleton]) -> list[Skeleton]:
    """Among handedness-equivalent skeletons, keep the RCC winner.

    Port of CG-2012's CGMotif::checkEquivalence5Gr selection logic.
    Equivalence is by full handedness signature (the binary's `handSum1`
    is a uint64 hash of this; using the tuple directly avoids hash
    collision issues). Within each equivalence class, the canonical
    representative is the one with:

      1. lowest `dist_sum` (pair-wise distance sum), epsilon 0.01;
      2. tied → lowest `layers` (unique-y count);
      3. tied → first in input order.

    Match to paper Appendix §1.2 RCC and to the decoded CG-2012 IL.
    """
    from collections import defaultdict
    by_hand: dict[tuple, list[tuple[int, Skeleton]]] = defaultdict(list)
    for idx, s in enumerate(skeletons):
        by_hand[handedness_signature(s)].append((idx, s))
    result: list[Skeleton] = []
    for group in by_hand.values():
        if len(group) == 1:
            result.append(group[0][1])
            continue
        # Compute RCC keys once per skeleton.
        scored = [(dist_sum(s), layers(s), idx, s) for idx, s in group]
        scored.sort(key=lambda t: (t[0], t[1], t[2]))
        result.append(scored[0][3])
    return result


def handedness_signature(skel: Skeleton) -> tuple[int, ...]:
    """Per-triple handedness signature, paper Appendix §1.1.2.

    Handedness of triple (i, j, k) = sign((p_i × p_j) · p_k) where
    each p_i is the 3D node vector (Euclidean xy from axial coords,
    z = +1 if orientations[i] else -1). The orientations alternate
    along the sequence starting from `skel.start_up`.

    Returns an integer sign in {-1, 0, +1} for each triple i<j<k,
    flattened in lexicographic order of (i, j, k). Length = C(n, 3).

    Per paper: two skeletons with the same handedness signature
    (label-by-label, same `(i,j,k)` triples) are equivalent. So this
    is the natural canonical-form invariant — finer than rotation
    equivalence (rotations preserve the signature) but coarser than
    full positional identity (different lattice arrangements with the
    same signature collapse).
    """
    n = skel.dim
    orientations = skel.orientations
    vecs = []
    for i in range(n):
        x, y = _axial_to_euclidean(skel.points[i])
        z = 1.0 if orientations[i] else -1.0
        vecs.append((x, y, z))
    out: list[int] = []
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                px, py, pz = vecs[i]
                qx, qy, qz = vecs[j]
                rx, ry, rz = vecs[k]
                cx = py * qz - pz * qy
                cy = pz * qx - px * qz
                cz = px * qy - py * qx
                t = cx * rx + cy * ry + cz * rz
                if abs(t) < 1e-9:
                    out.append(0)
                elif t > 0:
                    out.append(1)
                else:
                    out.append(-1)
    return tuple(out)


def canonical_key(skel: Skeleton) -> tuple:
    """Canonical-form key for a labeled skeleton.

    Symmetry group: translation × 6 hex-XY rotations.

    Reflections and z-flip are deliberately excluded — they flip
    handedness, and the paper treats mirror-image skeletons as
    distinct (they appear as L/R chirality pairs in CG-2012 panel
    records). Two skeletons get the same key iff they are the same
    labeled lattice arrangement up to translation and 60° rotation
    about any anchor, with the same `start_up` orientation and
    `chirality`.

    Phase C3 attempt (sequence-reversal + start_up-flip) was reverted:
    sequence-reversal flips handedness signs of every labeled triple
    (paper Appendix §1.1: handedness = (p×q)·r is sign-sensitive to
    the order of arguments), so per the paper's strict label-by-label
    equivalence definition it's not a valid quotient. Empirically
    it also collapsed Bent-A/Bent-B (which are sequence-reversals
    but enumerated separately by CG-2012) and CW/CCW triangle mirrors.

    Phase C4 investigation: using `handedness_signature` as the
    canonical key (paper's stated equivalence definition) over-
    collapses relative to CG-2012's oracle. S3 has only one triple
    (max 3 distinct signatures, including the "0" collinear case),
    but oracle has 11 S3 panels. Empirically: handedness-only gives
    S4=23 / S5=193, way below oracle 41 / 648; current canonical
    gives S4=84 (over) / S5=396 (under). Conclusion: CG-2012 enumerates
    per-labeling without applying the paper's handedness-equivalence
    rule, so our rotation-only canonical is closer to oracle behavior
    than handedness-equivalence would be. `handedness_signature` is
    retained as a utility for downstream chirality assignment.

    History:
      Phase A/B: 12 hex symmetries + z-flip (too coarse — collapsed
                 mirror pairs CG-2012 keeps distinct).
      Phase C1: rotation-only (translation × 6 rotations).
      Phase C2: + start_up flag in the key (REVERTED post-C4).
      Phase C2-revert: dropped start_up from the key because pre-C2
                 S4=42 matched oracle 41 exactly while post-C2 S4=84
                 overshot by 2x. The paper says UP-start and DOWN-start
                 SSPs differ in linker-plane assignment, but CG-2012's
                 oracle counts them as one skeleton. We retain
                 `Skeleton.start_up` and `Skeleton.orientations` for
                 downstream chirality / handedness work — the field is
                 no longer part of the canonical equivalence class.
    """
    pts = skel.points

    def translate_to_origin(seq: tuple[LatticePoint, ...]) -> tuple[LatticePoint, ...]:
        o = seq[0]
        return tuple(LatticePoint(p.q - o.q, p.r - o.r, p.z) for p in seq)

    def key_of(seq: tuple[LatticePoint, ...]) -> tuple:
        return tuple((p.q, p.r, p.z) for p in seq)

    cur = translate_to_origin(pts)
    candidates: list[tuple] = []
    for _ in range(6):
        candidates.append(key_of(cur))
        cur = tuple(_rotate60(p) for p in cur)
        cur = translate_to_origin(cur)
    return (min(candidates), skel.chirality)
