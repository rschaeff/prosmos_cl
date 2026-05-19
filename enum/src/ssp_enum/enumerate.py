"""Dimension-by-dimension SSP enumeration.

Paper Methods, "In our model of SSP generation, we combine two smaller
skeletons to obtain a larger one." S_n is grown from pairs (S_p, S_q) with
p + q = n; deduplication is via lattice symmetry and mirror-image pairing.

Currently implemented:
  - S1 (a single lattice point)
  - S2 (a seed pair: two adjacent points)
  - S3 *planar* (Z=0 only): the 4 base spatial-sequence arrangements that,
    when crossed with handedness via Z-displacement, yield CG-2012's
    11 S3 SSPs. Handedness extension is the natural next step.

Out of scope here:
  - General S_n for n >= 4 (needs full combine + canonical-form)
  - SSE-type assignment (H/E) and interaction-type assignment
  - SCC-2 (Appendix Fig. S2d) — kicks in for n >= 5
"""

from __future__ import annotations

from typing import Iterator

from .compactness import is_compact
from .lattice import LatticePoint
from .skeleton import Skeleton


def _canonical_planar_s3(skel: Skeleton) -> tuple:
    """Canonical-form key for a planar S3 skeleton.

    Two skeletons sharing this key represent the same spatial-sequence
    arrangement up to translation, hex-rotation, and reflection of the
    XY plane. The order of `points` is *preserved* (sequence labeling
    matters), so this is the smallest hashable representation that
    deduplicates over pure lattice symmetry.

    Strategy: translate so point[0] is at origin, then enumerate the
    12 hex isometries (6 rotations × 2 reflections), pick the
    lex-smallest tuple-of-tuples of coords.
    """
    pts = skel.points

    def translate(p, dq, dr):
        return LatticePoint(p.q + dq, p.r + dr, p.z)

    # Translate so pts[0] sits at origin
    o = pts[0]
    centered = tuple(translate(p, -o.q, -o.r) for p in pts)

    # 6-fold rotations + reflections in axial coords:
    # Rotation by 60°: (q, r) -> (-r, q+r)
    # Reflection across the q-axis: (q, r) -> (q+r, -r)
    def rotate60(p):
        return LatticePoint(-p.r, p.q + p.r, p.z)

    def reflect(p):
        return LatticePoint(p.q + p.r, -p.r, p.z)

    candidates = []
    cur = centered
    for _ in range(6):
        candidates.append(cur)
        candidates.append(tuple(reflect(p) for p in cur))
        cur = tuple(rotate60(p) for p in cur)

    def key(seq):
        return tuple((p.q, p.r, p.z) for p in seq)

    return min(key(c) for c in candidates)


def enumerate_s3_planar() -> Iterator[Skeleton]:
    """Yield the unique planar S3 skeletons (Z=0 throughout).

    Sequence position 1 is anchored at the origin (lattice translation
    symmetry). Positions 2 and 3 may be at *any* lattice point in a
    bounded region — sequence-consecutive SSEs need not be spatially
    adjacent (loops can be long), so `p1`–`p2` is not required to be
    a lattice edge. Compactness (`compactness.is_compact`) filters the
    result; deduplication is by S3 adjacency matrix, the canonical-form
    invariant for this dimension under lattice symmetry.

    For S3, the four distinct adjacency patterns are:
        Linear  (T,F,T) — p1-p2-p3 path with p1,p3 NOT adjacent
        Bent-A  (F,T,T) — p3 in spatial middle (jump from p1 to p2)
        Bent-B  (T,T,F) — p1 in spatial middle (jump from p2 to p3)
        Tri     (T,T,T) — all three lattice-adjacent

    These cross-product with handedness via Z-displacement to produce
    CG-2012's 11 S3 SSPs.
    """
    p1 = LatticePoint(0, 0, 0)
    box = [
        LatticePoint(q, r, 0)
        for q in range(-2, 3)
        for r in range(-2, 3)
        if (q, r) != (0, 0)
    ]
    seen: set[tuple] = set()
    for p2 in box:
        for p3 in box:
            if p3 == p2:
                continue
            skel = Skeleton(points=(p1, p2, p3))
            if not is_compact(skel):
                continue
            m = skel.adjacency_matrix()
            # For dim 3, the upper-triangle adjacency triple is a complete
            # canonical-form invariant under lattice symmetry: any two
            # planar S3 skeletons with the same (1,2), (1,3), (2,3) edge
            # pattern are equivalent under hex rotation/reflection.
            key = (m[0][1], m[0][2], m[1][2])
            if key in seen:
                continue
            seen.add(key)
            yield skel


def enumerate_dim(n: int) -> Iterator[Skeleton]:
    """Top-level: dispatch to dimension-specific generators where implemented."""
    if n == 1:
        yield Skeleton(points=(LatticePoint(0, 0, 0),))
        return
    if n == 2:
        yield Skeleton(points=(LatticePoint(0, 0, 0), LatticePoint(1, 0, 0)))
        return
    if n == 3:
        yield from enumerate_s3_planar()
        return
    raise NotImplementedError(f"enumerate_dim({n}); pending")
