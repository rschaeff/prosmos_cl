"""Dimension-by-dimension SSP enumeration.

Paper Methods, "In our model of SSP generation, we combine two smaller
skeletons to obtain a larger one." S_n is grown from pairs (S_p, S_q) with
p + q = n; deduplication is via lattice symmetry and mirror-image pairing.

Two enumeration paths live here:

  enumerate_dim(n)     - dimension dispatcher used by the existing test
                         suite. For n ∈ {1, 2, 3} it returns the same
                         results as before. (S3 still goes through the
                         hand-written `enumerate_s3_planar` + chirality
                         pipeline that yields 11 SSPs.)

  enumerate_skeletons(n) - new combine-pairs growth route. Returns
                         labeled skeletons (no chirality) by iteratively
                         growing from S_{n-1} via single-node addition.
                         Phase A scope: |s2| = 1 only. Phase B will add
                         multi-node combine (S4 = S2 + S2, S5 = S3 + S2).

Implemented:
  - S1, S2 (canonical seeds)
  - S3 *planar*: 4 base spatial-sequence arrangements (`enumerate_s3_planar`)
  - S3 *full*: 11 SSPs (`enumerate_s3`) matching CG-2012/S3/Stru.txt, with
    chirality {None, L, R} crossed with each base shape (Triangle chiral-only)
  - `enumerate_skeletons(n)`: single-node-combine growth from cached
    S_{n-1}; deduped by `combine.canonical_key`

Out of scope here:
  - Multi-node combine (S4 = S2 + S2, S5 = S3 + S2)
  - RCC tie-breaking on equivalent skeletons (paper Appendix § RCC)
  - Handedness-based equivalence (paper Appendix § "if their handedness
    is same in s1 and s2, then s1 and s2 are equivalent")
  - SSE-type assignment (H/E) and interaction-type assignment
"""

from __future__ import annotations

from typing import Iterator

from .combine import (
    canonical_key,
    combine_two_skeletons,
    combine_with_single_node,
    rcc_dedup,
)
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


def enumerate_s3() -> Iterator[Skeleton]:
    """Yield the 11 S3 SSPs (planar + chiral) matching CG-2012.

    Cross-product of the 4 base shapes from `enumerate_s3_planar` with
    chirality assignments:

      - Acyclic shapes (Linear, Bent-A, Bent-B): 3 variants each
        — unhanded (`chirality=None`), `'L'`, `'R'`.
      - Cyclic shape (Triangle, all three pairs adjacent): 2 variants
        only — `'L'`, `'R'`. There is no unhanded triangle in CG-2012;
        a closed 3-cycle is intrinsically chiral via the orientation of
        the sequence walk around the loop.

    Total: 3 × 3 + 2 = 11, matching CG-2012/S3/Stru.txt (panels 3-0..3-10).

    Panel mapping (CG-2012 panel → (shape, chirality) emitted here):
      3-0  Linear   None       3-5  Linear   L     3-6  Linear   R
      3-1  Triangle L          3-2  Triangle R
      3-3  Bent-A   None       3-7  Bent-A   L     3-8  Bent-A   R
      3-4  Bent-B   None       3-9  Bent-B   R     3-10 Bent-B   L

    Chirality is recorded as a label rather than via an explicit Z
    coordinate on lattice points; see `skeleton.py` for the rationale.
    """
    for skel in enumerate_s3_planar():
        m = skel.adjacency_matrix()
        is_cycle = m[0][1] and m[0][2] and m[1][2]  # Triangle: all three edges
        if not is_cycle:
            yield skel  # unhanded
        yield Skeleton(points=skel.points, chirality="L")
        yield Skeleton(points=skel.points, chirality="R")


_SEED_S1 = Skeleton(points=(LatticePoint(0, 0, 0),))
_SEED_S2 = Skeleton(points=(LatticePoint(0, 0, 0), LatticePoint(1, 0, 0)))


def enumerate_skeletons(n: int) -> list[Skeleton]:
    """Labeled skeletons at dimension `n` via combine-pairs growth.

    Iterates all (p, q) splits with p + q = n and p ≥ q ≥ 1. For each
    split, every (s1 ∈ Sp, s2 ∈ Sq) pair is combined (single-node when
    q == 1, multi-node otherwise). The union is deduplicated by
    `canonical_key` and filtered by `is_compact` (PCC ∧ SCC-1 ∧ SCC-2).

    Phase B (current) covers single-node and multi-node geometric combine.
    Deferred to Phase C:
      - RCC tie-breaking on equivalent skeletons (paper Appendix §1.2)
      - Handedness-based equivalence (paper Appendix §1.1: triple-wise
        handedness identity → equivalence)
      - SSE orientation tracking + the conflict-resolution 180° X-axis
        rotation when join orientations clash
    """
    if n == 1:
        return [Skeleton(points=(LatticePoint(0, 0, 0),))]
    if n == 2:
        return [Skeleton(points=(LatticePoint(0, 0, 0), LatticePoint(1, 0, 0)))]
    seen: dict[tuple, Skeleton] = {}
    for q in range(1, n // 2 + 1):
        p = n - q
        sp = enumerate_skeletons(p)
        if q == 1:
            for s1 in sp:
                # Each s1 already carries a start_up value (Sp is enumerated
                # with both variants), so single-node combine inherits it.
                for candidate in combine_with_single_node(s1):
                    candidate = Skeleton(points=candidate.points,
                                         chirality=candidate.chirality,
                                         start_up=s1.start_up)
                    if not is_compact(candidate):
                        continue
                    key = canonical_key(candidate)
                    if key not in seen:
                        seen[key] = candidate
        else:
            sq = enumerate_skeletons(q)
            for s1 in sp:
                for s2 in sq:
                    for candidate in combine_two_skeletons(s1, s2):
                        # Inherit s1's start_up. Paper Appendix §: if join
                        # orientation conflict, s2 is rotated 180° around X
                        # (= z-flip s2, flipping its orientations). Since we
                        # enumerate Sq with both start_up variants, the
                        # conflict-resolution variant is already covered.
                        candidate = Skeleton(points=candidate.points,
                                             chirality=candidate.chirality,
                                             start_up=s1.start_up)
                        if not is_compact(candidate):
                            continue
                        key = canonical_key(candidate)
                        if key not in seen:
                            seen[key] = candidate
    # No global RCC pass: per the decoded combine-driver logic in
    # CG-2012's Main() and CGMotif::extend(), RCC is applied LOCALLY
    # (per parent's extension batch via `eliminateEquivalentMotifs1`)
    # and GLOBALLY only on the "retained for higher dim" set (motifs
    # with scc=false OR forHigherDim=true). The output for each
    # dimension is *not* globally RCC-deduped, which is why oracle's
    # 41 S4 records include some handedness-equivalent variants from
    # different parents. Our enumeration is at the same level (modulo
    # the +1 discrepancy we never closed).
    return list(seen.values())


def enumerate_skeletons_rcc(n: int) -> list[Skeleton]:
    """Same as `enumerate_skeletons(n)` but applies RCC dedup as a final pass.

    RCC (paper Appendix §1.2 / decoded from CG-2012 `CGMotif::
    checkEquivalence5Gr`) collapses handedness-equivalent skeletons to a
    single canonical representative chosen by lowest (dist_sum, layers,
    original index).

    Empirically this is *more aggressive* than CG-2012's oracle: S4
    drops to 21 (oracle: 41), S5 to 117 (oracle: 648). The mechanical
    port matches the decoded IL precisely; the over-collapse suggests
    CG-2012 applies RCC only within local combine contexts (e.g.,
    siblings of the same parent during growth), not as a global final
    pass. We expose this as a separate function rather than baking it
    into `enumerate_skeletons` so callers can compare both views.
    """
    return rcc_dedup(enumerate_skeletons(n))


def enumerate_dim(n: int) -> Iterator[Skeleton]:
    """Top-level: dispatch to dimension-specific generators where implemented."""
    if n == 1:
        yield _SEED_S1
        return
    if n == 2:
        yield _SEED_S2
        return
    if n == 3:
        yield from enumerate_s3()
        return
    raise NotImplementedError(f"enumerate_dim({n}); pending")
