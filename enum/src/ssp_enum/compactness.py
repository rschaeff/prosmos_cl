"""Compactness criteria from Chitturi 2016 Methods + Appendix Fig. S2.

Three criteria define compactness (paraphrased from the paper):

  PCC ("primary compactness"): any lattice point that is *not* part of
       the skeleton can have at most 3 adjacent points in the skeleton.
       Counterexample: Appendix Fig. S2c(ii) shows a trapezoid where an
       interior non-skeleton point has 4 skeleton neighbors — rejected.

  SCC criterion-1 ("secondary compactness, criterion 1"): when a new
       node is added, it must either (a) extend an existing collinear
       set of two or more points (a "line" through the lattice), or (b)
       sit at a lattice position with at least two adjacent existing
       skeleton points. Counterexamples: Appendix Fig. S2c(iii-iv) —
       new node has only one neighbor and doesn't extend a collinear set.

  SCC criterion-2 ("secondary compactness, criterion 2"): the
       candidate's induced grid must match one of a predefined set of
       allowed grids for that dimension (Appendix Fig. S3). At S3, the
       allowed grids are 3-collinear and equilateral triangle; at S4, 3
       grids; at S5, 5 grids in the paper (4 under unlabeled-adjacency
       canonicalization, since the Appendix's Grid 5 / Grid 4 (S3-h)
       and Grid 4 (S3-g) share an unlabeled adjacency graph).
       Counterexample: Appendix Fig. S2d(ii) passes PCC and SCC-1 but
       has an off-whitelist induced grid → rejected.
"""

from __future__ import annotations

from .grids import is_in_whitelist
from .geom_scc2 import passes_geom_scc_2
from .lattice import LatticePoint
from .skeleton import Skeleton

# SCC-2 backend selector. The historic default is the *graph*-based whitelist
# (`grids.is_in_whitelist`), on which the current S5=198 darkness/negspace
# analysis rests. `"geometric"` switches to Fig-S3 hex-congruence
# (`geom_scc2.passes_geom_scc_2`), the paper-faithful test (S5 -> 140). Kept as
# a module toggle so the enumerator (which calls `is_compact` with no params
# through a recursion) can be flipped without threading a parameter, and so
# A/B comparison and regression stay one line apart.
_SCC2_MODE = "graph"  # "graph" | "geometric"


def set_scc2_mode(mode: str) -> str:
    """Set the SCC-2 backend ("graph" or "geometric"); returns the previous mode."""
    global _SCC2_MODE
    if mode not in ("graph", "geometric"):
        raise ValueError(f"unknown SCC-2 mode: {mode!r}")
    prev, _SCC2_MODE = _SCC2_MODE, mode
    return prev


def passes_pcc(skel: Skeleton) -> bool:
    """Perimeter Compactness Criterion (paper Appendix Fig. S2c(ii)).

    Every lattice point not in the skeleton must have ≤3 adjacent
    skeleton points. We check every distinct neighbor of every skeleton
    point — that's the only place a non-skeleton point can have any
    skeleton neighbors at all.
    """
    skel_set = set(skel.points)
    candidates: set[LatticePoint] = set()
    for p in skel.points:
        for n in p.neighbors():
            if n not in skel_set:
                candidates.add(n)
    for c in candidates:
        adj_count = sum(1 for n in c.neighbors() if n in skel_set)
        if adj_count > 3:
            return False
    return True


def _are_collinear(p1: LatticePoint, p2: LatticePoint, p3: LatticePoint) -> bool:
    """True iff three lattice points lie on a straight line on the hex lattice.

    On the axial-coord hex lattice, three points are collinear iff their
    Δq:Δr ratios are equal across all pairs (with Δz=0 throughout). Three
    integer points lie on a line iff the cross product of (p2-p1) x (p3-p1)
    is zero on each plane component.
    """
    if not (p1.z == p2.z == p3.z):
        return False
    d12q, d12r = p2.q - p1.q, p2.r - p1.r
    d13q, d13r = p3.q - p1.q, p3.r - p1.r
    return d12q * d13r - d12r * d13q == 0


def passes_scc_1(skel: Skeleton) -> bool:
    """Structural Compactness Criterion 1 (paper Appendix Fig. S2c(iii-iv)).

    For each node added after the seed pair (i.e., index >= 2 in 0-based,
    >= 3 in paper-1-based), require it to either:
      (a) be adjacent to ≥2 existing skeleton points, OR
      (b) extend a 2+-point collinear line in the existing skeleton.

    The seed pair (indices 0 and 1) is trivially compact: a single edge.
    """
    if skel.dim < 3:
        return True
    for i in range(2, skel.dim):
        added = skel.points[i]
        existing = skel.points[:i]
        adj_existing = [p for p in existing if added.is_adjacent(p)]
        if len(adj_existing) >= 2:
            continue
        if len(adj_existing) == 1:
            # Does `added` extend a 2+-point collinear line ending at adj_existing[0]?
            anchor = adj_existing[0]
            extends = False
            for q in existing:
                if q is anchor:
                    continue
                if anchor.is_adjacent(q) and _are_collinear(q, anchor, added):
                    extends = True
                    break
            if not extends:
                return False
        else:
            # Zero adjacent existing — not connected
            return False
    return True


def passes_scc_2(skel: Skeleton) -> bool:
    """Secondary Compactness Criterion 2 (Appendix § "secondary compactness criterion" + Fig. S3).

    The skeleton's induced grid (the unlabeled set of lattice positions
    it occupies) must be one of the allowed grids for its dimension.
    The whitelist is held in `grids.WHITELIST` per dimension. At dim ≤ 2
    the lattice arrangement is trivially unique; at dim ≥ 6 we don't
    have a whitelist defined yet, so SCC-2 vacuously passes.

    Backend is selected by `_SCC2_MODE` (see `set_scc2_mode`): the default
    `"graph"` uses the unlabeled-adjacency whitelist; `"geometric"` uses
    Fig-S3 hex-congruence (the paper-faithful test).
    """
    if _SCC2_MODE == "geometric":
        return passes_geom_scc_2(skel)
    return is_in_whitelist(skel)


def is_compact(skel: Skeleton) -> bool:
    """Combined compactness predicate: PCC ∧ SCC-1 ∧ SCC-2.

    SCC-2 is a no-op below dim 3 and above the dimensions for which we
    have a whitelist (currently 3, 4, 5). Behavior for dim 3 in our S3
    enumeration is unchanged: every shape we produce (P3 or K3) is on
    the S3 whitelist.
    """
    return passes_pcc(skel) and passes_scc_1(skel) and passes_scc_2(skel)
