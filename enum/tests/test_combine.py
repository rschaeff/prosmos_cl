"""Combine-pairs growth tests.

Phase A scope: single-node combine (`|s2| = 1`). These tests pin the
current behavior and document the gap to Phase B / oracle.
"""

from __future__ import annotations

from ssp_enum.combine import (
    adjacent_points,
    canonical_key,
    combine_two_skeletons,
    combine_with_single_node,
    valid_extension_points,
)
from ssp_enum.compactness import is_compact
from ssp_enum.enumerate import enumerate_s3_planar, enumerate_skeletons
from ssp_enum.lattice import LatticePoint
from ssp_enum.skeleton import Skeleton


def _adjacency_signature(skel):
    """Compact representation: upper-triangular adjacency as tuple of bools."""
    m = skel.adjacency_matrix()
    return tuple(
        m[i][j] for i in range(skel.dim) for j in range(i + 1, skel.dim)
    )


def test_valid_extension_points_for_s2():
    """For the canonical S2 (two adjacent points), valid(S2) should be:
    2 mutual-neighbors (triangle-completing) + 2 collinear extensions = 4.
    """
    s2 = Skeleton(points=(LatticePoint(0, 0, 0), LatticePoint(1, 0, 0)))
    valid = valid_extension_points(s2)
    # All four should be at z=0
    assert all(p.z == 0 for p in valid), valid
    # The 4 specific lattice points
    expected = {
        LatticePoint(-1, 0, 0),   # collinear left
        LatticePoint(2, 0, 0),    # collinear right
        LatticePoint(0, 1, 0),    # triangle above
        LatticePoint(1, -1, 0),   # triangle below
    }
    assert valid == expected


def test_adjacent_vs_valid_at_s2():
    """adjacent(s) is broader than valid(s): cross-layer points are
    adjacent but not valid (single neighbor, not collinear)."""
    s2 = Skeleton(points=(LatticePoint(0, 0, 0), LatticePoint(1, 0, 0)))
    adj = adjacent_points(s2)
    valid = valid_extension_points(s2)
    cross_layer_above = LatticePoint(0, 0, 1)
    assert cross_layer_above in adj
    assert cross_layer_above not in valid


def test_combine_with_single_node_yields_both_join_modes():
    """For each valid point, combine yields both end-join (new label last)
    and front-join (new label first; old labels shifted)."""
    s2 = Skeleton(points=(LatticePoint(0, 0, 0), LatticePoint(1, 0, 0)))
    candidates = list(combine_with_single_node(s2))
    # 4 valid points × 2 join modes = 8 raw candidates (before dedup)
    assert len(candidates) == 8
    # Confirm both end and front variants exist for some specific point
    p = LatticePoint(2, 0, 0)
    end_join = Skeleton(points=s2.points + (p,))
    front_join = Skeleton(points=(p,) + s2.points)
    assert end_join in candidates
    assert front_join in candidates


def test_canonical_key_translation_invariance():
    """Translating a skeleton in XY gives the same canonical key."""
    s = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(0, 1, 0),
    ))
    s_translated = Skeleton(points=(
        LatticePoint(5, -3, 0),
        LatticePoint(6, -3, 0),
        LatticePoint(5, -2, 0),
    ))
    assert canonical_key(s) == canonical_key(s_translated)


def test_canonical_key_rotation_invariance():
    """A 60° hex rotation of a skeleton gives the same canonical key."""
    s = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(2, 0, 0),
    ))
    # Rotate by 60°: (q,r) -> (-r, q+r)
    s_rotated = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(0, 1, 0),
        LatticePoint(0, 2, 0),
    ))
    assert canonical_key(s) == canonical_key(s_rotated)


def test_canonical_key_distinguishes_different_sequence_labels():
    """Linear (1—2—3 collinear) and Bent-B (1 in spatial middle) have
    the same UNLABELED grid but DIFFERENT labelings → distinct keys."""
    linear = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(2, 0, 0),
    ))
    bent_b = Skeleton(points=(
        LatticePoint(1, 0, 0),
        LatticePoint(0, 0, 0),
        LatticePoint(2, 0, 0),
    ))
    assert canonical_key(linear) != canonical_key(bent_b)


def test_enumerate_skeletons_s3_matches_base_shapes():
    """S3 via single-node combine should reproduce the 4 base shapes
    (Linear, Triangle, Bent-A, Bent-B) from enumerate_s3_planar."""
    combine_s3 = enumerate_skeletons(3)
    planar_s3 = list(enumerate_s3_planar())
    assert len(combine_s3) == 4
    assert len(planar_s3) == 4
    # Same adjacency-pattern set
    combine_patterns = {_adjacency_signature(s) for s in combine_s3}
    planar_patterns = {_adjacency_signature(s) for s in planar_s3}
    assert combine_patterns == planar_patterns


def test_enumerate_skeletons_s4_count():
    """Phase B regression pin: S4 = 10 via combine (single-node ∪ S2+S2).

    S2+S2 produces 4 unique skeletons that all turn out canonically
    equivalent to skeletons already reachable via S3+S1 (any compact
    4-point lattice arrangement is reachable by growing one node onto
    some compact 3-point arrangement). So the count stays at 10.

    Oracle has 41 — the remaining gap to Phase C: handedness/chirality
    variants (canonical_key currently quotients by reflections + z-flip,
    collapsing mirror pairs that CG-2012 keeps distinct).
    """
    assert len(enumerate_skeletons(4)) == 10


def test_enumerate_skeletons_s5_count():
    """Phase B regression pin: S5 = 41 via combine (single-node ∪ S3+S2).

    Up from 31 in Phase A — S3+S2 adds 10 S5 skeletons unreachable
    from S4+S1 alone. Oracle has 648 distinct skel_ids; the gap is
    Phase C (handedness equivalence + mirror/z-flip variants kept
    distinct).
    """
    assert len(enumerate_skeletons(5)) == 41


def test_s2_plus_s2_produces_compact_candidates():
    """Sanity: combine_two_skeletons on S2+S2 must emit candidates that
    pass is_compact. Guards against silent geometric/positioning bugs."""
    s2 = enumerate_skeletons(2)[0]
    raw = list(combine_two_skeletons(s2, s2))
    assert len(raw) > 0
    compact = [c for c in raw if is_compact(c)]
    assert len(compact) >= 4, f"only {len(compact)} compact S2+S2 candidates"


def test_s3_plus_s2_emits_new_s5_skeletons():
    """S3+S2 should contribute skeletons not reachable from S4+S1 alone.

    We can't easily isolate "S3+S2 contribution" through enumerate_skeletons
    (it iterates all splits), so check the lower-bound: total S5 (41)
    exceeds single-node-only (31)."""
    # The 31 figure is the Phase A baseline (single-node growth from S4).
    # If S3+S2 contributed nothing, total would still be 31.
    assert len(enumerate_skeletons(5)) > 31


def test_enumerate_skeletons_all_compact():
    """Every skeleton produced by combine must satisfy PCC + SCC-1 + SCC-2."""
    for n in (3, 4, 5):
        for s in enumerate_skeletons(n):
            assert is_compact(s), f"S{n} skeleton {s} not compact"
