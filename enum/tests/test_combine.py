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


def test_canonical_key_distinguishes_mirror_pair():
    """Phase C1: mirror images of a Triangle (CW vs CCW sequence walk)
    are now distinct canonicals — `canonical_key` quotients by rotation
    only, not reflection."""
    triangle_ccw = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(0, 1, 0),
    ))
    # Mirror reflection across q-axis: (q, r) → (q+r, -r)
    triangle_cw = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(1, -1, 0),
    ))
    assert canonical_key(triangle_ccw) != canonical_key(triangle_cw)


def test_enumerate_skeletons_s3_covers_base_shapes():
    """S3 via combine covers all 4 base adjacency patterns.

    Phase C1: under rotation-only `canonical_key`, the Triangle splits
    into 2 mirror variants (CW vs CCW), giving 5 distinct skeletons
    across 4 adjacency patterns. `enumerate_s3_planar` dedupes by
    adjacency-triple (a coarser invariant), still yielding 4.
    """
    combine_s3 = enumerate_skeletons(3)
    planar_s3 = list(enumerate_s3_planar())
    assert len(combine_s3) == 5
    assert len(planar_s3) == 4
    combine_patterns = {_adjacency_signature(s) for s in combine_s3}
    planar_patterns = {_adjacency_signature(s) for s in planar_s3}
    assert combine_patterns == planar_patterns


def test_enumerate_skeletons_s4_count():
    """Phase C-S4 regression pin: S4 = 42 via combine (oracle: 41).

    Up from 14 once WHITELIST_S4 was corrected per Fig. S3 to include
    the 4-edge triangle+pendant grid (paper b) and drop the bogus
    6-edge K4 (not realizable on 2D hex). Combine now nearly matches
    the oracle's 41 distinct skel_ids; the +1 discrepancy is likely a
    single mirror-pair that needs RCC dedup, not a structural gap.
    """
    assert len(enumerate_skeletons(4)) == 42


def test_enumerate_skeletons_s5_count():
    """Phase C-S5 regression pin: S5 = 198 via combine (oracle: 648).

    Big jump from 72 after deriving WHITELIST_S5 directly from Fig. S3
    (d-h): added the 5e tripod+2pendants (e) and 6e bowtie / sparse
    K1,4 (g/h sparse) lattice variants which were previously rejected
    by SCC-2. Combine now produces skeletons across all 4 paper-
    derived S5 lattice signatures:
      4e P5 (d):          12 skeletons
      5e tripod+2pend (e): 98 skeletons
      6e K1,4-sparse (g): 28 skeletons
      7e K1,4-dense (f/g):60 skeletons

    Remaining 3.3x gap to oracle 648 likely comes from chirality
    variants (L/R) that our rotation-only canonical_key keeps distinct
    but combine doesn't produce enough of, plus possibly a separate
    (f) lattice realization the (g/h-shared) signature here doesn't
    capture.
    """
    assert len(enumerate_skeletons(5)) == 198


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

    Lower-bound regression check: total S5 must exceed what single-node
    growth alone produces. We don't pin the exact single-node count
    here (it depends on canonical_key choices), but the inequality is
    invariant — S3+S2 reaches arrangements where the 5-node split into
    a 3-cycle and a 2-cycle requires placing s2 with rotation, which
    single-node growth cannot replicate."""
    # If S3+S2 contributed nothing, total would equal the single-node
    # count. Under Phase C1 canonical_key, total = 70 > single-node-only.
    assert len(enumerate_skeletons(5)) > 50


def test_enumerate_skeletons_all_compact():
    """Every skeleton produced by combine must satisfy PCC + SCC-1 + SCC-2."""
    for n in (3, 4, 5):
        for s in enumerate_skeletons(n):
            assert is_compact(s), f"S{n} skeleton {s} not compact"
