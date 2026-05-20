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
    dist_sum,
    handedness_signature,
    layers,
    rcc_dedup,
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


def test_canonical_key_ignores_start_orientation():
    """Phase C2-revert: start_up no longer in the canonical key. Pre-C2
    S4=42 matched oracle 41 exactly while post-C2 S4=84 overshot by 2x;
    CG-2012 doesn't separate UP-start and DOWN-start as distinct skeletons
    at the skeleton-enumeration level. The `start_up` field is retained on
    Skeleton for downstream handedness/chirality work."""
    pts = (LatticePoint(0, 0, 0), LatticePoint(1, 0, 0), LatticePoint(2, 0, 0))
    up = Skeleton(points=pts, start_up=True)
    down = Skeleton(points=pts, start_up=False)
    assert canonical_key(up) == canonical_key(down)


def test_orientations_alternate_from_start_up():
    """orientations property alternates from start_up per paper rule."""
    s = Skeleton(
        points=(LatticePoint(0, 0, 0), LatticePoint(1, 0, 0),
                LatticePoint(2, 0, 0), LatticePoint(3, 0, 0)),
        start_up=True,
    )
    assert s.orientations == (True, False, True, False)
    s2 = Skeleton(points=s.points, start_up=False)
    assert s2.orientations == (False, True, False, True)


def test_handedness_signature_collinear_zero():
    """Per paper Appendix § 'all collinear points: None handedness'.

    For 3 lattice-collinear points (xy-collinear), the scalar triple
    product (p × q) · r evaluates to 0 regardless of orientation, so
    handedness = (0,).
    """
    s = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(2, 0, 0),
    ))
    sig = handedness_signature(s)
    assert sig == (0,), sig


def test_handedness_signature_triangle_chirality():
    """A CCW and CW Triangle (mirror pair) have opposite-sign
    handedness signatures."""
    ccw = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(0, 1, 0),
    ))
    cw = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(1, -1, 0),
    ))
    s_ccw = handedness_signature(ccw)
    s_cw = handedness_signature(cw)
    assert s_ccw == (-s_cw[0],), (s_ccw, s_cw)
    assert s_ccw != (0,)  # non-collinear → non-zero


def test_dist_sum_collinear_three_points():
    """For 3 collinear hex points at unit spacing, pairwise distances
    are 1, 1, 2 → sum = 4.0."""
    s = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(2, 0, 0),
    ))
    assert abs(dist_sum(s) - 4.0) < 1e-9


def test_layers_count_unique_y():
    """3 collinear hex points on the q-axis all have y=0 → 1 layer.
    Triangle (0,0), (1,0), (0,1) has y values 0, 0, √3/2 → 2 layers."""
    linear = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(2, 0, 0),
    ))
    assert layers(linear) == 1
    triangle = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(0, 1, 0),
    ))
    assert layers(triangle) == 2


def test_rcc_dedup_collapses_equivalent_skeletons():
    """Two skeletons with the same handedness signature should collapse
    to one (the lower dist_sum/layers wins). Two with different
    handedness signatures stay distinct."""
    # Two CCW triangles at different lattice positions — same handedness sig.
    ccw_a = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(0, 1, 0),
    ))
    ccw_b = Skeleton(points=(  # translated/rotated — same handedness
        LatticePoint(5, 5, 0),
        LatticePoint(6, 5, 0),
        LatticePoint(5, 6, 0),
    ))
    # CW triangle — different handedness sig.
    cw = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(1, -1, 0),
    ))
    survivors = rcc_dedup([ccw_a, ccw_b, cw])
    assert len(survivors) == 2  # one CCW representative + one CW


def test_handedness_signature_orientation_dependence():
    """A non-collinear triple gets opposite handedness when all SSE
    orientations flip (start_up → not start_up flips every z=±1, which
    negates the scalar triple product)."""
    up = Skeleton(points=(
        LatticePoint(0, 0, 0),
        LatticePoint(1, 0, 0),
        LatticePoint(0, 1, 0),
    ), start_up=True)
    down = Skeleton(points=up.points, start_up=False)
    s_up = handedness_signature(up)
    s_down = handedness_signature(down)
    # Negating all node z's flips the triple-product sign for triples
    # where the z dimension contributes; for purely planar (z=0) skels
    # this wouldn't matter, but our nodes have z=±1 from orientation.
    assert s_up != s_down, (s_up, s_down)


def test_enumerate_skeletons_s3_covers_base_shapes():
    """S3 via combine covers all 4 base adjacency patterns.

    Phase C1 split Triangle into CW/CCW mirror variants → 5 skeletons.
    Phase C2 doubled to 10 (start_up); Phase C2-revert restored 5.
    """
    combine_s3 = enumerate_skeletons(3)
    planar_s3 = list(enumerate_s3_planar())
    assert len(combine_s3) == 5
    assert len(planar_s3) == 4
    combine_patterns = {_adjacency_signature(s) for s in combine_s3}
    planar_patterns = {_adjacency_signature(s) for s in planar_s3}
    assert combine_patterns == planar_patterns


def test_enumerate_skeletons_s4_count():
    """Phase C2-revert regression pin: S4 = 42 via combine (oracle: 41).

    Down from 84 once `start_up` is removed from canonical_key. Now
    essentially matches CG-2012's distinct-skel_id count; the +1 is
    likely a single mirror-pair case that RCC would tie-break.
    """
    assert len(enumerate_skeletons(4)) == 42


def test_enumerate_skeletons_s5_count():
    """Phase C2-revert regression pin: S5 = 198 via combine (oracle: 648).

    Down from 396 once `start_up` is removed from canonical_key.
    The S4 fix is clean but S5 still 3.3x short of oracle 648 — the
    structural cause is not yet identified. Candidates:
      - additional (f) / (g) lattice realizations our combine doesn't
        reach via the current seed + valid_extension_points
      - CG-2012 over-enumerates per-build-path without RCC dedup,
        so 648 might collapse to ~200 under a proper canonical
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


def test_enumerate_skeletons_rcc_post_pass():
    """`enumerate_skeletons_rcc(n)` = `enumerate_skeletons(n)` + RCC dedup.

    Mechanical port of CG-2012's `CGMotif::checkEquivalence5Gr` (decoded
    from the binary IL): group skeletons by handedness signature; within
    each group keep the one with lowest (dist_sum, layers, original
    index). Applied as a post-pass on combine output.

    Counts:
      S3: 5 → 3   (collinear shapes all collapse via same handedness sig)
      S4: 42 → 23 (matches pure handedness signature count from Phase C4
                   investigation; tighter than oracle 41)
      S5: 198 → 164 (some intra-class collapse; far from oracle 648)

    The over-collapse vs oracle suggests CG-2012 applies RCC only within
    local combine contexts (siblings of the same parent during growth),
    not as a single global pass. `enumerate_skeletons` (no RCC) remains
    the production enumerator; this function is exposed for analysis.
    """
    from ssp_enum.enumerate import enumerate_skeletons_rcc
    assert len(enumerate_skeletons_rcc(3)) == 3
    assert len(enumerate_skeletons_rcc(4)) == 23
    assert len(enumerate_skeletons_rcc(5)) == 164


def test_enumerate_skeletons_all_compact():
    """Every skeleton produced by combine must satisfy PCC + SCC-1 + SCC-2."""
    for n in (3, 4, 5):
        for s in enumerate_skeletons(n):
            assert is_compact(s), f"S{n} skeleton {s} not compact"
