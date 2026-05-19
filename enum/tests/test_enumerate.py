"""Enumeration tests: our fresh implementation vs. CG-2012 expectations."""

from __future__ import annotations

from collections import Counter

from ssp_enum.enumerate import enumerate_dim, enumerate_s3, enumerate_s3_planar


def test_s1_count():
    assert len(list(enumerate_dim(1))) == 1


def test_s2_count():
    assert len(list(enumerate_dim(2))) == 1


def test_s3_planar_count():
    """The 4 distinct planar S3 skeletons (linear + triangle + 2 bents).

    CG-2012 reports 11 S3 SSPs total; that count includes handedness
    variations via Z-displacement (3-1, 3-2 = triangle L/R; 3-5, 3-6 =
    linear L/R; etc.). The planar (Z=0) base set is the 4 spatial-sequence
    arrangements documented in CG-2012/S3/Stru.txt: linear, triangle,
    bent-A (sequence-3 in spatial middle), bent-B (sequence-1 in middle).
    """
    skeletons = list(enumerate_s3_planar())
    assert len(skeletons) == 4


def test_s3_planar_topologies():
    """The 4 planar S3 skeletons each have a distinct adjacency pattern."""
    skeletons = list(enumerate_s3_planar())
    # Encode each as the upper-triangular adjacency: (1,2), (1,3), (2,3)
    adj_patterns = set()
    for s in skeletons:
        m = s.adjacency_matrix()
        pattern = (m[0][1], m[0][2], m[1][2])
        adj_patterns.add(pattern)

    expected = {
        (True, False, True),    # linear: 1-2 adj, 2-3 adj, 1-3 NOT
        (True, True, True),     # triangle: all three adjacent
        (False, True, True),    # bent-A: 1-2 NOT, 1-3 adj, 2-3 adj
        (True, True, False),    # bent-B: 1-2 adj, 1-3 adj, 2-3 NOT
    }
    assert adj_patterns == expected


def test_s3_full_count_matches_cg2012():
    """CG-2012/S3/Stru.txt reports 11 S3 SSPs (panels 3-0 .. 3-10).

    Three planar/unhanded (Linear, Bent-A, Bent-B) plus eight chiral
    (4 shapes × {L, R}) = 11. The Triangle has no unhanded variant.
    """
    ssps = list(enumerate_dim(3))
    assert len(ssps) == 11
    # enumerate_s3 should agree with enumerate_dim(3)
    assert len(list(enumerate_s3())) == 11


def test_s3_chirality_breakdown():
    """3 unhanded + 4 L + 4 R = 11 (Linear/Bent-A/Bent-B each get None+L+R;
    Triangle gets only L+R)."""
    ssps = list(enumerate_s3())
    counts = Counter(s.chirality for s in ssps)
    assert counts[None] == 3
    assert counts["L"] == 4
    assert counts["R"] == 4


def test_s3_per_shape_chirality():
    """Each acyclic shape yields {None, L, R}; the cyclic Triangle yields {L, R}.

    Maps directly to CG-2012's panel layout (Stru.txt):
      Linear:   3-0 (None), 3-5 (L), 3-6 (R)
      Triangle: 3-1 (L), 3-2 (R)               — no unhanded
      Bent-A:   3-3 (None), 3-7 (L), 3-8 (R)
      Bent-B:   3-4 (None), 3-9 (R), 3-10 (L)
    """
    by_shape: dict[tuple[bool, bool, bool], set[str | None]] = {}
    for s in enumerate_s3():
        m = s.adjacency_matrix()
        shape = (m[0][1], m[0][2], m[1][2])
        by_shape.setdefault(shape, set()).add(s.chirality)

    triangle = (True, True, True)
    open_shapes = {
        (True, False, True),    # Linear
        (False, True, True),    # Bent-A
        (True, True, False),    # Bent-B
    }
    assert by_shape[triangle] == {"L", "R"}
    for shape in open_shapes:
        assert by_shape[shape] == {None, "L", "R"}, shape
