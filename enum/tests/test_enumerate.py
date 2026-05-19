"""Enumeration tests: our fresh implementation vs. CG-2012 expectations."""

from __future__ import annotations

from ssp_enum.enumerate import enumerate_dim, enumerate_s3_planar


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
