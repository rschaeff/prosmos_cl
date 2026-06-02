"""Tests for ssp_enum.assignment — Skeleton -> SSPRecord adapter."""

from itertools import product

from ssp_enum.assignment import skeletons_to_records, _find_sheets
from ssp_enum.enumerate import enumerate_skeletons
from ssp_enum.prosmos import write_query


# ---------------------------------------------------------------------------
# Cardinality: full Cartesian SSE typing yields 2**n records per skeleton.
# ---------------------------------------------------------------------------

def test_cardinality_s3():
    skels = enumerate_skeletons(3)
    records = list(skeletons_to_records(skels))
    assert len(records) == len(skels) * (2 ** 3)


def test_cardinality_s4():
    skels = enumerate_skeletons(4)
    records = list(skeletons_to_records(skels))
    assert len(records) == len(skels) * (2 ** 4)


def test_cardinality_s5():
    skels = enumerate_skeletons(5)
    records = list(skeletons_to_records(skels))
    assert len(records) == len(skels) * (2 ** 5)


# ---------------------------------------------------------------------------
# Type coverage: every (skeleton, type-assignment) pair appears exactly once.
# ---------------------------------------------------------------------------

def test_every_type_assignment_present_per_skeleton():
    skels = enumerate_skeletons(3)
    expected_types = set(product(("H", "E"), repeat=3))
    by_skel: dict[int, set] = {}
    for r in skeletons_to_records(skels):
        by_skel.setdefault(r.skeleton_id, set()).add(r.sse_types)
    assert len(by_skel) == len(skels)
    for s_id, type_set in by_skel.items():
        assert type_set == expected_types


# ---------------------------------------------------------------------------
# Matrix shape: n x n, diagonal '*', strict upper triangle non-empty, lower
# triangle empty (the convention `prosmos.write_query` depends on).
# ---------------------------------------------------------------------------

def test_matrix_well_formed():
    skels = enumerate_skeletons(4)
    for r in skeletons_to_records(skels):
        assert len(r.matrix) == 4
        for i, row in enumerate(r.matrix):
            assert len(row) == 4
            assert row[i] == "*", f"diagonal not * at {i}: {row[i]!r}"
            for j in range(i):
                assert row[j] == "", f"lower triangle non-empty at ({i},{j}): {row[j]!r}"
            for j in range(i + 1, 4):
                assert row[j] != "", f"upper triangle empty at ({i},{j})"
                assert row[j] in ("c", "t", "u", "v", "C", "T", "X", "-"), \
                    f"unexpected code {row[j]!r}"


# ---------------------------------------------------------------------------
# Sheet detection: explicit small case.
# ---------------------------------------------------------------------------

def test_find_sheets_simple_chain():
    # 4 nodes in a row, all E. Adjacency: 1-2, 2-3, 3-4. One sheet (1,2,3,4).
    types = ("E", "E", "E", "E")
    adj = [
        [False, True, False, False],
        [True, False, True, False],
        [False, True, False, True],
        [False, False, True, False],
    ]
    sheets = _find_sheets(types, adj)
    assert sheets == [[1, 2, 3, 4]]


def test_find_sheets_two_disjoint():
    # 4 E nodes split into two pairs: 1-2 adj, 3-4 adj, 2-3 not adj.
    types = ("E", "E", "E", "E")
    adj = [
        [False, True, False, False],
        [True, False, False, False],
        [False, False, False, True],
        [False, False, True, False],
    ]
    sheets = _find_sheets(types, adj)
    assert sheets == [[1, 2], [3, 4]]


def test_find_sheets_ignores_helices():
    # Mixed types: only Es participate.
    types = ("E", "H", "E", "E")
    adj = [
        [False, True, True, False],
        [True, False, True, False],
        [True, True, False, True],
        [False, False, True, False],
    ]
    # Although 1 is adjacent to 3 (both E) and 3 to 4 (both E), 1-2/2-3 are
    # E-H/E-H so the helix doesn't link sheets. Sheet 1: {1,3,4}.
    sheets = _find_sheets(types, adj)
    assert sheets == [[1, 3, 4]]


# ---------------------------------------------------------------------------
# All-H skeletons: no sheets, no sheetS/sheetD.
# ---------------------------------------------------------------------------

def test_all_helix_no_sheet_directives():
    skels = enumerate_skeletons(4)
    for r in skeletons_to_records(skels):
        if r.sse_types == ("H", "H", "H", "H"):
            assert r.same_sheet == ()
            assert r.diff_sheet == ()


# ---------------------------------------------------------------------------
# Round-trip through write_query — no exception, query has the expected
# structural lines.
# ---------------------------------------------------------------------------

def test_write_query_roundtrip_s3():
    skels = enumerate_skeletons(3)
    for r in skeletons_to_records(skels):
        out = write_query(r)
        lines = out.splitlines()
        assert lines[0] == "1 2 3"
        assert lines[1] == " ".join(r.sse_types)
        # Per-SSE length lines must appear.
        n_length = sum(1 for ln in lines if ln.startswith("length"))
        assert n_length == 3


def test_handedness_emission_tracks_signature():
    """Records inherit per-triple handedness from `handedness_signature`,
    which is non-zero whenever the SSE-direction alternation gives a
    non-coplanar triple. We should see both empty-handedness records
    (coplanar / collinear shapes) and non-empty ones."""
    from ssp_enum.combine import handedness_signature
    skels = enumerate_skeletons(3)
    seen_with_hand = False
    seen_without_hand = False
    for r in skeletons_to_records(skels):
        sig = handedness_signature(skels[r.skeleton_id])
        if any(s != 0 for s in sig):
            assert r.handedness, f"sig {sig} non-zero but no handedness emitted"
            seen_with_hand = True
        else:
            assert not r.handedness, f"sig all-zero but handedness emitted {r.handedness}"
            seen_without_hand = True
    assert seen_with_hand
    assert seen_without_hand
