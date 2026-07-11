"""Phase 3 — paper-faithful S5 query generation (geometric SCC-2 basis).

Locks the record count + per-grid breakdown, the grid-e two-variant emission,
the f conditional handedness, and the core invariants (no handedness line on a
coplanar triple; well-formed matrices; deterministic labeling).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from itertools import combinations

from ssp_enum.compactness import set_scc2_mode
from ssp_enum.enumerate import enumerate_skeletons
from ssp_enum.assignment import skeletons_to_records
from ssp_enum.combine import handedness_signature
from ssp_enum.geom_scc2 import grid_of, reference_labeling, mandatory_handedness_gap

TRI = list(combinations(range(1, 6), 3))
IDX = {t: i for i, t in enumerate(TRI)}


def _geom_s5():
    try:
        set_scc2_mode("geometric")
        return enumerate_skeletons(5)
    finally:
        set_scc2_mode("graph")


def test_reference_labeling_is_bijection_for_all_survivors():
    for s in _geom_s5():
        phi = reference_labeling(s)
        assert phi is not None
        assert sorted(phi.keys()) == [1, 2, 3, 4, 5]
        assert sorted(phi.values()) == [1, 2, 3, 4, 5]


def test_paper_faithful_count_and_grid_breakdown():
    s5 = _geom_s5()
    gid = {i: grid_of(s) for i, s in enumerate(s5)}
    recs = list(skeletons_to_records(s5, paper_faithful=True))
    assert len(recs) == 708
    from collections import Counter
    assert dict(Counter(gid[r.skeleton_id] for r in recs)) == \
        {"f": 480, "e": 160, "gh": 56, "d": 12}


def test_no_handedness_line_on_a_coplanar_triple():
    s5 = _geom_s5()
    sig = {i: handedness_signature(s) for i, s in enumerate(s5)}
    for r in skeletons_to_records(s5, paper_faithful=True):
        for (i, j, k, _lr) in r.handedness:
            assert sig[r.skeleton_id][IDX[(i, j, k)]] != 0


def test_matrices_well_formed():
    s5 = _geom_s5()
    for r in skeletons_to_records(s5, paper_faithful=True):
        for i in range(5):
            assert r.matrix[i][i] == "*"
            for j in range(i + 1, 5):
                assert r.matrix[i][j] != ""


def test_grid_e_emits_two_differing_variants():
    s5 = _geom_s5()
    gid = {i: grid_of(s) for i, s in enumerate(s5)}
    from collections import defaultdict
    groups = defaultdict(dict)
    for r in skeletons_to_records(s5, paper_faithful=True):
        if gid[r.skeleton_id] == "e":
            groups[(r.skeleton_id, r.sub_first)][r.sub_second] = r
    assert groups  # grid e present
    for variants in groups.values():
        assert set(variants) == {0, 1}
        assert variants[0].matrix != variants[1].matrix          # disjunction differs
        assert variants[0].handedness == variants[1].handedness  # handedness identical


def test_f_conditional_handedness_activates_on_HHHEE_in_reference():
    s5 = _geom_s5()
    gid = {i: grid_of(s) for i, s in enumerate(s5)}
    saw_9 = False
    for r in skeletons_to_records(s5, paper_faithful=True):
        if gid[r.skeleton_id] != "f":
            continue
        phi = reference_labeling(s5[r.skeleton_id])
        inv = {v: k for k, v in phi.items()}
        ref_types = tuple(r.sse_types[inv[a] - 1] for a in range(1, 6))
        cond = ref_types[:3] == ("H", "H", "H") and ref_types[3:] == ("E", "E")
        gap = len(mandatory_handedness_gap(s5[r.skeleton_id]))
        # conditional adds 2 triples; gaps remove some. Without gaps: 7 or 9.
        if cond and gap == 0:
            assert len(r.handedness) == 9
            saw_9 = True
        elif not cond and gap == 0:
            assert len(r.handedness) == 7
    assert saw_9


def test_gap_survivors_are_the_28_residual():
    s5 = _geom_s5()
    gaps = [grid_of(s) for s in s5 if mandatory_handedness_gap(s)]
    from collections import Counter
    assert dict(Counter(gaps)) == {"f": 21, "gh": 7}
