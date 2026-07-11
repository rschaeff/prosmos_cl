"""Geometric SCC-2 (Chitturi 2016 Appendix, "secondary compactness criterion 2").

The paper's SCC-2 test is *geometric*: a candidate skeleton passes iff its
induced grid (the unlabeled set of hex-lattice positions it occupies)
**"matches (or is a symmetric variation of) one of these predefined grids"**
(Fig. S3) — i.e. is congruent, under the hex plane's isometry group (6
rotations x mirror) composed with translation, to a Fig-S3 grid.

`grids.py`'s `unlabeled_grid_signature` used a *coarser* test — the unlabeled
adjacency graph. A graph can have several non-congruent hex realizations
(e.g. the straight Fig-S3(e) grid AND a bent 1-2-3-4 path with the same
pendant, which is graph-isomorphic to (e) but not the Fig-S3 shape). The
graph test admits all of them; the geometric test admits only the Fig-S3
shape. This is the root of the 198-vs-paper discrepancy: the extra bent
geometries have no paper grid, hence no paper handedness — a paper
handedness triple becomes geometrically coplanar for them.

This module holds the Fig-S3 grid point-sets (S5 decoded; S3/S4 pending),
a hex-isometry canonical form for point sets, and the geometric SCC-2 test.
"""

from __future__ import annotations

from .skeleton import Skeleton

# ---------------------------------------------------------------------------
# Hex axial-coord isometries (D6 = 6 rotations x mirror).  z is not involved:
# an induced grid is a 2D point set in the hex plane.
# ---------------------------------------------------------------------------

Point = tuple[int, int]
PointSet = tuple[Point, ...]


def _rot60(q: int, r: int) -> Point:
    return (-r, q + r)


def _reflect(q: int, r: int) -> Point:
    return (q + r, -r)


def _isometry_images(pts: list[Point]) -> list[list[Point]]:
    """The 12 D6 images (6 rotations x mirror) of a point set."""
    out: list[list[Point]] = []
    cur = list(pts)
    for _ in range(6):
        out.append(list(cur))
        out.append([_reflect(q, r) for q, r in cur])
        cur = [_rot60(q, r) for q, r in cur]
    return out


def geometric_canonical(pts) -> PointSet:
    """Canonical form of an unlabeled hex point set under D6 x translation.

    Returns the lex-min, translation-normalized, sorted tuple of (q, r)
    over all 12 hex isometries. Two point sets are hex-congruent iff their
    canonical forms are equal.
    """
    plist = [(int(q), int(r)) for q, r in pts]
    best: PointSet | None = None
    for img in _isometry_images(plist):
        s = sorted(img)
        q0, r0 = s[0]
        t = tuple(sorted((q - q0, r - r0) for q, r in img))
        if best is None or t < best:
            best = t
    return best if best is not None else ()


# ---------------------------------------------------------------------------
# Fig-S3 grids.  S5 decoded from the SI figure (page 11, gs -r400); S3/S4
# still to be read off the same figure (panels a-c) — see scc2_geometric_plan.
# Each grid: reference labeling as axial (q, r), node i at index i-1.
# ---------------------------------------------------------------------------

FIG_S3_S5: dict[str, list[Point]] = {
    # (d) P5 path: 1-2-3-4-5 collinear.                         4 edges.
    "d": [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0)],
    # (e) 1-2-3-4 collinear + node 5 in the 2-3 corner.         5 edges.
    "e": [(0, 0), (1, 0), (2, 0), (3, 0), (1, 1)],
    # (f) dense: 1-2-3 collinear, 4 & 5 on the next row.        7 edges.
    "f": [(0, 0), (1, 0), (2, 0), (1, 1), (0, 1)],
    # (g/h) star: centre 2 with 4 leaves 1,3,4,5.               6 edges.
    "gh": [(-1, 0), (0, 0), (1, 0), (0, 1), (0, -1)],
}


# ---------------------------------------------------------------------------
# Per-grid MINIMAL handedness (Chitturi 2016, "minimum number of handedness"
# — Fig-S3 caption philosophy applied to handedness). Triples are 1-based in
# the reference labeling above. Read from the CG-2012 oracle (reference/
# IA-S5.txt, which reproduces the paper's per-count distribution
# {0:45, 2:440, 6:1378, 7:826, 9:118} exactly) cross-checked against our own
# geometric handedness signature. See scc2_handedness_stepB.md for derivation.
#
# The listed triples are the ones that carry an L/R constraint; the L/R value
# itself is geometric (emitted per-skeleton from `combine.handedness_signature`
# at query-build time). A grid's full geometric non-coplanar set can be larger
# than its minimal list — the paper minimizes to the smallest set that keeps
# SSPs from *distinct* grids from sharing hits.
#
#   d  : 0 constraints.
#   e  : 2 (all contain the pendant node 5). Three symmetry-equivalent minimal
#        pairs exist in the oracle; the reference-labeling pick is recorded as
#        MANDATORY below (residual: exact geometric block-match to our labeling
#        is still open — any of the three is a valid size-2 disambiguator).
#   g/h: 6 = the complete geometric non-coplanar set (minimal == full).
#   f  : 7 MANDATORY + 2 CONDITIONAL. The conditional pair {(1,4,5),(3,4,5)}
#        turns on ONLY for the typing 1-2-3 = H, 4-5 = E (nodes 4,5 form a
#        strand sheet) — verified in oracle block skel=41: typings EEEHH/EEEEH/
#        EEEHE/EEEEE/HHHHH/HHHEH/HHHHE -> 7 lines; typing HHHEE -> 9 lines.
#        Matches SI ("1,2,3 helices & 4,5 strands -> sheet").

Triple = tuple[int, int, int]

HANDEDNESS_S5_MANDATORY: dict[str, tuple[Triple, ...]] = {
    "d": (),
    "e": ((1, 4, 5), (2, 3, 5)),
    "gh": ((1, 2, 4), (1, 2, 5), (1, 3, 4), (1, 3, 5), (2, 3, 4), (2, 3, 5)),
    "f": ((1, 2, 4), (1, 2, 5), (1, 3, 4), (1, 3, 5),
          (2, 3, 4), (2, 3, 5), (2, 4, 5)),
}

# Conditional handedness: {grid: (triples, condition)}. `condition` names the
# SSE-typing that activates the extra triples (grid f only, at S5).
HANDEDNESS_S5_CONDITIONAL: dict[str, tuple[tuple[Triple, ...], str]] = {
    "f": (((1, 4, 5), (3, 4, 5)), "types==HHHEE (1-2-3 helix, 4-5 strand-sheet)"),
}


# ---------------------------------------------------------------------------
# Per-grid MINIMAL interactions (Chitturi 2016 Fig-S3 line styles). Read
# directly off Fig S3 (SI page 11); see scc2_interactions_stepB.md. Verified:
# for every grid, REQUIRED ∪ DISJUNCTION ∪ OPTIONAL-X partitions the grid's
# full lattice-edge set exactly.
#
# Fig-S3 caption: "Solid lines indicate an interaction. In Fig 5(e) two broken
# lines indicate that at least one of the interactions is mandatory. Through
# experimentation we determined the minimum number of (specific) interactions
# for each grid so that SSPs derived from skeletons (that induce distinct grids)
# do not have common hits."
#
#   REQUIRED  (solid) : concrete interaction demanded in the SSP query.
#   DISJUNCTION (broken): "≥1 mandatory" set — at least one must hold (grid e).
#   optional (X)      : lattice-adjacent but not required (the remaining grid
#                       edges); marked X/wildcard in the query.
#
# Edges are 1-based in the reference labeling (`FIG_S3_S5`), which matches the
# Fig-S3 node labels for every S5 grid.

INTERACTIONS_S5_REQUIRED: dict[str, tuple[tuple[int, int], ...]] = {
    "d": ((1, 2), (2, 3), (3, 4), (4, 5)),               # P5 path
    "e": ((1, 2), (2, 3), (3, 4)),                       # collinear run
    "f": ((1, 2), (2, 3), (3, 4), (4, 5), (1, 5)),       # 5-cycle ring
    "gh": ((1, 2), (2, 3), (2, 4), (2, 5)),              # K1,4 star (centre 2)
}

# "at least one of these is mandatory" (broken lines). Grid e only, at S5.
INTERACTIONS_S5_DISJUNCTION: dict[str, tuple[tuple[int, int], ...]] = {
    "e": ((2, 5), (3, 5)),
}

# Lattice-adjacent but optional (X). = grid lattice edges minus required minus
# disjunction. f: the two spokes into the centre of the ring; g/h: leaf-leaf.
INTERACTIONS_S5_OPTIONAL: dict[str, tuple[tuple[int, int], ...]] = {
    "d": (),
    "e": (),
    "f": ((2, 4), (2, 5)),
    "gh": ((1, 5), (3, 4)),
}


def _build_whitelist(grids: dict[str, list[Point]]) -> dict[PointSet, str]:
    return {geometric_canonical(v): k for k, v in grids.items()}


GEOM_WHITELIST: dict[int, dict[PointSet, str]] = {
    5: _build_whitelist(FIG_S3_S5),
}


def grid_of(skel: Skeleton) -> str | None:
    """The Fig-S3 grid name whose shape `skel`'s induced grid is congruent to.

    Returns None if the skeleton's point set matches no Fig-S3 grid at its
    dimension, or if no whitelist is defined for that dimension.
    """
    wl = GEOM_WHITELIST.get(skel.dim)
    if wl is None:
        return None
    return wl.get(geometric_canonical((p.q, p.r) for p in skel.points))


def passes_geom_scc_2(skel: Skeleton) -> bool:
    """Geometric SCC-2: dims with a Fig-S3 whitelist require congruence;
    dims without one pass vacuously (same convention as the graph test)."""
    if skel.dim not in GEOM_WHITELIST:
        return True
    return grid_of(skel) is not None


# ---------------------------------------------------------------------------
# Congruence labeling: map a survivor skeleton's sequence labels onto the
# reference-grid labels. This is what lets the per-grid reference constraints
# (handedness triples, interaction edges) be transported onto any survivor.
# ---------------------------------------------------------------------------

def _ordered_isometry_images(pts: list[Point]) -> list[list[Point]]:
    """Like `_isometry_images` but the caller relies on element order being
    preserved (each image is the isometry applied element-wise)."""
    return _isometry_images(pts)


def _valid_labelings(skel: Skeleton, grid: str) -> list[dict[int, int]]:
    """All distinct sequence-label -> reference-label maps from a hex isometry
    + translation carrying `skel`'s point set onto `FIG_S3_S5[grid]`."""
    ref = FIG_S3_S5[grid]
    ref_set = set(ref)
    ref_index = {pt: i + 1 for i, pt in enumerate(ref)}
    plist = [(p.q, p.r) for p in skel.points]
    seen: set[tuple] = set()
    out: list[dict[int, int]] = []
    for img in _ordered_isometry_images(plist):
        ax, ay = img[0]
        for rx, ry in ref:
            dq, dr = rx - ax, ry - ay
            shifted = [(q + dq, r + dr) for q, r in img]
            if set(shifted) == ref_set:
                mapping = {i + 1: ref_index[shifted[i]] for i in range(len(img))}
                key = tuple(sorted(mapping.items()))
                if key not in seen:
                    seen.add(key)
                    out.append(mapping)
    return out


def reference_labeling(skel: Skeleton) -> dict[int, int] | None:
    """Map `skel`'s 1-based sequence label -> the Fig-S3 reference label.

    A survivor's grid may have automorphisms → several isometries land its
    point set on the reference, giving several valid labelings. Because a
    triple's handedness coplanarity depends on the survivor's *sequence-parity*
    pattern (not just positions), these labelings are NOT interchangeable for
    handedness: some transport a paper-mandatory triple onto a coplanar
    survivor triple. We therefore pick, deterministically:

      1. the labeling maximizing the number of paper-mandatory handedness
         triples that are non-coplanar in `skel` (so paper handedness applies);
      2. lexicographically smallest label map among those tied.

    The same labeling is used for interactions (automorphism-equivalent, so any
    is valid there). Returns None if the skeleton matches no Fig-S3 grid. For a
    residual set of survivors (S5: 28 — 21 f, 7 g/h) no labeling makes *all*
    mandatory triples non-coplanar; the max-clean one is returned and the
    coplanar mandatory triples are simply absent from that SSP's handedness
    (a geometric fact, flagged by `mandatory_handedness_gap`).
    """
    grid = grid_of(skel)
    if grid is None:
        return None
    from .combine import handedness_signature
    from itertools import combinations
    tri_index = {t: i for i, t in enumerate(combinations(range(1, skel.dim + 1), 3))}
    sig = handedness_signature(skel)
    mandatory = HANDEDNESS_S5_MANDATORY.get(grid, ())

    def n_noncoplanar(mapping: dict[int, int]) -> int:
        inv = {v: k for k, v in mapping.items()}         # ref label -> skel label
        c = 0
        for (a, b, ck) in mandatory:
            skt = tuple(sorted((inv[a], inv[b], inv[ck])))
            if sig[tri_index[skt]] != 0:
                c += 1
        return c

    best: dict[int, int] | None = None
    best_key: tuple | None = None
    for mapping in _valid_labelings(skel, grid):
        lex = tuple(mapping[i] for i in range(1, skel.dim + 1))
        key = (-n_noncoplanar(mapping), lex)             # max non-coplanar, then lex-min
        if best_key is None or key < best_key:
            best_key, best = key, mapping
    return best


def paper_query_spec(skel: Skeleton, types: tuple[str, ...]) -> dict | None:
    """Assemble the paper-faithful directive spec for one (skeleton, typing).

    Transports the reference grid's minimal interactions + minimal handedness
    onto `skel` via `reference_labeling`, in `skel`'s own 1-based labels:

      grid            : Fig-S3 grid name.
      required_edges  : frozenset{i,j} that must carry a concrete interaction.
      optional_edges  : frozenset{i,j} lattice-adjacent but optional (-> 'X').
      disjunction     : (edgeA, edgeB) for grid e's ">=1 mandatory" pair, or None.
                        Callers emit two query variants (A required / B required).
      hand_triples    : sorted (i,j,k) to emit an L/R handedness line for
                        (mandatory, plus conditional when the typing condition
                        holds), restricted to triples non-coplanar in `skel`.
      hand_gap        : mandatory triples dropped because coplanar in `skel`.

    Returns None if the skeleton matches no Fig-S3 grid.
    """
    grid = grid_of(skel)
    phi = reference_labeling(skel)
    if grid is None or phi is None:
        return None
    inv = {v: k for k, v in phi.items()}                 # ref label -> skel label

    def edge(a: int, b: int) -> frozenset:
        return frozenset((inv[a], inv[b]))

    required = {edge(a, b) for a, b in INTERACTIONS_S5_REQUIRED.get(grid, ())}
    optional = {edge(a, b) for a, b in INTERACTIONS_S5_OPTIONAL.get(grid, ())}
    disj_ref = INTERACTIONS_S5_DISJUNCTION.get(grid)
    disjunction = (edge(*disj_ref[0]), edge(*disj_ref[1])) if disj_ref else None

    # handedness: mandatory always; conditional when the ref-labeled typing
    # meets the grid's condition (f: 1-2-3 helix, 4-5 strand-sheet -> HHHEE).
    ref_types = tuple(types[inv[a] - 1] for a in range(1, skel.dim + 1))
    triples = list(HANDEDNESS_S5_MANDATORY.get(grid, ()))
    cond = HANDEDNESS_S5_CONDITIONAL.get(grid)
    if cond and grid == "f" and ref_types[:3] == ("H", "H", "H") and ref_types[3:] == ("E", "E"):
        triples += list(cond[0])

    from .combine import handedness_signature
    from itertools import combinations
    tri_index = {t: i for i, t in enumerate(combinations(range(1, skel.dim + 1), 3))}
    sig = handedness_signature(skel)
    hand_triples: list[tuple[int, int, int]] = []
    hand_gap: list[tuple[int, int, int]] = []
    for (a, b, c) in triples:
        skt = tuple(sorted((inv[a], inv[b], inv[c])))
        (hand_triples if sig[tri_index[skt]] != 0 else hand_gap).append(skt)

    return {
        "grid": grid,
        "required_edges": required,
        "optional_edges": optional,
        "disjunction": disjunction,
        "hand_triples": sorted(hand_triples),
        "hand_gap": sorted(hand_gap),
    }


def mandatory_handedness_gap(skel: Skeleton) -> list[tuple[int, int, int]]:
    """Paper-mandatory handedness triples (in `skel`'s labels) that are
    geometrically coplanar in `skel` under `reference_labeling` — i.e. cannot
    carry an L/R constraint. Empty for all but the residual degenerate set."""
    grid = grid_of(skel)
    phi = reference_labeling(skel)
    if grid is None or phi is None:
        return []
    from .combine import handedness_signature
    from itertools import combinations
    tri_index = {t: i for i, t in enumerate(combinations(range(1, skel.dim + 1), 3))}
    sig = handedness_signature(skel)
    inv = {v: k for k, v in phi.items()}
    gaps = []
    for (a, b, c) in HANDEDNESS_S5_MANDATORY.get(grid, ()):
        skt = tuple(sorted((inv[a], inv[b], inv[c])))
        if sig[tri_index[skt]] == 0:
            gaps.append(skt)
    return gaps
