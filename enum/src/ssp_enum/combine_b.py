"""Model B (combinatorial) enumeration of S_n skeletons.

Parallel to `combine.py` (Model A, geometric). The two models differ in
how "adjacency" is determined for a skeleton:

  Model A (combine.py):
    Adjacency comes from lattice positions — two SSEs are adjacent iff
    their hex points are at unit distance. The full geometric adjacency
    graph of any placement is automatically the skeleton's edge set.

  Model B (this module):
    Adjacency is declared explicitly during combine. Each combine step
    adds one new node and declares one new edge (the "join" edge to the
    existing skeleton). Geometric coincidences — hex-adjacencies that
    arise from the chosen placement but aren't declared — are recorded
    as `X` (optional non-edge) in the resulting ProSMoS query rather
    than as edges.

Why both:
  Model A respects paper Appendix §1.1 strictly (adjacency = hex unit
  distance). It cannot produce certain oracle skeletons (K1,4 star,
  C5 cycle) because those graphs aren't realizable as induced subgraphs
  on the 2D hex lattice. Model B matches CG-2012's enumeration, which
  treats skeletons as abstract adjacency graphs and admits these
  structures. See `project_s5_oracle_gap.md` in session memory.

What this module currently implements (Phase 2 minimum):
  * `enumerate_skeletons_combinatorial(n)` — Sn via single-node combine
    only (multi-node combine for Sn = Sp + Sq with q ≥ 2 is deferred).
  * SCC-1 is relaxed to "the new node is geometrically adjacent to ≥1
    existing node" — sufficient for connectivity. Paper's "≥2 declared
    adj" is too strict to reach sparse oracle topologies via combine.
  * PCC still applied geometrically.
  * Each new node declares one join edge to a chosen existing node.
    Other geometric adjacencies of the new node to existing nodes are
    *not* declared (they'd be marked X in the downstream ProSMoS query).

The combinatorial enumeration emits skeletons with explicit `edges`
sets; the geometric adjacency available via `points` is still recorded
(it's what PCC uses) but is not used for the skeleton's edge graph.
"""

from __future__ import annotations

from typing import Iterator

from .combine import canonical_key, valid_extension_points
from .compactness import passes_pcc, passes_scc_2
from .grids import is_in_whitelist
from .lattice import LatticePoint
from .skeleton import Skeleton


# Seed skeletons with explicit edges.
_SEED_S1_B = Skeleton(points=(LatticePoint(0, 0, 0),), edges=frozenset())
_SEED_S2_B = Skeleton(
    points=(LatticePoint(0, 0, 0), LatticePoint(1, 0, 0)),
    edges=frozenset({(0, 1)}),
)


def _adjacent_existing(new_point: LatticePoint, s1: Skeleton) -> list[int]:
    """Indices of s1 points geometrically adjacent to `new_point`."""
    return [i for i, p in enumerate(s1.points) if new_point.is_adjacent(p)]


def combine_b_single_node(s1: Skeleton) -> Iterator[Skeleton]:
    """Yield Sn+1 candidates by adding one node to s1.

    For each candidate position adjacent to ≥1 existing point, and for
    each choice of which adjacent existing point to declare-edge to,
    emit two candidates (end-join and front-join). Each emission adds
    *one* declared edge.

    SCC-1 relaxation: the paper's "≥2 adjacent existing" requirement
    is replaced with "≥1 adjacent existing" so that K1,4-style sparse
    arrangements are reachable. The geometric `valid_extension_points`
    (which enforces the paper's strict rule) is *not* used here. PCC
    is still enforced via `passes_pcc` on the geometric arrangement.
    """
    s1_set = set(s1.points)
    base_edges = s1.edges if s1.edges is not None else frozenset()

    # Candidate positions: any hex-adj-to-skeleton point not already in skeleton.
    candidates: set[LatticePoint] = set()
    for p in s1.points:
        for nbr in p.neighbors():
            if nbr not in s1_set:
                candidates.add(nbr)

    n = s1.dim
    for new_point in candidates:
        adj_existing = _adjacent_existing(new_point, s1)
        if not adj_existing:
            continue  # not adjacent to anything — disconnected
        for join_idx in adj_existing:
            # End-join: new node gets label n (0-indexed); existing labels unchanged.
            end_points = s1.points + (new_point,)
            end_edges = base_edges | {(join_idx, n)}
            yield Skeleton(points=end_points, edges=frozenset(end_edges))
            # Front-join: new node gets label 0; existing labels shift +1.
            front_points = (new_point,) + s1.points
            shifted = {(a + 1, b + 1) for (a, b) in base_edges}
            # Join edge: (0, join_idx + 1) in new labeling.
            new_join_a, new_join_b = sorted((0, join_idx + 1))
            front_edges = shifted | {(new_join_a, new_join_b)}
            yield Skeleton(points=front_points, edges=frozenset(front_edges))


def _is_b_compact(skel: Skeleton) -> bool:
    """Combinatorial-mode compactness check.

    PCC is geometric and unchanged (paper §1.1: based on the induced
    grid on the hex lattice). SCC-1 is *not* enforced here — the
    relaxed "≥1 adj" rule is built into `combine_b_single_node`'s
    candidate generation, so any skeleton reaching this point already
    satisfies the relaxed connectivity. SCC-2 (whitelist of induced
    grids) is applied against the *geometric* induced grid since the
    paper's grid whitelist describes lattice arrangements, not
    declared-edge graphs.

    Note: applying geometric SCC-2 here means combinatorial skeletons
    whose lattice arrangement doesn't match a paper grid are rejected
    even though their declared-edge graph might match an oracle
    skeleton. This is conservative — relax later if the design-target
    work demands sparser lattice placements.
    """
    if not passes_pcc(skel):
        return False
    # For SCC-2 we check the geometric induced grid (lattice positions),
    # since the whitelist is paper-grid-based and grids are lattice
    # arrangements regardless of declared-edge interpretation.
    geom_view = Skeleton(points=skel.points, chirality=skel.chirality,
                         start_up=skel.start_up, edges=None)
    return is_in_whitelist(geom_view)


def enumerate_skeletons_combinatorial(n: int) -> list[Skeleton]:
    """Sn via Model B (combinatorial) single-node growth.

    Phase 2 scope: single-node only. Multi-node combine in this model
    is more involved (declared-edge merging during join) and is
    deferred until we exercise this against the 14 design targets.
    """
    if n == 1:
        return [_SEED_S1_B]
    if n == 2:
        return [_SEED_S2_B]
    prev = enumerate_skeletons_combinatorial(n - 1)
    seen: dict[tuple, Skeleton] = {}
    for s1 in prev:
        for candidate in combine_b_single_node(s1):
            if not _is_b_compact(candidate):
                continue
            # canonical_key uses adjacency_matrix() which now reflects
            # the *declared* edges (Model B). So two skeletons with the
            # same lattice but different declared graphs get distinct keys.
            key = canonical_key(candidate)
            if key not in seen:
                seen[key] = candidate
    return list(seen.values())
