#!/usr/bin/env python3
"""Targeted S6 queries: every one-SSE extension of a given S5 cell.

WHY TARGETED. The full S6 grid is 2,372 skeletons x 2**6 typings = 151,808
queries -- 24x the S5 grid, a multi-day sweep. Seeding from a specific S5 cell
and enumerating only the S6 motifs that CONTAIN it costs ~29 queries (about 15
skeleton extensions x 2 for the new SSE being H or E): a ~5,000x reduction.

WHY IT IS ALSO BETTER SCIENCE. The ~15 extensions of one seed form a SIBLING
SET: same parent motif, same sampling depth, differing only in where the sixth
SSE is placed and what type it is. That is an internally controlled comparison.
Absolute emptiness at S6 means little -- every added SSE multiplies the
constraint count, so hit counts fall mechanically and *some* S6 cell is always
empty. The interpretable question is conditional: does THIS one-SSE addition
destroy realizability when its siblings survive? Report occupancy relative to
the sibling set, never as a bare "0 hits => dark" claim.

SEED NAMING mirrors the S5 grid (`s5-<skel>-<typing>`), where skel indexes
`enumerate_skeletons(5)` and typing indexes `product(('H','E'), repeat=5)`
(H=0, E=1, MSB = SSE 1). So `171:7` is sk171 ty07 = HHEEE -- the cell holding
the Ruczinski topology-23 five-SSE motif, which has 0 experimental instances.

MATRIX CONVENTION. Default is `--matrix typed`: the fully typed graph198 encoding
(`v C T t c u` + handedness), so the seed's 5x5 submatrix matches the CANONICAL
graph198 S5 cell cell-for-cell (verified: the seed block of a typed s6t-0171-0007
query equals s5-0171-0007.query exactly) and the S6 extension is comparable to the
typed S5 baseline that actually identified sk171 ty07 as the AFDB-only survivor.
The whole project's canonical grid is graph198; the orientation-blind `adjonly`
grid is a DIFFERENT question (`queries_adjonly`, matches ~20x more) and its sweeps
were quarantined as the wrong query set. `--matrix adjonly` remains available for a
deliberately looser, seed-matched sibling comparison, but it is NOT the canonical
baseline and must never be mixed with the typed grids.

CONTAINMENT IS NOT AUTOMATIC. A one-node extension does not always preserve the
parent's sheet topology: adding a strand can re-partition the sheets so the
seed's own strands end up split across two. Observed on s5-0171-0007, where two
prepend variants turned `sheetS 3 4 5` into `sheetS 4 6` + `sheetS 1 5` -- every
adjacency still matched, but the seed's three-strand sheet was broken, so the S6
motif no longer CONTAINS the S5 one. Such candidates are rejected: `contains_seed`
requires the seed's adjacency submatrix to match cell-for-cell AND each seed sheet
to remain a subset of some S6 sheet (sheets may grow, never split).

JOIN SIDE MATTERS. `combine_with_single_node` both appends and prepends, so the
seed's five SSEs land at positions 1-5 (append) or 2-6 (prepend). The emitted
typing pins the seed's types at whichever positions it actually occupies and
varies only the new SSE, so a query always contains the seed motif intact.

Output:
  <out>/s6t-<skel>-<typing>-<NNN>.query   the queries
  <out>/manifest.tsv                      query -> seed, join side, new SSE
                                          position/type, sibling group
"""
from __future__ import annotations

import argparse
import sys
from itertools import product
from pathlib import Path

_THIS = Path(__file__).resolve().parent
sys.path.insert(0, str(_THIS.parent / "src"))

from ssp_enum.assignment import skeletons_to_records  # noqa: E402
from ssp_enum.combine import canonical_key, combine_with_single_node  # noqa: E402
from ssp_enum.compactness import is_compact  # noqa: E402
from ssp_enum.enumerate import enumerate_skeletons  # noqa: E402
from ssp_enum.prosmos import write_query  # noqa: E402
from ssp_enum.skeleton import Skeleton  # noqa: E402

# Specific interaction codes collapse to the adjacency wildcard; '*', '-' and
# 'X' already carry the same meaning in both conventions.
_TYPED_CODES = frozenset("uvctCT")


def to_adjonly(rec):
    """Rewrite a typed SSPRecord into the adjacency-only convention."""
    import dataclasses
    m = [[("x" if c in _TYPED_CODES else c) for c in row] for row in rec.matrix]
    return dataclasses.replace(
        rec,
        matrix=tuple(tuple(r) for r in m),
        handedness=(),           # adjonly grid carries no handedness directives
    )


def contains_seed(rec, seed_rec, seed_pos) -> bool:
    """True iff `rec` genuinely contains `seed_rec` at positions `seed_pos`.

    Two conditions, both necessary (see CONTAINMENT IS NOT AUTOMATIC above):
      1. every seed pair's interaction code is reproduced exactly, and
      2. every seed sheet survives as a subset of some sheet in `rec` --
         a sheet may absorb the new SSE, but must not be split.
    """
    n = len(seed_pos)
    for i in range(n):
        for j in range(i, n):
            if rec.matrix[seed_pos[i]][seed_pos[j]] != seed_rec.matrix[i][j]:
                return False
    remap = {i + 1: seed_pos[i] + 1 for i in range(n)}     # 1-based seed -> 1-based rec
    rec_sheets = [set(sh) for sh in rec.same_sheet]
    for sheet in seed_rec.same_sheet:
        want = {remap[x] for x in sheet}
        if not any(want <= have for have in rec_sheets):
            return False
    return True


def seed_types(typing_idx: int, dim: int = 5) -> tuple[str, ...]:
    """Typing index -> SSE type tuple, matching the S5 grid's convention."""
    return list(product(("H", "E"), repeat=dim))[typing_idx]


def extensions(seed: Skeleton):
    """Yield (s6_skeleton, join_side, seed_positions) for each compact
    one-node extension. seed_positions are 0-based indices in the S6 skeleton
    that correspond to the seed's SSEs, in order."""
    n = seed.dim
    for cand in combine_with_single_node(seed):
        cand = Skeleton(points=cand.points, chirality=cand.chirality)
        if not is_compact(cand):
            continue
        # combine_with_single_node emits append then prepend, preserving the
        # seed's point order in both cases.
        if cand.points[:n] == seed.points:
            yield cand, "append", tuple(range(n))
        elif cand.points[1:] == seed.points:
            yield cand, "prepend", tuple(range(1, n + 1))


def build(seed_skel_id: int, seed_typing: int, out: Path, dry_run: bool = False,
          matrix: str = "adjonly"):
    s5 = enumerate_skeletons(5)
    if not 0 <= seed_skel_id < len(s5):
        raise SystemExit(f"skeleton {seed_skel_id} out of range (|S5|={len(s5)})")
    seed = s5[seed_skel_id]
    stypes = seed_types(seed_typing)
    tag = f"{seed_skel_id:04d}-{seed_typing:04d}"
    print(f"seed s5-{tag}  types={''.join(stypes)}")

    # the seed's own record, in the same convention, as the containment reference
    seed_rec = next(r for r in skeletons_to_records([seed]) if r.sse_types == stypes)
    if matrix == "adjonly":
        seed_rec = to_adjonly(seed_rec)

    rows, emitted, seen, rejected = [], 0, set(), 0
    for cand, side, seed_pos in extensions(seed):
        ckey = canonical_key(cand)
        # One record per typing of this S6 skeleton; keep only those whose
        # seed positions carry the seed's types (the new SSE is free).
        for rec in skeletons_to_records([cand]):
            if tuple(rec.sse_types[p] for p in seed_pos) != stypes:
                continue
            new_pos = ({*range(cand.dim)} - {*seed_pos}).pop()
            key = (ckey, rec.sse_types)
            if key in seen:                    # same skeleton+typing via both joins
                continue
            out_rec = to_adjonly(rec) if matrix == "adjonly" else rec
            if not contains_seed(out_rec, seed_rec, seed_pos):
                rejected += 1                  # extension broke the seed's sheet
                continue
            seen.add(key)
            name = f"s6t-{tag}-{emitted:03d}"
            if not dry_run:
                (out / f"{name}.query").write_text(write_query(out_rec))
            rows.append((name, f"s5-{tag}", "".join(stypes), side,
                         new_pos + 1, rec.sse_types[new_pos],
                         "".join(rec.sse_types), str(ckey)))
            emitted += 1

    if not dry_run:
        out.mkdir(parents=True, exist_ok=True)
        with open(out / "manifest.tsv", "w") as fh:
            fh.write("query\tseed_cell\tseed_types\tjoin_side\tnew_sse_pos\t"
                     "new_sse_type\ts6_types\ts6_canonical_key\n")
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")
    print(f"  {len(seen)} distinct (skeleton, typing) -> {emitted} queries"
          f"  [{rejected} rejected: extension broke the seed's sheet]")
    return emitted


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("seeds", nargs="+", metavar="SKEL:TYPING",
                    help="S5 seed cells, e.g. 171:7 (sk171 ty07)")
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--matrix", choices=("adjonly", "typed"), default="typed",
                    help="typed (default) = canonical graph198 encoding; adjonly is "
                         "the orientation-blind grid (NOT canonical, do not mix)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if not a.dry_run:
        a.out.mkdir(parents=True, exist_ok=True)
    total = 0
    for s in a.seeds:
        sk, ty = (int(x) for x in s.split(":"))
        total += build(sk, ty, a.out, a.dry_run, a.matrix)
    print(f"total: {total} queries -> {a.out}")


if __name__ == "__main__":
    main()
