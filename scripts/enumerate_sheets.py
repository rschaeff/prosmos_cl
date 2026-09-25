#!/usr/bin/env python3
"""Enumerate n-strand beta-sheet topologies and write them as ProSMoS queries.

The four-strand census used 96 topologies = 12 strand orders x 8 orientations. That
factorization is the quotient of all n! * 2^n (order, orientation) assignments by a group
of order 4:

  reversal  -- read the sheet from the other end:  spatial_seq -> spatial_seq[::-1]
  flip      -- invert every strand's direction:    orientation -> complement

In this encoding the two act on DIFFERENT coordinates and so are independent, which is why
the count factors exactly: (n!/2) * (2^n/2). Canonical representative, matching
ruczinski/queries/motifs.tsv: spatial_seq is the lexicographic min of itself and its
reverse, and strand 1 points up.

ENCODING TRAP, verified against all 96 published rows (2026-09-24):
`orientation` is indexed by SEQUENCE STRAND -- orientation[i] is the direction of strand
i+1 -- NOT by spatial slot. `ruczinski/verify_panels.py:feats()` uses a per-slot
`ori_drawn` instead, which is the drawn-panel convention (the `panel_ori_drawn` column) and
a different thing. Indexing by slot reproduces only 42 of the 96 code matrices.

Two strands interact iff they are adjacent in the sheet (|position difference| == 1); the
code is `c` when they point the same way (parallel) and `t` when they do not
(antiparallel). Every other pair is `-`. The resulting pairing graph is a path, so interior
strands have exactly two c/t neighbours -- which is the definition acceptance.py applies.

Self-checks (all asserted at run time):
  * dedup by code matrix and dedup by the group action agree
  * n=4 reproduces exactly the 96 code matrices in the published motifs.tsv
  * counts: 96 at n=4, 960 at n=5, 11520 at n=6
"""

from __future__ import annotations

import argparse
import csv
import os
from itertools import permutations, product

MIN_LEN_E = 5
MAX_LEN = 1000


def pairs(n):
    """Upper-triangular pair order: (1,2),(1,3),...,(1,n),(2,3),... -- as motifs.tsv."""
    return [(a, b) for a in range(1, n + 1) for b in range(a + 1, n + 1)]


def code_matrix(spatial_seq, orientation, n):
    """Ten characters at n=5, six at n=4. orientation is indexed by sequence strand."""
    pos = {int(spatial_seq[p]): p for p in range(n)}
    ori = {i + 1: orientation[i] for i in range(n)}
    return "".join(
        ("c" if ori[a] == ori[b] else "t") if abs(pos[a] - pos[b]) == 1 else "-"
        for a, b in pairs(n)
    )


def jumps(spatial_seq, n):
    """Consecutive strands in sequence that are not adjacent in the sheet."""
    pos = {int(spatial_seq[p]): p for p in range(n)}
    return sum(1 for i in range(1, n) if abs(pos[i] - pos[i + 1]) != 1)


def enumerate_topologies(n):
    """Canonical representatives, in a stable order."""
    orders = sorted({min(s, s[::-1]) for s in ("".join(map(str, p)) for p in permutations(range(1, n + 1)))})
    orients = ["u" + "".join(o) for o in product("ud", repeat=n - 1)]
    return [(s, o) for s in orders for o in sorted(orients)]


def check_quotient(n, topologies):
    """The group action and the code matrix must induce the same equivalence."""
    by_group = {
        (min(s, s[::-1]), o if o[0] == "u" else "".join("u" if c == "d" else "d" for c in o))
        for s in ("".join(map(str, p)) for p in permutations(range(1, n + 1)))
        for o in ("".join(x) for x in product("ud", repeat=n))
    }
    by_code = {code_matrix(s, o, n) for s, o in by_group}
    assert len(by_group) == len(topologies), (len(by_group), len(topologies))
    assert len(by_code) == len(topologies), (
        f"code matrix is not a complete invariant at n={n}: "
        f"{len(by_code)} matrices for {len(topologies)} topologies"
    )
    expected = (len(list(permutations(range(n)))) // 2) * (2 ** n // 2)
    assert len(topologies) == expected, (len(topologies), expected)


def write_query(spatial_seq, orientation, n):
    code = code_matrix(spatial_seq, orientation, n)
    it = iter(code)
    lines = [
        " ".join(str(i + 1) for i in range(n)),
        " ".join("E" for _ in range(n)),
    ]
    grid = {}
    for a, b in pairs(n):
        grid[(a, b)] = next(it)
    for i in range(1, n + 1):
        cells = ["*"] + [grid[(i, j)] for j in range(i + 1, n + 1)]
        lines.append(" " * (2 * (i - 1)) + " ".join(cells))
    lines.append("sheetS " + " ".join(str(i + 1) for i in range(n)))
    for i in range(1, n + 1):
        lines.append(f"length {i} E {MIN_LEN_E} {MAX_LEN}")
    return "\n".join(lines) + "\n"


def validate_against_published(n, topologies, tsv, qdir):
    """n=4 must reproduce the published set exactly -- matrices and query text."""
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    pub_codes = {r["code_matrix"] for r in rows}
    ours = {code_matrix(s, o, n) for s, o in topologies}
    assert ours == pub_codes, (
        f"code matrix sets differ: {len(ours - pub_codes)} ours-only, "
        f"{len(pub_codes - ours)} published-only"
    )
    by_code = {r["code_matrix"]: r["motif_id"] for r in rows}
    mismatched = []
    for s, o in topologies:
        c = code_matrix(s, o, n)
        path = os.path.join(qdir, "motif_%02d.query" % int(by_code[c]))
        if open(path).read() != write_query(s, o, n):
            mismatched.append(path)
    assert not mismatched, f"{len(mismatched)} query files differ, e.g. {mismatched[:3]}"
    for r in rows:
        s, o = r["spatial_seq"], r["orientation"]
        assert code_matrix(s, o, n) == r["code_matrix"], r["motif_id"]
        assert str(jumps(s, n)) == r["jumps"], r["motif_id"]
    return len(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n", type=int, default=5, help="strands per sheet (default 5)")
    p.add_argument("--out", default=".",
                   help="output root; writes motifs<n>.tsv and queries/ (default: cwd)")
    p.add_argument("--validate-n4", metavar="RUCZINSKI_DIR",
                   default="/home/rschaeff/work/prosmos_2026/ruczinski",
                   help="published four-strand set to validate the construction against")
    p.add_argument("--dry-run", action="store_true", help="check and count, write nothing")
    a = p.parse_args()

    tops = enumerate_topologies(a.n)
    check_quotient(a.n, tops)
    print(f"n={a.n}: {len(tops)} topologies "
          f"({len({s for s, _ in tops})} strand orders x {len({o for _, o in tops})} orientations)")

    if a.validate_n4:
        n4 = enumerate_topologies(4)
        check_quotient(4, n4)
        k = validate_against_published(
            4, n4,
            os.path.join(a.validate_n4, "queries", "motifs.tsv"),
            os.path.join(a.validate_n4, "queries"),
        )
        print(f"validated: reproduces all {k} published four-strand topologies, "
              f"code matrices and query text byte-for-byte")

    if a.dry_run:
        return

    qdir = os.path.join(a.out, "queries")
    os.makedirs(qdir, exist_ok=True)
    tsv = os.path.join(a.out, f"motifs{a.n}.tsv")
    with open(tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["motif_id", "spatial_seq", "orientation", "code_matrix", "jumps"])
        for i, (s, o) in enumerate(tops, start=1):
            w.writerow([i, s, o, code_matrix(s, o, a.n), jumps(s, a.n)])
            with open(os.path.join(qdir, "motif_%04d.query" % i), "w") as q:
                q.write(write_query(s, o, a.n))
    print(f"wrote {tsv} and {len(tops)} queries to {qdir}")


if __name__ == "__main__":
    main()
