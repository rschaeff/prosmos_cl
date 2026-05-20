"""ProSMoS query.txt writer.

The ProSMoS search tool consumes query files in this format:

    1 2 3 4 5
    H E E E E
    * T C - -
      * t c -
        * - c
          * -
            *
    sheetS 2 3 4 5
    length 1 H 8 1000
    length 2 E 5 1000
    length 3 E 5 1000
    length 4 E 5 1000
    length 5 E 5 1000
    handedness 1 2 3 R

  Line 1: SSE indices 1..n (space-separated, 1-based).
  Line 2: SSE type per index — 'E' (strand), 'H' (helix), or 'X' (any).
  Lines 3..(n+2): upper-triangular interaction matrix, row i has 2*i
    leading spaces then cells separated by single spaces. Cells:
      '*'  diagonal
      '-'  no interaction (also marks non-adj pair of strands in same sheet)
      'X'  optional interaction (non-lattice-adj generally)
      'c'  parallel hydrogen-bond
      't'  anti-parallel hydrogen-bond
      'u'  parallel non-H-bond
      'v'  anti-parallel non-H-bond
      'C'  parallel helix-strand contact
      'T'  anti-parallel helix-strand contact
  Then directives (in any order, repeatable):
      sheetS i j k...   SSEs i, j, k... are in the same sheet
      sheetD i j        SSEs i, j explicitly in different sheets
      length i T min max   SSE i (of type T) has min..max residues
      handedness i j k L|R  triple has the specified chirality

This writer consumes the `SSPRecord` produced by `oracle.parse()` —
that record already has every field the format needs.

Per lab convention (see [[feedback-prosmos-query-length]] in session
memory), every SSE in the query gets a length constraint, not just
those where the minimum is load-bearing. Defaults: β-strand min 5
residues, α-helix min 8 (per paper Methods); max 1000.
"""

from __future__ import annotations

from typing import Iterable, Mapping

from .oracle import SSPRecord


# Paper Methods minimum lengths.
DEFAULT_MIN_LENGTH: Mapping[str, int] = {'E': 5, 'H': 8}
DEFAULT_MAX_LENGTH: int = 1000


def write_query(
    record: SSPRecord,
    *,
    min_lengths: Mapping[str, int] | None = None,
    max_length: int = DEFAULT_MAX_LENGTH,
) -> str:
    """Render `record` in ProSMoS query.txt format.

    Returns the file contents as a string (no trailing newline).
    Caller writes to disk via Path.write_text or similar.

    `min_lengths` overrides per-type minimums (defaults: E=5, H=8).
    `max_length` is the same for every SSE.
    """
    if min_lengths is None:
        min_lengths = DEFAULT_MIN_LENGTH

    n = record.dim
    lines: list[str] = []

    # Line 1: SSE indices.
    lines.append(" ".join(str(i + 1) for i in range(n)))

    # Line 2: SSE types.
    lines.append(" ".join(record.sse_types))

    # Lines 3..(n+2): upper-triangular matrix.
    for i in range(n):
        leading = " " * (2 * i)
        cells = []
        for j in range(i, n):
            if j == i:
                cells.append("*")
            else:
                cell = record.matrix[i][j]
                # If the cell is empty in the parsed matrix (lower triangle
                # of source IA.txt), default to '-' (shouldn't happen for j > i
                # in a well-formed record, but defensive).
                cells.append(cell if cell else "-")
        lines.append(leading + " ".join(cells))

    # Directives.
    # sheetS / sheetD before length, per the example convention.
    for sheet in record.same_sheet:
        lines.append("sheetS " + " ".join(str(x) for x in sheet))
    for sheet in record.diff_sheet:
        lines.append("sheetD " + " ".join(str(x) for x in sheet))

    # length per SSE.
    for idx, sse_type in enumerate(record.sse_types, start=1):
        # Default to 5 for unknown types (e.g., 'Z' in the source).
        # 'Z' shouldn't appear in well-formed records; if it does, the
        # downstream search will reject it anyway.
        type_for_query = sse_type if sse_type in min_lengths else "E"
        mn = min_lengths.get(sse_type, min_lengths.get("E", 5))
        lines.append(f"length {idx} {type_for_query} {mn} {max_length}")

    # handedness.
    for i, j, k, lr in record.handedness:
        lines.append(f"handedness {i} {j} {k} {lr}")

    return "\n".join(lines) + "\n"


def find_record(records: Iterable[SSPRecord], full_id: str) -> SSPRecord:
    """Look up a record by its full id like '5-141-7-7'.

    `full_id` parses as `<dim>-<skel>-<third>-<sub_first>`. The
    sub_second is ignored (it's not in the design-target filename
    convention; the design-target panels are `5-141-7-7.bmp` where
    the trailing `-7` is sub_first and the implicit sub_second is
    whatever the record has).
    """
    parts = full_id.split("-")
    if len(parts) < 4:
        raise ValueError(f"expected dim-skel-third-sub, got {full_id!r}")
    dim, skel, third, sub_first = (int(p) for p in parts[:4])
    for r in records:
        if (r.dim, r.skeleton_id, r.third_idx, r.sub_first) == (
            dim, skel, third, sub_first,
        ):
            return r
    raise KeyError(f"no record matches {full_id!r}")
