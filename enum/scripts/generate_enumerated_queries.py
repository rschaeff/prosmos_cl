"""Generate ProSMoS query.txt files for the full enumerated S3-S5 SSPs.

For each dimension N in {3, 4, 5}, this runs `enumerate_skeletons(N)` to get
the hex-induced skeleton set (sizes 5, 42, 198 respectively), then crosses
each skeleton with the full 2**N Cartesian product of `H`/`E` SSE-type
assignments via `assignment.skeletons_to_records`, and writes each record to
`example/ssp_enumerated/queries_typed/sN/sN-<skel>-<type>.query` using the
`prosmos.write_query` writer.

S3 alone is 11 SSPs when you include the L/R/None chirality variants, but the
combine pipeline (paper Methods, §1.1) operates on planar shapes only; the
chirality split is documented as added downstream and is captured here only
through `handedness_signature`-derived per-triple directives.

Total emitted at S3-S5:
  S3: 5 * 8   =     40
  S4: 42 * 16 =    672
  S5: 198 * 32 =  6336
  Total                7048
"""

from __future__ import annotations

import argparse
from pathlib import Path

import sys

# When run directly (not `python -m`), add src/ to the path.
_THIS = Path(__file__).resolve().parent
_ENUM_ROOT = _THIS.parent
sys.path.insert(0, str(_ENUM_ROOT / "src"))

from ssp_enum.assignment import skeletons_to_records  # noqa: E402
from ssp_enum.enumerate import enumerate_skeletons  # noqa: E402
from ssp_enum.prosmos import write_query  # noqa: E402


DEFAULT_OUT = _ENUM_ROOT.parent / "example" / "ssp_enumerated" / "queries_typed"


def emit(out_root: Path, dims: list[int], dry_run: bool = False) -> None:
    out_root.mkdir(parents=True, exist_ok=True)
    grand_total = 0
    for n in dims:
        sub = out_root / f"s{n}"
        sub.mkdir(parents=True, exist_ok=True)
        skels = enumerate_skeletons(n)
        count = 0
        for rec in skeletons_to_records(skels):
            # Stable filename: dimension, skeleton index, type index.
            fname = f"s{n}-{rec.skeleton_id:04d}-{rec.sub_first:04d}.query"
            if not dry_run:
                (sub / fname).write_text(write_query(rec))
            count += 1
        grand_total += count
        print(f"S{n}: {len(skels)} skeletons -> {count} queries -> {sub}")
    print(f"total: {grand_total} queries")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--out", type=Path, default=DEFAULT_OUT,
        help=f"output root (default: {DEFAULT_OUT})",
    )
    p.add_argument(
        "--dims", type=int, nargs="+", default=[3, 4, 5],
        help="dimensions to emit (default: 3 4 5)",
    )
    p.add_argument(
        "--dry-run", action="store_true",
        help="count only, don't write files",
    )
    args = p.parse_args()
    emit(args.out, args.dims, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
