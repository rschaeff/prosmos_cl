"""Generate ProSMoS query.txt files for the design targets reachable
from the CG-2012 oracle (9 of 14).

Run from the `enum/` directory:
    PYTHONPATH=src python3 scripts/generate_design_target_queries.py

Output goes to ../example/ssp_design_targets/queries_enum/ (alongside
the existing hand-translated queries/ directory). Files are named to
match the panel naming used elsewhere in the repo.

The 5 unreachable targets (5-283-1-2, 5-307-1-2, 5-243-1-2, 5-265-7-7,
5-234-7-7) are from a later post-2012 enumeration on the website and
have no corresponding CG-2012 IA-S5.txt record; they are skipped.
"""

from __future__ import annotations

from pathlib import Path

from ssp_enum.oracle import parse
from ssp_enum.prosmos import find_record, write_query

# (Panel number, full id) per the design-target index.
TARGETS = [
    (1, "5-269-0-0"),
    (2, "5-311-0-0"),
    (3, "5-289-0-0"),
    (4, "5-288-0-0"),
    (5, "5-280-0-0"),
    (6, "5-282-0-0"),
    (7, "5-306-0-0"),
    (8, "5-309-0-0"),
    (36, "5-283-1-2"),
    (40, "5-307-1-2"),
    (41, "5-243-1-2"),
    (56, "5-265-7-7"),
    (57, "5-234-7-7"),
    (58, "5-141-7-7"),
]

ENUM_DIR = Path(__file__).resolve().parent.parent
REFERENCE = ENUM_DIR / "reference"
OUT_DIR = ENUM_DIR.parent / "example" / "ssp_design_targets" / "queries_enum"


def main() -> None:
    records = list(parse(REFERENCE / "IA-S5.txt"))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    reachable, unreachable = [], []
    for panel, target_id in TARGETS:
        try:
            r = find_record(records, target_id)
        except KeyError:
            unreachable.append((panel, target_id))
            continue
        path = OUT_DIR / f"{panel:02d}-{target_id}.query"
        path.write_text(write_query(r))
        reachable.append((panel, target_id, path))

    print(f"Wrote {len(reachable)} queries to {OUT_DIR}:")
    for panel, target_id, path in reachable:
        print(f"  panel {panel:>2}: {target_id}  →  {path.name}")
    if unreachable:
        print(f"\nSkipped {len(unreachable)} targets not in CG-2012 IA-S5.txt:")
        for panel, target_id in unreachable:
            print(f"  panel {panel:>2}: {target_id}  (post-2012 website enumeration)")


if __name__ == "__main__":
    main()
