# Design-target ProSMoS queries generated from CG-2012 oracle

These 9 `.query` files were auto-generated from `~/chalam/CG-2012/S5/IA.txt`
via `enum/scripts/generate_design_target_queries.py` (using
`enum/src/ssp_enum/prosmos.py`).

These supersede the hand-translated `../queries/` for the panels that exist
in CG-2012: where the oracle and the website BMP diverge, the oracle
is authoritative (it's the byte-source CG-2012 emitted; the website was
a later rendering that doesn't match the underlying enumeration for 13
of 14 panels — see `../queries/README.md` for the diagnosis).

## Coverage

| Panel | Generated | Source |
|---:|:---:|---|
| 1–8 (group i, 5-strand β-sheets) | ✓ all 8 | CG-2012 IA-S5.txt |
| 36, 40, 41 (group iii) | ✗ | post-2012 website enumeration; not in CG-2012 |
| 56, 57 (group vi) | ✗ | post-2012 website enumeration |
| 58 (`5-141-7-7`) | ✓ | byte-identical between CG-2012 and website |

For the 5 missing panels we'd need either:
- The post-2012 enumeration source (not currently available locally), or
- Hand-construction from each BMP using `Skeleton.edges` + ProSMoS writer,
  then validating against the BMP visually.

## How to regenerate

```
cd ~/dev/prosmos_cl/enum
PYTHONPATH=src python3 scripts/generate_design_target_queries.py
```

The generator reads from `enum/reference/IA-S5.txt` (a symlink into
`~/chalam/CG-2012/S5/`), so this only works on a machine where the
chalam archive is mounted. On a clean checkout these query files are
the durable artifact.

## Format conventions

Per `enum/src/ssp_enum/prosmos.py`:
- Matrix uses paper §1.1.1 symbol convention: `c`/`t`/`u`/`v`/`C`/`T` for
  the six lattice-adj interaction types, `X` for optional non-adj,
  `-` for explicit non-interaction.
- One `length` directive per SSE (lab convention — every SSE gets a
  length constraint, not just helices). Defaults: E=5..1000, H=8..1000.
- `sheetS`/`sheetD` directives reproduce the oracle's `same_sheet` and
  `diff_sheet` membership lists.
- `handedness i j k L|R` directives reproduce the oracle's per-triple
  handedness specifications.
