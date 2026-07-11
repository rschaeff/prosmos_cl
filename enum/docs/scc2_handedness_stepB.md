# S5 minimal handedness lists (per Fig-S3 grid, reference labeling)

Reads the per-grid **minimal handedness** sets that Phase 3 will apply, for the
geometric-SCC-2 reference labelings (`geom_scc2.FIG_S3_S5`). Source: the CG-2012
oracle `reference/IA-S5.txt` (which reproduces the paper's per-count
distribution `{0:45, 2:440, 6:1378, 7:826, 9:118}` exactly), cross-checked
against our own geometric handedness signature (`combine.handedness_signature`).

The listed triples carry an L/R constraint; the **L/R value itself is geometric**
(emitted per-skeleton at query-build). A grid's minimal list can be *smaller*
than its full geometric non-coplanar set — the paper minimizes to the smallest
set that keeps SSPs from distinct grids from sharing hits ("the minimum number
of (specific) interactions/handedness ... through experimentation", Fig-S3
caption). Constants live in `geom_scc2.HANDEDNESS_S5_{MANDATORY,CONDITIONAL}`.

## The lists

| grid | edges | geom non-coplanar | paper minimal | list (reference labeling) |
|---|---|---|---|---|
| **d** | 4 | 0 | 0 | — |
| **e** | 5 | 6 | 2 | `(1,4,5) (2,3,5)` * |
| **g/h** | 6 | 6 | 6 | `(1,2,4) (1,2,5) (1,3,4) (1,3,5) (2,3,4) (2,3,5)` |
| **f** | 7 | 9 | 7 (+2 cond) | `(1,2,4) (1,2,5) (1,3,4) (1,3,5) (2,3,4) (2,3,5) (2,4,5)` |

- **g/h** is unambiguous: minimal == the complete geometric non-coplanar set
  (6 = 6). Present in the oracle as a single count-6 set ×186.
- **f** is 7 **mandatory** + 2 **conditional**. The conditional pair
  `{(1,4,5),(3,4,5)}` turns on **only** for the SSE typing **`HHHEE`** (nodes
  1-2-3 helix, nodes 4-5 strand → sheet). Verified decisively in oracle block
  `skel=41` (which carries our reference f labeling): typings
  `EEEHH/EEEEH/EEEHE/EEEEE/HHHHH/HHHEH/HHHHE` → 7 lines; typing `HHHEE` → 9 lines
  (= 7 + `(1,4,5),(3,4,5)`). Matches SI ("1,2,3 helices & 4,5 strands→sheet").
  The full 9-set is our exact geometric non-coplanar set (oracle count-9 ×12).

### The g/h and f derivations (the ask)

```
f  block skel=41 (reference f labeling):
   types EEEHH EEEEH EEEHE EEEEE HHHHH HHHEH HHHHE -> 7 mandatory
   types HHHEE                                     -> 9 (7 + conditional {(1,4,5),(3,4,5)})
g/h  single oracle count-6 set, ×186 = full geometric non-coplanar 6-set
```

## Residual (e only, does not affect f / g-h)

Grid e's minimal-2 has **three** symmetry-equivalent all-node-5 candidates in the
oracle — `{(1,4,5),(2,3,5)}` (×24), `{(1,2,5),(3,4,5)}` (×32),
`{(1,3,5),(2,4,5)}` (×32) — each a valid size-2 disambiguator (all contain the
pendant node 5, all ⊆ our e-6). The exact match to *our* reference e labeling
was not pinned: the oracle stores minimal *interactions*, not the full lattice
grid (lattice-adjacent pairs can be X-marked), so node-5 adjacency can't be read
back cleanly for a geometric block-match. `(1,4,5),(2,3,5)` is recorded as the
provisional pick (consistent with the earlier plan decode). To finalize: match
an oracle e-block's full geometry to `FIG_S3_S5["e"]` via a lattice
reconstruction, or reproduce the paper's disambiguation optimization on our 4
grids. Any of the three is correct for disambiguation.

## Provenance / reproduce

`work/prosmos_2026/handedness_lists_derive.py` (prints the tables + the skel=41
block). Oracle: `enum/reference/IA-S5.txt` → `/home/rschaeff/chalam/CG-2012/S5/IA.txt`.

## Feeds Phase 3

`HANDEDNESS_S5_MANDATORY` / `HANDEDNESS_S5_CONDITIONAL` in `geom_scc2.py`. Phase 3
emits, per skeleton (mapped to its grid via `grid_of` + the congruence
labeling): the mandatory triples always, the conditional triples when the
skeleton's SSE typing meets the condition — each with its geometric L/R sign.
Still pending for Phase 3: the **minimal interaction** sets (Fig-S3 solid =
required, broken = ≥1-mandatory), read off the figure line styles, and the
congruence-labeling map from arbitrary survivor → reference labels.
