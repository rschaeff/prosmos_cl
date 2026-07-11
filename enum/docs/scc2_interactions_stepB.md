# S5 (and S4) minimal interaction sets — read off Fig S3

The **minimal interactions** (Rule 3) read directly off Chitturi 2016 **Fig S3**
(SI page 11) line styles, for the geometric-SCC-2 reference labelings. Constants
in `geom_scc2.INTERACTIONS_S5_{REQUIRED,DISJUNCTION,OPTIONAL}`.

Fig-S3 caption: *"Solid lines indicate an interaction. In Fig 5(e) two broken
lines indicate that at least one of the interactions is mandatory. Through
experimentation we determined the minimum number of (specific) interactions for
each grid so that SSPs derived from skeletons (that induce distinct grids) do not
have common hits."*

Line-style → query semantics:
- **solid** = REQUIRED: a concrete interaction demanded (letter code in the SSP).
- **broken** = DISJUNCTION: "at least one mandatory" (grid e's two dashed lines).
- **neither** = OPTIONAL (`X`): lattice-adjacent but not required.

The Fig-S3 node labels coincide with our reference labeling (`FIG_S3_S5`) for
every S5 grid, so the sets below are already in our labeling. **Consistency
check passed**: for each grid, REQUIRED ∪ DISJUNCTION ∪ OPTIONAL = exactly the
grid's full hex lattice-edge set (no leftovers, no overlaps).

## S5 (the analysis target)

| grid | lattice edges | required (solid) | disjunction (≥1) | optional (X) |
|---|---|---|---|---|
| **d** P5 path | 4 | `(1,2)(2,3)(3,4)(4,5)` | — | — |
| **e** collinear+corner | 5 | `(1,2)(2,3)(3,4)` | **`(2,5)‖(3,5)`** | — |
| **f** ring | 7 | `(1,2)(2,3)(3,4)(4,5)(1,5)` | — | `(2,4)(2,5)` |
| **g/h** star | 6 | `(1,2)(2,3)(2,4)(2,5)` | — | `(1,5)(3,4)` |

- **d**: full path is required; nothing optional.
- **e**: the 1-2-3-4 run is required solid; node 5's two attachments `(2,5)`,`(3,5)`
  are the **broken pair** — at least one must hold. No plain-optional edges.
- **f**: required = the **5-cycle ring** `1-2-3-4-5-1`; the two spokes from the
  ring into the "waist" node, `(2,4)` and `(2,5)`, are optional X. (Node 2 in the
  figure connects only to 1 and 3 by solid lines.)
- **g/h**: required = the **K1,4 star** (centre = node 2, the 4 spokes); the two
  leaf-leaf lattice edges (`(1,5)`,`(3,4)` in our reference labeling) are optional
  X. **g and h have the identical required set** — they differ only in leaf
  geometry (which leaf-leaf edges exist), which is exactly why they collapse to
  our single "gh" grid.

## S4 (secondary; for later generalization)

| grid | required (solid) | optional (X) |
|---|---|---|
| **a** P4 path | `(1,2)(2,3)(3,4)` | — |
| **b** triangle+pendant | `(1,2)(2,3)(2,4)(3,4)` | — |
| **c** C4 cycle | `(1,2)(2,3)(3,4)(1,4)` | short rhombus diagonal |

## The broken-line disjunction (grid e) — Phase-3 encoding question

The "≥1 of `(2,5)`,`(3,5)` mandatory" is a **disjunction** the plain ProSMoS
query matrix can't express directly (each cell is one code). Options for Phase 3,
in preference order:
1. Verify whether searchmatrix supports an OR/disjunction convention (unlikely).
2. **Emit two query variants** per e-SSP — one with `(2,5)` required + `(3,5)` X,
   one with `(3,5)` required + `(2,5)` X — and union the hits. Faithful.
3. Fallback: require **both** (a slight over-specification) — documented.

Recommend option 2 (faithful, cheap: ×2 queries only on grid-e SSPs).

## Reproduce

Fig S3 render: `gs -r500 -dFirstPage=11 -dLastPage=11` on
`doc/1-s2.0-S0022283616302893-mmc1 (1).pdf`; grids cropped per-panel. Partition
check: the snippet in this task / `geom_scc2` constants vs
`skeleton.adjacency_matrix()` lattice edges.

## Phase 3 now has all its inputs

- typing rule — **landed** (`assignment.paper_typing`, 744).
- geometric SCC-2 / grids — **landed** (`geom_scc2`, S5=140).
- minimal handedness — **read** (`scc2_handedness_stepB.md`).
- minimal interactions — **read** (this doc).

Remaining for Phase 3 implementation: the **congruence-labeling map** (arbitrary
survivor → reference labels) to apply these per-grid sets, plus the grid-e
two-variant emission. Then Phase 4's impact probe.
