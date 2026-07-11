# Phase 3 result — paper-faithful S5 query generation

Wires the congruence-labeling map + the per-grid minimal handedness / minimal
interactions (Steps A/B) into query generation. `skeletons_to_records(...,
paper_faithful=True)` now emits **paper-faithful S5 SSP records** from the
geometric-SCC-2 skeleton set.

## The labeling map

`geom_scc2.reference_labeling(skel)` maps a survivor's 1-based sequence labels
onto its Fig-S3 grid's reference labels, via the hex isometry + translation that
carries the point set onto the reference. Grids have automorphisms (multiple
valid labelings); because a triple's handedness coplanarity depends on the
survivor's **sequence-parity** pattern (not just positions), those labelings are
not interchangeable for handedness. The map therefore picks, deterministically:
1. the labeling maximizing the count of paper-mandatory handedness triples that
   are **non-coplanar** in the survivor; then
2. lexicographically smallest.

Validated: bijection for all 140 survivors; a fully clean labeling exists for
**112/140**.

## Output

`paper_faithful=True` (implies the typing filter) → **708 records** from the 140
survivors:

| grid | survivors | valid typings/skel | variants | records |
|---|---|---|---|---|
| d | 12 | 1 (EEEEE) | 1 | 12 |
| e | 40 | 2 (run-4 → E; node 5 free) | **2** | 160 |
| g/h | 28 | 2 (two crossing 3-runs → all-E/all-H) | 1 | 56 |
| f | 60 | 8 (3-run E/H × 4,5 free) | 1 | 480 |

Per record:
- **interactions** — required grid edges get a concrete code, optional lattice
  edges get `X`, non-lattice same-sheet strand pairs get `-`.
- **grid e disjunction** — the `≥1 of {(2,5),(3,5)}` broken pair is emitted as
  **two variants** (`sub_second` 0/1): each requires one edge and X's the other.
  Union of the two variants' hits = the paper's "at least one" semantics.
- **handedness** — mandatory triples always; the **f conditional** pair
  `{(1,4,5),(3,4,5)}` added when the reference-labeled typing is `HHHEE`
  (1-2-3 helix, 4-5 strand-sheet); each with its geometric L/R sign.

## Invariants (tested)

- No handedness line is ever emitted on a coplanar (sign-0) triple.
- Every matrix is well-formed (diagonal `*`, full upper triangle).
- Grid-e variants differ only in the disjunction cells; handedness identical.
- f: clean records carry 7 handedness lines, or **9** when the conditional fires.

## Known residual — the 28 handedness-degenerate survivors

For **28 survivors (21 f, 7 g/h)** *no* labeling makes all paper-mandatory
triples non-coplanar: the survivor's sequence-parity pattern forces one mandatory
triple coplanar, so it cannot carry an L/R line. These records come out with
1–2 fewer handedness lines (the 4/5/8 buckets in the count distribution) and are
flagged by `geom_scc2.mandatory_handedness_gap`.

This is almost certainly an **enumeration artifact**, not a paper phenomenon:
the CG-2012 oracle never has a grid-f block with < 7 handedness lines, i.e. the
paper's enumeration does not produce these parity-degenerate labelings. Our
`enumerate.py` explicitly **defers RCC + handedness-equivalence dedup** (see its
docstring); those 28 are the expected fallout and should collapse once that dedup
lands. Tracked for **Phase 5** (oracle reconciliation). Phase 3 emits them
best-effort and flags them rather than dropping silently.

## What's wired vs pending

- typing filter, geometric grids, minimal handedness, minimal interactions,
  labeling map, grid-e variants, f conditional — **wired**.
- Pending: **Phase 4** — write the 708 to ProSMoS `.query` files + the cheap
  impact probe (sieve / occupancy / negspace on a sample, old-198 vs new-140) →
  the go/no-go on a full darkness/negspace re-run. **Phase 5** — RCC/handedness-
  equivalence dedup that should resolve the 28 residual and reconcile the 648.

## Reproduce

```python
from ssp_enum.compactness import set_scc2_mode
from ssp_enum.enumerate import enumerate_skeletons
from ssp_enum.assignment import skeletons_to_records
set_scc2_mode("geometric")
s5 = enumerate_skeletons(5)
set_scc2_mode("graph")
recs = list(skeletons_to_records(s5, paper_faithful=True))   # 708
```
