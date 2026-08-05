# Plan — geometric SCC-2 (reconcile our enumeration with the paper's grids)

**Goal.** Replace our *graph-based* secondary-compactness check (SCC-2) with a
*geometric* one, so our enumerated skeleton set becomes the paper-faithful hex set
on which the Chitturi 2016 Fig-S3 **minimal handedness** and **minimal interaction**
sets apply cleanly. This closes Rules 2 & 3 (blocked) and resolves the standing
198-vs-648 discrepancy at its root.

## Root cause (established)
- `grids.py:unlabeled_grid_signature` + `WHITELIST_S{3,4,5}` test a candidate's
  induced grid by its **unlabeled adjacency graph**.
- The paper (SI §, after Fig S2): a skeleton passes SCC-2 iff its induced grid
  **"matches (or is a symmetric variation of) one of these predefined grids"** —
  i.e. hex-isometry **geometric congruence** to a Fig-S3 grid, *not* graph iso.
- A graph can have several non-congruent hex realizations. Our graph-check admits
  all of them (e.g. straight grid-e **and** a bent 5-1-3-4 path with a pendant that
  is graph-isomorphic to e but not the Fig-S3 shape). The paper admits only the
  Fig-S3 shape. Those extra geometries have **no paper grid → no paper handedness**,
  which is exactly why ~20% of our 198 don't map (a paper handedness triple becomes
  geometrically coplanar for them). Confirmed: skeleton with edges
  {(1,2),(1,3),(1,5),(2,3),(3,4)} → paper (1,4,5) maps to a coplanar our-triple.

## Assets already in hand
- **Fig-S3 S5 grids decoded** (axial coords, reference labelings):
  - (d) line: 1(0,0) 2(1,0) 3(2,0) 4(3,0) 5(4,0) — 4 edges, handedness None.
  - (e): 1(0,0) 2(1,0) 3(2,0) 4(3,0) 5(1,1) — 5 edges, H {(1,4,5),(2,3,5)}.
  - (f): 1(0,0) 2(1,0) 3(2,0) 4(1,1) 5(0,1) — 7 edges, H 7-list + cond {(1,4,5)*,(3,4,5)*}.
  - (g/h): 1(-1,0) 2(0,0) 3(1,0) 4(0,1) 5(0,-1) — 6 edges, H 6-list.
- Fig-S3 renders: `gs -r400` page 11 of the SI PDF (poppler/imagemagick blocked;
  ghostscript works). S3/S4 grids (panels a-c) still to be read off the same figure.
- Working label-alignment code (hex-isometry match + reference labeling) that already
  reproduces paper counts d=0/e=2/g=6/f=7 for the geometrically-matching majority.

## Phases

### Phase 0 — Encode the Fig-S3 grids geometrically (the new whitelist)
- Read S3 (a,b) and S4 (a,b,c) grids off Fig S3 the same way as S5; encode all as
  canonical hex point-sets (axial) + reference labelings + per-grid handedness lists
  + **solid-line (required-interaction) sets** and **broken-line (≥1-mandatory)
  disjunctions** (Rule 3 data — read from the figure line styles).
- Define `geometric_canonical(points)`: the lex-min point set over the 12 hex
  isometries (6 rotations × reflection) after translating to a canonical anchor.
- Whitelist = {geometric_canonical(grid) : grid ∈ Fig-S3}. (S5: 4 entries; the
  count of *grids* is unchanged — 2 S3, 3 S4, 4 S5 — only the *test* changes.)

### Phase 1 — Reimplement SCC-2 as geometric congruence
- New `passes_scc_2(skel)`: `geometric_canonical(skel.points) in GEOM_WHITELIST[n]`.
- Keep PCC and SCC-1 unchanged (SCC-1 already forces e.g. straight paths).
- Re-verify RCC still selects one representative per equivalence class; geometric
  SCC-2 may change which candidates survive to RCC.
- Keep the old graph-whitelist behind a flag for A/B comparison and regression.

### Phase 2 — Re-enumerate and validate the new skeleton set
- Run enumeration with geometric SCC-2 → new S5 skeleton count `N'` (expected ≠ 198;
  likely fewer — bent variants dropped — but confirm empirically, both directions).
- **Diff** new set vs old 198: list the skeletons dropped (and any added). Confirm the
  dropped ones are exactly the non-Fig-S3 geometries (bent-e etc.).
- **Acceptance test:** every surviving skeleton's induced grid is geometrically one of
  the Fig-S3 grids, AND the Phase-3 handedness application yields **exactly** the paper
  count per grid ({0,2,6,7} + 9 for the type-conditional f) for **all** skeletons —
  the clean-mapping property that failed at 20% before.

### Phase 3 — Apply minimal handedness + minimal interactions (Rules 2 & 3)
- With every skeleton mapped to a Fig-S3 grid + reference labeling, emit:
  - handedness = per-grid list (geometric L/R from our signature), incl. the
    **type/sheet-conditional** f-triples (SI line 25: 1,2,3 helices & 4,5 strands→sheet).
  - interactions = **required** = solid-line set (letters); everything else optional
    `X`; non-adjacent same-sheet strands `-` (already correct). Handle the `≥1
    mandatory` broken-line disjunction (grid e) — needs a query-encoding representation
    (may require a small ProSMoS-matrix convention; verify searchmatrix supports it).
- Fold both into `assignment.py` behind `paper_faithful=True` (alongside the already-
  landed `paper_typing`). Unit-test per-grid counts.

### Phase 4 — Regenerate queries + impact assessment
- Generate the fully paper-faithful S5 query set: `N'` skeletons × valid typings
  (typing rule) × minimal handedness × minimal interactions.
- **This changes the analysis basis.** Before committing to a full re-run, do a
  *cheap impact probe*: re-run the sieve / occupancy / negative-space on a sample and
  compare headline numbers (dark fraction, D1/D2 split, negspace count) old-198 vs
  new-N'. Decide whether the darkness conclusions survive the basis change.
- If they move materially, schedule the full re-run of the darkness/negspace/rarefaction
  pipeline on the new basis (dataset_analysis_plan.md, paper_section.md, deck figures).

### Phase 5 — Reconcile with oracle 648 / document the lineage
- Confirm the CG-2012 binary's 648 are abstract-graph over-generations (K₁,₄/C₅ not
  hex-realizable) — i.e. the *binary* diverged from the *paper*; our new geometric set
  is the paper-faithful one. Update `project_s5_oracle_gap` memory + cg2012_internals.md.

## Validation strategy (no clean external reference)
- The Fig-S3 grids ARE the geometric reference. Faithfulness = (a) every surviving
  skeleton's grid ∈ Fig-S3 by congruence; (b) handedness counts hit the paper's fixed
  per-grid numbers for 100% (not 80%); (c) spot-check a handful of skeletons' emitted
  queries against the oracle's SSPs *where a genuine geometric match exists*.
- Regression: keep the graph-based path; assert the geometric set ⊆ graph set and that
  the symmetric-difference is exactly the enumerated bent/degenerate cases.

## Risks & the big decision
1. **Invalidates the current basis.** All darkness/negspace/D2/sieve results are on the
   198. Phase 4's probe gates whether a full re-run is needed — this is the load-bearing
   cost and should be an explicit go/no-go with the user.
2. **Figure-reading fidelity.** S3/S4 grids and the solid/broken interaction lines must
   be read off Fig S3 accurately; a mis-read grid changes the whitelist. Mitigate by
   cross-checking grid edge counts + the reproduced handedness counts.
3. **Broken-line disjunction** (`≥1 mandatory`) may need a ProSMoS query convention we
   haven't used; confirm searchmatrix semantics before relying on it (fallback: treat as
   all-required, a slight over-spec, documented).
4. **RCC coupling.** Geometric SCC-2 changes the candidate pool feeding RCC; re-verify
   the representative selection is stable and deterministic.
5. **Scope creep to S3/S4/S6.** Do S5 first (the analysis target); generalize after.

## Suggested sequencing
- **Step A (cheap, decisive):** Phase 0 (S5 only) + Phase 1 + Phase 2 → get `N'` and the
  100%-clean-mapping proof. This alone answers "what is the paper-faithful S5 count" and
  whether the fix is sound — a few hours, no analysis re-run.
- **Step B:** Phase 3 (handedness + interactions) + Phase 4 impact probe. Decide on full
  re-run.
- **Step C (conditional):** full analysis re-run + Phase 5 documentation.
