# Session handoff 2026-08-01 — the geometric census (H-bond-overlap), motif 17, /embedding

Continues the Ruczinski four-strand census work (`~/work/prosmos_2026`, scoped work
repo on `master`; browse in `~/dev/prosmos_inspect` on `corrected-sweep-reexport`).
The session turned the Figure-2 census from a strand-*count* accounting into a
geometrically honest one, and made motif 17 the headline. **Running record:
`~/work/prosmos_2026/paper/FACTS.md` §3.1b/§3.2.** Memory:
[[project_geometric_census]].

## What landed

1. **The geometric filter → `census_hbond` is now definitive.** One hit = 4
   pairwise-disjoint strands in one β-sheet **AND** passing Nick's H-bond-overlap
   test: for every interior strand, its two flanking c/t neighbours must overlap
   along its axis (proxy for H-bonding an overlapping stretch → opposite faces → a
   real contiguous sheet). Removes matches that satisfy the pairwise code but aren't
   a sheet. `ruczinski/census_hbond.py → census_hbond.json`.

   | measure | realized/96 | dark |
   |---|---|---|
   | complete (=4) | 74 | 22 |
   | **core (≤6)** | **89** | **7** |
   | raw (any) | 96 | 0 |

2. **The vector-product "both-sides" test is REJECTED as twist-sensitive** — the
   key methodological lesson. Full-perpendicular version *over*-rejects (killed
   motif 43 to a false raw-0, inflated core-dark to 9; `census_geom`/`census_bothsides`,
   both dropped). Fitted-sheet-normal version *under*-rejects (passes 8vqe, motif 17's
   own failing example). H-bond-overlap is coordinate-based, needs no normal, and
   reproduces Nick's 8vqe acid test. Do not resurrect the vector product.

3. **The 7 core-dark residual, led by motif 17.** Motifs 17 (panel 26), 19 (28),
   25 (41), 26 (45), 29 (46), 41 (43), 72 (34). Split by *why*: **geometry-driven**
   17/25/26/41 (their ≤6 matches fail H-bond-overlap) vs **sheet-size** 19/29/72
   (only ever in ≥7 sheets). Motif 43 & 61 are near-boundary *survivors* (core 1
   each) the over-strict test had wrongly killed. All 96 realize at raw — no absentee.

4. **Core rarefaction re-derived with the geometric filter** (`rarefaction_core.py`):
   complete S_obs 73 / Chao1 82; **core S_obs 87 / Chao1 90 / Good's 0.9985**
   (near-saturated, CI lower bound ≈ S_obs); **raw S_obs 96 / Chao1 96 / Good's 1.000,
   f1=0 — provably closed.** "Census not sampling-complete" is a complete-sheet artifact.

5. **8vqe figure — motif 17's geometric failure** (`figures/fig_8vqe_motif17.py`,
   from e8vqeB3's real strand-axis endpoints). (A) β1∥β3∥β4 a genuine parallel
   3-strand sheet (|axis·axis| 0.87–1.0), β2 ~perpendicular (≤0.19). (B) the test:
   interior β4 FAILS (neighbours β2 [-1.2,1.9] & β3 [7.2,18.6] → 5.3 Å gap), interior
   β3 PASSES (β1 & β4 overlap +12.1 Å). Cited in `draft_results_3.2_core.md` §3.2.

6. **`/embedding` current end-to-end** (prosmos_inspect). Repointed at census_hbond;
   `geom_valid` threaded through exemplar selection (panel-ordered positions + 3D
   coords), smallest-host-first, index leads with motif 17. Added `coreDarkReason`
   ('geometry'|'sheet-size') to export+types; MotifHeader branches on it (names 8vqe
   for 17, notes the low-pLDDT AFDB counterpoint). Table page swapped to
   `fig2_triptych.svg`, prose rewritten (22→7→0), dead `fig2_dual.svg` dropped.
   Verified live: `/embedding`→307→`/embedding/17`; typecheck clean. **Montages
   regenerated** (all 22 `.pse` from the clean exemplars; gitignored, served from disk).

## Key numbers to trust

- Census (experimental): **74 / 89 / 96**; dark **22 / 7 / 0**; raw provably complete.
- 7 core-dark: 17,19,25,26,29,41,72 — 4 geometry (17/25/26/41), 3 sheet-size (19/29/72).
- Rarefaction core Chao1 **90** (Good's 0.9985); raw **96** (Good's 1.000, f1=0).
- 8vqe: β4-interior overlap **−5.3 Å (FAIL)**, β3-interior **+12.1 Å (PASS)**.

## Earlier in this session (context for the paper)

- **Deposition/method analysis** (from the pre-consolidation part): forbidden-realized
  topologies are ~47% cryo-EM (3.4× baseline), median 30 domains/PDB (EM 58 vs X-ray
  12) — realized in large multidomain/complex contexts. Two clusters: cryo-EM-new
  (e.g. motif 3, one family) vs never-really-absent (motif 56, diverse, pre-2002). The
  Dunbrack/PISCES 2002 *literal* set is not available; deposition ≤2001 used as proxy.
  Candidate story to tell alongside prediction: some topologies only realized in
  certain multidomain/complex contexts.

## Open / not done

1. **Nick's featured domains into /embedding** — H0UML0_D7, R5FFZ3_D3, A0A1H8P2A6_D4,
   A0A2N6P4J1_D5 (his picks). Not yet added as featured pins / fold-composition axis.
2. **Merge `draft_results_3.2_core.md` into MANUSCRIPT.md** — the draft ends with a
   downstream-edit checklist (abstract "not a saturated space", conclusion "lone
   absentee motif 19", §3.2b/c recompute on the geometric core set, abundance medians).
3. **§3.2b/c** (X-group monotypy, predicted-only) still use the complete-sheet census;
   recompute on the geometric core hit set.
4. Optional: real PyMOL cartoon panel for the 8vqe figure (structure-file search on
   `/data/ecod/...` was slow; the coordinate schematic is self-contained for now).
5. Decide whether the cryo-EM / multidomain-context finding goes in the paper.

## Commits this session

prosmos_2026 (master): `5e9f3f7` (triptych consolidation), `8fb7a28` (core rarefaction),
`1145469` (geometric H-bond-overlap census), `77092fd` (8vqe figure). Dropped:
census_geom/census_bothsides (twist-sensitive dead ends).

prosmos_inspect (corrected-sweep-reexport): `b71ddc7` (consolidate), `a5f673e`
(geometric filter, exemplars real, lead with 17), `0436ae1` (prose+figure+coreDarkReason).
Montages regenerated on disk (not committed — gitignored).
