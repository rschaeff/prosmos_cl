# Session handoff — 2026-07-16

Corrected full-scale S5 sweeps; why β-sandwiches are dark; the Ig rescue; the
PDB arm; and three claims of mine that the data killed.

---

## 1. AFDB arm — corrected, complete

The published darkness was a **reader bug**, not biology. Same queries, same DB,
hardened reader + loop inversion:

| | bugged | corrected |
|---|---|---|
| distinct domains lit | 15,327 (0.31%) | **1,191,270 (24.2%)** |
| | | **78×** |

- 800/800 chunks, 4,921,931 domains, 2,057 core-hours, 1.53 s/record.
- 8,710,541 **match incidences** ≠ distinct domains: a lit domain matches **7.3**
  queries on average. The dark-fraction denominator is the distinct union.
- Per-cell enrichment (`s5_grid/s5_grid_3panel_enrichment.png`): **zero cells lost
  signal**, median cell gained 4×, 30% gained ≥16×, and **1,733 cells (27%) went
  0 → lit**. A strictly additive fix — strong evidence the old result was
  truncation, not a different-but-valid measurement.

### Eligibility funnel (A/B/C/D sieve, full DB)
| bucket | | share |
|---|---|---|
| A | <5 SSEs — cannot host an S5 motif | 43.0% |
| B | ≥5 SSEs but <5 pass length (H≥8, E≥5) | 13.1% |
| C | sparse — no connected 5-SSE motif | 1.1% |
| D | **compact, has a connected 5-SSE motif** | 42.1% |

**D-dark = 892,305** — passes every structural prerequisite and still matches
nothing. That is the real dark matter. Of eligible domains, 43% hit.

---

## 2. What the dark domains are — full-coverage fold enrichment

`s5_grid/fold_full.txt`. 1.1M annotated D-bucket cores (471K dark / 630K hit),
annotated via lotta `afdb_200m.ecod_cluster_summary`, named via dione
`ecod_rep.cluster`.

**Dark-enriched:** Ig β-sandwich **7.8× (110,008 dark = 23.4% of the whole dark
set)**, jelly-roll 6.9×, glycosyl-hydrolase-like 6.6×, β-propeller 4.1×,
α-β plaits 3.5×, β-grasp 2.9×, SH3 2.8×.

**Hit-enriched:** repetitive α-hairpins **0.04×**, flavodoxin-like 0.29×,
HTH 0.51×, **Rossmann-like 0.51×**, **P-loop 0.53×**.

> **The split is not α vs β — it's flat vs bilayer.** Rossmann, flavodoxin and
> P-loop are α/β and they *hit*: one central sheet with helices either side lays
> flat on the lattice. What goes dark is the **two-sheet** architecture.

Earlier sampled numbers (Ig 25×, glyoxalase 36×) were the `n_pass==5`
subpopulation — right for that population, wrong to headline. Use the full
D-bucket numbers above.

---

## 3. Mechanism — resolved

- The enumerator emits only **4 base 5-SSE topologies** (Fig-S3 grids d/e/f/g-h;
  4/5/6/7 required edges). The densest is the **gem** `[2,2,3,3,4]`, 7 edges.
- A β-sandwich 5-core is **9 edges**, `[3,3,4,4,4]` = **K5-minus-one** = the gem
  **plus the two skip contacts** (red–blue, purple–green in the hand drawing).
  It is **70%** of dark Ig n_pass=5 cores; 91.5% of those cores are all-β.
- **Root cause = the cross-sheet crossing angle.** Using the metamatrix `sheet`
  lines (`sheet <id> <n> <1-based SSE idx…>`): 89% of dark Ig cores span exactly
  2 sheets; within-sheet contacts median **21°** (tight, near-parallel — the grid
  can type these), cross-sheet median **35°** (broad, to 60°). The angled crossing
  fans each strand across several strands of the other sheet → the extra contacts
  (breaks adjacency, >7) **and** an angle the parallel/antiparallel type codes
  can't name (breaks contact type). One geometric fact, both failures.

### Size dependence — it is NOT β-blindness
Same Ig population, eligible-hit rate: **S3 93% → S4 26% → S5 14%**
(baselines ~58 / 27 / 43%). It **crosses over**: β-sandwiches are the *best*-
covered fold at S3 and the worst at S5, because a 3-strand sub-motif fits in one
sheet but a 5-strand sub-core cannot be picked without spanning both.

Ig overall at S5: **144,564 searched, 13.3% hit, 86.7% dark** — and not a size
artifact (93.9% are eligible; 14% eligible-hit vs a 43% baseline).

---

## 4. The Ig rescue — bounded

- The gate is the **sheet declaration**, proven: the *identical* 9-edge matrix
  hits **0/1503** dark cores under `sheetS 1 2 3 4 5` (one sheet) and **18%**
  under a two-sheet declaration.
- **11 merged two-sheet templates cover 80%** of dark Ig 2-sheet cores
  (measured through searchmatrix; 31 *exact* templates → 11 permissive ones).
  Curve: 1→17%, 4→52%, 11→80%, 16→86%.
- Scope: ~10–15 templates per β architecture → low hundreds for β space.

### The rule this breaks — and it is faithful, not a bug
**Corrected claim.** ProSMoS *does* allow two sheets at S5: 13 of the 198 all-β
skeletons declare two (and cross-sheet contacts exist in both the CG-2012 oracle,
15% of blocks, and our enumeration, 6%). My earlier "contact ⟹ same sheet / no
cross-sheet contacts" was **wrong**.

The real rule is the **sparse-required-contact cap**, present in *both*:
oracle ≤5 required edges, ours ≤7. A 9-contact β-sandwich core is outside both.
The empty quadrant is **dense × two-sheet** — the hex lattice hosts sparse
two-sheet fine, but cannot realise dense cross-sheet packing. So it's a
fundamental limit of the ProSMoS SSP paradigm, and the rescue is a deliberate
**vocabulary extension** to bilayer motifs.

---

## 5. PDB arm

- **57.1% of ECOD experimental domains have an S5 hit** (partial, 68% of DB) vs
  AFDB's 24.2% — because PDB is 86.5% eligible vs AFDB's 56.3%.
- **Runtime is O(N⁵), not a regression.** PDB mean **14.0** SSEs vs AFDB **6.2**;
  ≥20 SSEs: 24.2% vs 1.8%. Mean C(N,5) work/record 78,422 vs 7,934 → **predicted
  9.9×**, **observed 8.5×** (13.73 vs 1.62 s/record). Chunks are not comparable —
  PDB chunks hold 5× fewer records *because* each costs more.
- **Where the work goes** (`s5_grid/work_distribution.png`): AFDB's top **1%** of
  domains do **95%** of the work (a few giants; the rest is free); PDB's top 1%
  do **57%**, needing 10% to reach 93% — the cost is spread across *most*
  domains, not outliers you could special-case.

---

## 6. PALSSE is not the leak (meeting concern #1)

- **DSSP on 1,377 identical AFDB structures** (`dsspcmbi` from
  `~/src/Dali_v5/DaliLite.v5/bin/`): **97%** of PALSSE's <5-SSE calls confirmed
  <5 by DSSP. The bucket-A exclusion is real.
- **The bias runs opposite to the concern**: mean DSSP−PALSSE = **−1.27**
  (DSSP finds *fewer*). PALSSE over-segments. Our eligible set is if anything
  slightly generous — the conservative direction.
- **Worked example A0A0K0DT80 nD3** (res 176–375): of 17 Mol*/DSSP-apparent SSEs,
  **16 (94%)** are PALSSE-detected and pass the length filter; **1** genuine miss.
  PALSSE emits **22** where DSSP draws 17. On the *full chain* it's 78%, but the
  gap is **DPAM's domain parse** (res 1–20, 46–66, 161–175, 376–413 are in no
  domain and were never handed to PALSSE) — not SSE detection.
  Figures: `s5_grid/nd3_v1.png`, `nd3_v2.png`.

---

## 7. AFDB vs PDB — NOT YET INTERPRETABLE (was: "the comparison that survives")

> ⚠ **REVISED 2026-07-16 (later).** The fold-level numbers here were built on a
> corrupted AFDB array and a figure that drew empty cells as agreement. Full
> post-mortem: `enum/docs/s5_grid_defects_2026-07-16.md`. Summary:
> - **670 AFDB cells were silently zeroed** (53,397 folds = 23% of the AFDB
>   total) by a sharding race in the grid assembly. PDB was unaffected, so the
>   bias was **AFDB-only and downward**.
> - **42% of cells are dead in both** datasets. They evaluate to
>   `log2(NA/NP) = +0.0004` and rendered as neutral grey — *pixel-identical to
>   perfect agreement*. Being the largest block, they **pinned the all-cell
>   median at exactly 0.00**.
> - **640 queries are physically unsatisfiable** (a sheet strand needing ≥3
>   lateral neighbours) and must be masked, not counted as darkness.

**Structure-level was 100% artifact.** Normalising to rate-per-domain fixes size
but *not* redundancy. That part stands.

| | structure-level | fold-level (BUGGY — do not use) | **fold-level (CORRECTED)** |
|---|---|---|---|
| median PDB/AFDB rate ratio, all cells | +3.86 log₂ = **14.5×** | +0.00 = 1.00× ❌ | −0.585 log₂ |
| median, data-bearing satisfiable cells | — | — | **−1.22 log₂ = 0.43×** |
| AFDB-only cells | 968 | 620 ❌ | **720** |
| PDB-only cells | 65 | 698 ❌ | **135** |

Both sets are ~30× redundant for *different* reasons — PDB 34.2× (crystallography
effort: antibodies, lysozyme, kinases), AFDB 30.8× (sequence-cluster reps; mean
84 sequences each, **max 1,607,304**). Denominators are near-identical (**3,656 vs
3,655 folds**), which is what makes the fold view fair. That rationale is intact;
only the numbers were wrong.

> **RETRACTED:** *"Experimental and predicted structure space contain S5 motifs at
> indistinguishable fold-level rates."* That 1.00× was the pseudo-count of 2,677
> empty cells, not biology.

> **DO NOT simply substitute 0.43×.** The PDB sweep is **68% complete**, and
> partial coverage depresses PDB counts in exactly that direction — it is
> confounded with the effect. The honest statement today: **the fold-level
> comparison is not interpretable until `624866` finishes.** This is an
> independent reason to hold the prosmos_inspect re-export (§8).

Figure: `s5_grid/afdb_vs_pdb_fold_level_final.png` (script `plot_nr_final.py`,
array `grid_afdb_nT_rebuilt.npy`, mask `impossible_mask.npy`). Do **not** show
`afdb_vs_pdb_fold_level.png` (buggy) or the structure-level `afdb_vs_pdb_grid.png`
(redundancy artifact).

---

## 8. prosmos_inspect

**Shipped** (merged to `main`, `68908e4`): click any exemplar → its **full base
PALSSE assignment** (every SSE, length, pass/fail the H≥8/E≥5 filter, matched
slots tinted in the schematic's colours). `export/export_palsse.py` →
`data/{ds}/palsse/{did}.json`; new API route; `PalssePanel.tsx`.
Keyed by `did` (PALSSE is a domain property, not cell-dependent like contacts;
`afdb_unassigned` has `ecodUid = 0`). Self-validating: asserts the matched
segments really are SSEs of the resolved record.
Coverage: manual_pdb 5,865 (3 unresolved) · pdb_exp 11,864 (0) ·
afdb_assigned 4,320 (43) · afdb_unassigned 382 (0).

**⚠ THE APP IS BUILT ON BUGGED DATA.** `export_cells.py` reads
`s5_full_afdb/hits` + `s5_full_pdb_exp/hits` — the pre-fix sweeps. Every cell
count, exemplar list and `_nr` view is pre-correction.

Re-export scope, **established by test, not assumption**:
| dataset | verdict | evidence |
|---|---|---|
| `manual_pdb` | **CLEAN — no re-export** | hardened binary on `ecod_db_manual_v4`: **3,083 = 3,083** (exact match to the app) |
| `pdb_exp` | **BUGGED** | corrected sweep gives **52,634** for cell 90/10 from *68%* of the DB vs the app's **21,521** |
| `afdb_assigned` | **BUGGED** | its hitroot *is* the run that gave 15,327 / 0.31% |

Re-export chain (starts upstream — `full_spec.json` is itself built from hits):
`find hits → {ds}_hits.txt → build_promiscuity.py (+ uid_class.tsv) →
build_full_spec.py → export_cells.py → export_contacts / export_palsse /
export_compare / build_nr_view`.

**HOLD BOTH — do them together once the PDB sweep (`624866`) finishes.**
(Decision 2026-07-16.) Re-exporting AFDB alone leaves a window where AFDB is
corrected and PDB is stale, silently breaking every cross-dataset comparison and
`compare/` matrix in the app; and re-exporting `pdb_exp` at 68% would bake in
partial counts that *look* authoritative. §7's defects are a second, independent
reason: the fold-level comparison is not interpretable until PDB is complete.

---

## 9. Claims of mine the data killed

Recorded because each nearly reached a slide:
1. **"Handedness is the filter."** Stripping handedness rescues **2%** of dark
   cores (16/738) and **0** of the β / path-P5 cores. I called it from hit-*file*
   counts mid-run; those 198 files were 16 records × ~12 queries.
2. **"Contact ⟹ same sheet."** False — cross-sheet contacts are in both the
   oracle (15%) and our enumeration (6%). The real rule is the sparse cap.
3. ~~**"The 968-vs-65 AFDB-only asymmetry will survive."** It didn't — 620 vs 698
   at fold level.~~ **UN-RETRACTED 2026-07-16.** The asymmetry *does* survive:
   **720 AFDB-only vs 135 PDB-only**. The "symmetric 620/698" was the sharding
   bug — **563 of those 698 "PDB-only" cells were AFDB cells that had been
   zeroed**. I retracted a correct claim on the strength of corrupted data.
   (Still provisional: at 68% PDB, AFDB-only is overstated and PDB-only
   understated, so the true asymmetry is smaller than 720/135.)
4. **"1.00× — indistinguishable fold-level rates."** The median of a population
   that was 42% empty cells, each evaluating to +0.0004. See §7.
5. **"The enumeration forces `-` on non-adjacent same-sheet strands, making those
   queries self-contradictory."** No — 84% of same-sheet non-adjacent pairs in
   *real* records are coded `-` (4805/5700); the queries match reality.
6. **"The generator only emits `t`, so parallel sheets can never match."** No —
   sheet-internal codes are `t` 5,168 / `c` 3,444 / `-` 2,919.

The pattern worth internalising: **4, 5 and 6 were all mechanisms I found
convincing before testing.** The two that survived contact with data (§7's
defects) were both found by *looking at a picture* and asking why it had
structure.

---

## 10. Running / next

| job | state |
|---|---|
| PDB sweep `624866` | 280/400 (70%), 57.1% lit, 0 fails, ETA ~15h → then merge/archive/retry |
| 8 Å cutoff `626145` | quantifies how much of the 57% adjacency-mismatch is 11 Å vs 8 Å |
| ~~expanded Ig `626154`→`626411`~~ | **VOID — resubmitted as `627042`** (see below) |
| expanded Ig `627042` | 17 two-sheet templates × 800 chunks × full AFDB: recovery + off-target |

**The first expanded-Ig run (`626154`) was void, not a result.** Its classifier
reported *"0 records hit, 0% Ig recovery"* — which would have contradicted the
offline 80% coverage. Cause: `manifest.txt` was built with `ls`, **aliased to long
format**, so every line read `-rw-r--r-- 1 rschaeff lab 195 … tmpl_0.query`.
`basename` died 11× (once per template), **zero queries loaded**, and 800 chunks
searched an empty query set in 0s each — rc=0 throughout. searchmatrix exits 0 on
an empty query set, so nothing downstream could tell it from a real negative.
Fixed: manifest built with `find`; `array_inverted.sbatch` now validates every
manifest line and exits 3 (commit `deb7be0`); smoke-tested one chunk (167 hits)
before fanning out. Re-run uses **all 17** templates — `merged_tmpl.py` emits them
in greedy coverage order, so `tmpl_0..10` is the 80% set and `11..16` extends the
curve.

Then: (a) **re-export `afdb_assigned` + `pdb_exp` together** once PDB completes
(§8); (b) PDB negspace verdict; (c) rebuild the fold-level comparison at 100% PDB
using `plot_nr_final.py` — it is **not interpretable before then** (§7).

**Contact cutoff.** `generateMatrix/src/external.h:32` —
`DISTANCE_DEFAULT 11.0 // 8 A this parament need to be adjusted` (contact =
`overlap > 2.5 && min CA-CA <= 11`, computed over the **full CA trace**, so the
DB's 2-endpoint coords can't reproduce it — the 8 Å test needs the real pipeline).
Rebuilt binaries: **must be `-O0`** (`-O2` segfaults; the 11 Å rebuild is
byte-identical to production).

**DB access** (do not guess usernames — the classifier blocks credential probing):
user `ecod`, password in `~/.pgpass`. Annotations lotta:45000 `ecod_protein`
schema `afdb_200m` (`ecod_cluster_summary`: cluster_name → uniprot_acc,
domain_idx, t_group_id, cluster_size). Fold names dione:45000 `ecod_protein`
schema `ecod_rep` table `cluster` (`type='X'`). PDB uid → T-group:
`pdb_exp_build/uid_class.tsv` (`pdb_analysis.domain` is empty on lotta, absent on
dione).

Artifacts in `work/prosmos_2026/s5_grid/`. Deck (10 slides):
https://claude.ai/code/artifact/b9626532-ac7f-45f7-88b9-a5ee045704bb
