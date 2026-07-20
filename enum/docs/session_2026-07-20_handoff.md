# Session handoff — 2026-07-18 → 07-20

Continues `session_2026-07-16_handoff.md`. That one covers the corrected AFDB
sweep, the eligibility funnel, the Ig rescue and the grid defects. This one covers
the completed PDB sweep, the AFDB-only deliverable and its controls, the app
re-export, and the scale (n=3/4/5) question.

---

## 1. Sweeps — BOTH COMPLETE

**PDB `624866` finished: 400/400 chunks, 496,359 records, 295,360 domains lit
(59.5%), 13,824,959 motif rows, 3,180/6,336 queries hit.** Archive at
`s5_pdb_inv/archive/{hits.tsv.zst, distinct_hitters.txt, query_counts.tsv}`.

**chunk_0059 rescue — the reusable lesson.** It is the single heaviest chunk
(median 35 SSEs/record, rank 1/400) and timed out at 24h. A plain retry is
futile: deterministic, same wall, and `retry_missing.sh` defaults to 24h, so the
driver would loop it forever. Fix: split into 8 COST-balanced sub-chunks
(`s5_pdb_inv/split59/`, 15.6e9 each vs the chunk's 124.6e9), ~3.6h each, then a
fold job unioned the hits and wrote a synthetic `parts/60.tsv`. The 8 sub-hitparts
were consolidated into one canonical `chunk_0059.db.tsv.zst` so the merge's
completeness guard (want==have vs chunks.list) passed. **If a chunk times out
again: split, don't retry.**

---

## 2. The #1 deliverable — AFDB-only cells

`enum/docs/afdb_only_cells_2026-07-17.md`, deck section "the depth control".

Naive: **619 AFDB-only cells** at 100% PDB (fold level, satisfiable only;
3,624 AFDB-lit / 3,173 PDB-lit / 3,005 shared / 168 PDB-only).

**Size-matched rarefaction is the load-bearing control.** AFDB searched 4.92M
records vs PDB 496k (10x; 4x in lit structures). Thin AFDB and recount:

| AFDB depth | AFDB-only | PDB-only |
|---|---|---|
| 100% | 619 | 168 |
| 50% | 426 | 302 |
| **25%** (≈PDB's lit depth) | **284** | **501** |
| 10% | 148 | 842 |

**At equal depth the asymmetry reverses.** Sanity gate: rebuilding the grid from
the (cell,record) pairs reproduced `grid_afdb_nT_rebuilt` at 100.0% pattern
agreement, r=1.0000. Persistence over 20 matched-depth draws leaves **90 robust**
cells, 407 middle, 122 pure artifacts.

Deliverable artifacts: 60 rendered exemplars across 4 legible plates
(`s5_grid/renders/gallery_{1..4}.png`), 3 detail cards, adjacency context
(`afdb_context.txt`). **55 of 60 sit one H<->E flip from a PDB-populated named
fold** (mean 2.8 such neighbours).

---

## 3. Controls, all of which held

- **Fisher (asked for by Qian) + why it misleads.** Per-cell 2x2 on FOLD counts,
  BH. It calls 1,457/3,792 cells (38%) significant — the null is wrong. Worse,
  it is directionally crippled: PDB-only cells carry 1-5 folds so
  `[[0,3656],[4,3931]]` can never reach significance, while AFDB-only cells reach
  17 folds *because AFDB was sampled deeper*. "PDB=0 & AFDB significant" is
  findable (42 cells); "AFDB=0 & PDB significant" is **structurally impossible**
  (0). Replacement offered in the same currency: thin AFDB to PDB depth; a fold
  with k members survives with probability 1-(1-q)^k, so the thinned count is
  **Poisson-binomial** — exact, no resampling, BH-correctable. Drops 42 -> 5.
  **Four cells significant under both:** sk16/ty12 (17 folds), sk171/ty7 (11),
  sk158/ty27 (10), sk34/ty0 (8). `fisher_cells.py`, `fisher_zero_cells.json`.
- **pLDDT.** Read from the B-factor of the CHOPPED DPAM domain (not
  `ted_domain.plddt` — TED draws different boundaries). AFDB-only 83.6 vs shared
  84.1, Mann-Whitney p=0.41. **Orthogonality result:** cells the depth control
  REJECTS average 84.6, *higher* than those it keeps (84.0). Confidence and
  sampling depth are independent failure modes.
- **Context (EM / viral), fold level.** ZERO of 1,725 cells enriched for either.
  Domain level said 334 EM-enriched + 328 viral-enriched + 60 viral-only (top
  n=168, p~1e-150) — all redundancy: every viral-only cell is ONE T-group
  (hepatitis C protease, poliovirus capsid, coronavirus spike) solved hundreds of
  times. ~100% artifact rate. `enum/docs/context_stratification_2026-07-20.md`.
- **PDB vintage.** Rebuild the PDB side by release year: AFDB-only goes
  2,668 (1995) -> 619 (2026); **77% of the 1995 set filled in**. But the rate
  collapsed 22x (159/yr -> 7.2/yr). Fit `553 + 2264*exp(-0.116(t-1995))`,
  half-life 6.0y, floor ~553. **Falsifiable: ~593 AFDB-only cells in 2030.**
  Consistent with depth-limitation, not a rival to it. `pdb_vintage.py`.
- **Taxonomy.** Eukaryote enrichment declines monotonically with strand count
  (40% of all-helix cells enriched -> 0% of 5-strand) and SURVIVES controlling for
  each cell's own fold composition. The larger residual is *prokaryote* enrichment
  in strand-rich cells (167 cells vs 113). `tax_cells.py`.

---

## 4. The instrument's ceiling — the most transferable result

**β-flower test.** Durairaj/…/Pereira, *Nature* 2023 (doi 10.1038/s41586-023-06622-3,
PMC10584680) found a genuinely novel fold, the **β-flower**: `A0A494VZL1`
(prototype, a singleton), `A0A0S7BXY3` (circular permutant), Pfam PF21784-6,
clan CL0395. We ran it through our exact pipeline (PALSSE -> generateMatrix ->
all 6,336 queries). It lights 5 cells, ALL shared; the prototype lights exactly
one cell, `sk132/ty31`, which already holds **366 PDB folds**.

**TED at scale.** 28,902 CATH-unassigned domains: **99.66% land in PDB-occupied
cells** (0.34% vs 0.18% touch an AFDB-only cell; OR 1.83, p=4.5e-4 — real and
substantively nil). And **52.5% of our lit AFDB domains are TED-novel**, so
"CATH-unassigned" describes half of predicted structure space and is a coverage
statement about a PDB-derived classification, not a discovery class.

**Conclusion:** novelty at the fold level and novelty at the 5-SSE topology level
are different quantities. A new fold is not a new 5-SSE topology — five strands
from a 13-strand barrel look like five strands from anything else. S5 motif search
is not a fold-discovery instrument. This is the answer to "why didn't ProSMoS find
new folds", and it is a category statement, not a failure of execution.

**Singleton caveat (RS's, and it is correct).** Our AFDB set is non-singleton
cluster reps: **10.6M of 15.5M clusters are singletons and were excluded — more
than twice what we searched**, and that is where novelty concentrates. So all
"no novel topology" claims are conditioned on sequence space that already has
relatives. BUT the β-flower (itself a singleton) shows the filter is not why we
missed it: searching every singleton would still have put it in a saturated cell.
Two separate ceilings — sampling and vocabulary — and only the second is fatal.
**Singleton structures are NOT available locally** (grey's store is explicitly
`non_singleton_4p9m_structures`; `set12` is ECOD symlinks), so a singleton arm is
fetch-bound: ~50k EBI downloads minimum. Recommended as a bounded, pre-registered
diagnostic if at all, not a discovery run.

---

## 5. prosmos_inspect re-export — DONE (branch, not merged)

Branch `corrected-sweep-reexport`, commit `34efb7c`. **NOT merged to main.**

- **pdb_exp**: 1.35M -> 5.31M hits, 3,174 cells, 77,071 contacts, 25,211 PALSSE.
  The corrected sweep is local-hits-mode so there is NO hit tree: exemplars carry
  pre-resolved segments in `full_spec` (from hitparts TSVs, 100%) and `hitroots`
  is deliberately empty.
- **afdb_assigned needed a REDESIGN.** The old dataset keyed on numeric
  `ecod_af2_pdb` uids; only **1.5%** of the corrected sweep's accessions exist
  there. Identity moved to **UniProt accession + DPAM domain**, following
  `afdb_unassigned`. 3,952 cells, 71,029 exemplars, 71,029 contacts, 32,788 PALSSE.
- **Structures: 32,788 / 32,788 (100%)** extracted from
  `/home/grey/afdb.200m/non_singleton_4p9m_structures` into
  `work/prosmos_2026/afdb_domain_struct/` (3.5 GB). Already chopped to the DPAM
  range (verified against `ecod_domain_range` AND against an independently fetched
  EBI v6 model — coordinates identical, so the v4-searched/v6-shown disclosure
  stays accurate). Index: `/home/grey/resources/afdb200M/356M_pdb_folder_maps/*.list`
  (dir, tar, accession) with **tar name = accession[-4:-2]**. Copied out rather
  than symlinked so the app does not depend on another user's tree.
- Code: new `AFDB_DOMAIN_DIR`; structure route serves afdb_assigned by `did`;
  contacts route guard widened to `\w{1,32}`; `export_contacts` dataset-aware;
  `export_palsse` pointed at `metamatricesDB.clean` with a direct-name mode
  (it had been silently resolving **0 of 32,788** — old regex expected `_nD<n>`,
  ours are `_D<n>`); viewer keys AFDB by did/accession.
- `datasets.json` synced: pdb_exp nOccupied 2466->3174, afdb_assigned 2074->3952,
  pdb_exp_nr N None->3935, afdb_assigned_nr N 1086->3656 (1086 was the *bugged*
  fold count — the _nr view had been normalising against a denominator 3.4x small).
- All 15 compare matrices regenerated.
- Backup of the pre-re-export app data: `work/prosmos_2026/inspect_data_backup_20260718/`
  and `data/afdb_assigned_bugged/`.

---

## 6. RUNNING / IMMEDIATE NEXT

**Scale question — "why n=5?"** The only justification was operational
(unsaturated but computable) and the comparability-to-Chitturi argument is void
(our geometric enumeration gives 140 vs their 198). Running n=3 and n=4 on a 10%
chunk subsample of both DBs to get a three-point curve.

- queries: `work/prosmos_2026/queries_n34/{s3,s4}` (40 and 672) from
  `enum/scripts/generate_enumerated_queries.py --dims 3 4`. **Verified the same
  generator reproduces all 6,336 swept S5 queries byte-identically**, so the
  three points are measured under identical constraints (incl. handedness — 5,952
  of 6,336 S5 queries carry chirality directives; an inconsistency there would
  have made the whole trend an artifact).
- jobs: `628683` n3-afdb, `628684` n3-pdb (**both COMPLETE, 80/80 and 40/40,
  rc=0**), `628685` n4-afdb, `628686` n4-pdb (running).
- record maps for like-for-like n=5 restriction:
  `s5_grid/{afdb,pdb}_sub_records.pkl` — 491,963 AFDB and 49,640 PDB records.

**Three quantities to compute across n=3,4,5 on the identical record sets:**
1. **Sharing** — mean distinct T-groups per lit cell. The criterion argued for:
   the largest n at which motifs are still REUSED across unrelated folds rather
   than being per-domain fingerprints. At n=5 sharing is heavy (median 5 T-groups
   per cell; the big cells span 973 X-groups), which suggests 5 may be *below* the
   diagnostic optimum.
2. **Saturation** — fraction of n-space lit (n=5 is 64% of satisfiable and still
   climbing with depth: 2,479 -> 3,624 over 10% -> 100%, so NOT saturated).
3. **AFDB-vs-PDB agreement** vs n. Flat => the alphabet claim is scale-robust.
   Narrowing => divergence begins, and that locates the alphabet/grammar boundary
   and would justify the ~15x cost of n=6.

RS expects S3 to be a boundary condition (or to expose an enumerator problem) —
useful as a curve anchor either way, and it was cheap.

---

## 7. Other open threads

- **TED novel-fold search**: `work/prosmos_2026/ted_novel/metamatricesDB.clean`
  (6,681 records from 7,427 TED novel models, built 2026-07-05) was being searched
  against all 6,336 S5 queries; result not yet collected. This is a *direct* test
  (no range-overlap join) of the TED conclusion in §4.
- **8 Å cutoff `626145`** TIMED OUT at 12h with no result. Needs a longer limit or
  a split before it can say how much of the 57% adjacency mismatch is cutoff-driven.
- `export_compare` uses a two-proportion z-test; Fisher would be better on the
  sparse cells (see §3).
- The app branch is unmerged and `data/` is gitignored, so nothing is live.
- Deck: https://claude.ai/code/artifact/0aaa2106-745b-487e-b4fc-5153054a03d4
  (rebuild with `s5_grid/build_deck.py`, then republish the same file path).

---

## 8. Traps hit this session (all silent, all cost real time)

1. **`ls` is aliased to long format** — building a manifest with it produced
   `-rw-r--r-- ... file.query` lines, `basename` died 11x, **0 queries loaded**,
   and 800 chunks searched an empty query set reporting "0 hits" with rc=0. Now
   guarded: `array_inverted.sbatch` validates every manifest line (commit
   `deb7be0`). Use `find`, never `ls`.
2. **`ls dir/*.pdb | wc -l` overflows the arg list** and reports 0 — I twice
   mistook a full extraction for a failure. Use `find`.
3. **Relative paths in a chunks.list** — jobs `--chdir` elsewhere, so all 240
   tasks failed instantly. The existence guard caught it (zero compute wasted).
4. **`afdb_200m.taxonomy` contains NO viruses** (AlphaFold-derived; AFDB v4
   excluded viral proteomes). Routing PDB domains through it returns 0 viral *by
   construction* — the first viral run reported "global viral rate 0.00%".
   Correct source: **`ecod_commons.domain_taxonomy` on dione** (PDB-native, keyed
   by `ecod_uid`, 3.1M rows, 159,444 viral).
5. **Current UniProt resolves only 34% of our AFDB accessions**, and survivors are
   51% eukaryotic against a true 27% — it would have *manufactured* the expected
   eukaryote signal. Correct source: **`afdb_200m.protein`** (214.7M rows,
   `uniprot_acc` -> `organism_tax_id`, 99% coverage, contemporaneous).
   Found only after I truncated my own column search at 20 tables and wrongly
   concluded the path was dead — do not truncate schema searches.
6. **Several afdb_200m tables are EMPTY scaffolding**: `eda_domain_base`,
   `deck_s8_cluster_tax_profile`, `domain_plddt`. Check row counts before designing
   around a table.
7. **Domain-level counting inflates everything.** It produced 334 EM-enriched,
   328 viral-enriched and 60 viral-only cells that all vanish at fold level. I
   walked into this immediately after warning about it in the Fisher section.
   **Count folds, not domains.**
8. `scipy.optimize` hits an MKL threading error (`__kmpc_global_thread_num`) on
   this box; numpy-only fits work.
