# Dataset input analysis — resolving the rarefaction puzzle

**Status:** active fork (opened 2026-07-09). Parks the confusing rarefaction
slide; replaces it with a completeness-audited, confound-controlled comparison
of motif-hit propensity across the three search DBs.

## The problem
The rarefaction curve made it look like the AFDB sweep was incomplete
(4.9M structures → 12,511 hitting). Three hypotheses:
- **(a)** AFDB genuinely low on motif hits.
- **(b)** PDB genuinely high (redundancy / crystallization bias).
- **(c)** something incomplete (build loss, search truncation, resolution loss).

Not mutually exclusive. Plan is ordered to kill (c) first, then decompose the
real gap into (b) vs (a).

## Reference numbers (established this session)
| stage | manual | pdb_exp | afdb |
|---|---:|---:|---:|
| DB entries (`grep -c .ssd`) | ? | 496,359 | 4,921,931 |
| ≥5 SSE (capable) | ? | 445,838 (89.8%) | 2,797,352 (56.8%) |
| distinct hitting structures | 11,590 | 72,180 | 11,971 |
| hit rate (of capable) | ? | 16.2% | 0.43% |
Per-cell head-to-head (per capable structure): all-α HHHHH ~10× PDB-favored;
α/β Rossmann (HEHEH) ~99× PDB-favored → cell-type-differential = fingerprint of
a compositional/sampling difference, not a uniform artifact.

---

## Phase 0 — Completeness & provenance audit (rule out c)

- **0.1 Build-loss rate** — ✅ DONE. pdb_exp: 500,313 input domains → 496,798 raw
  → 496,359 clean (**99.2% retained**). afdb: ~4.9M input (106,341 tarballs ×~46
  models) → 4,921,936 raw → 4,921,931 clean (**~100%**; `.clean` removed only 5).
  Both builds are essentially lossless; malformed-entry concern is negligible.
  (Note: afdb entries are DPAM domains/whole-model reps — characterize the exact
  unit in Phase 1, it affects denominator semantics.)
- **0.2 Sweep coverage/runtime** — ✅ DONE. afdb & pdb_exp each cover exactly
  6,336 distinct queries, 0 rc≠0, 0 dupes. Runtimes: afdb median 724s/max 1353s;
  pdb_exp median 1382s/max 2078s — all ≪ 1h wall, no timeouts. afdb runs *faster*
  than pdb_exp despite a ~10× larger DB → consistent with fewer real matches, not
  skipped entries.
- **0.3 searchmatrix capacity audit** — ✅ DONE. `searchControl.h:277`
  `while(inputfile1.getline(intMline,...))` streams the DB line-by-line; no
  entry array, no entry-count cap. Cannot silently truncate the scan. (Caveat:
  1MB `intMline` buffer vs 10MB getline arg — only >1MB matrix lines, i.e.
  ~>1400-SSE entries, risk overflow; rare, and .clean + rc=0 argue it didn't
  bite. Worth a spot check for giant entries.)
- **0.4 Positive control** — ✅ satisfied by the A0A044UCZ4 trace: a specific
  afdb structure was followed end-to-end to a *correct* match (cell sk133/ty0,
  SSE3↔SSE5 coded `u`, min CA ≤11 Å per generateMatrix), and the permissive
  HHHHH cell (s5-0132-0000) drew 7,372 afdb hits. The afdb search path
  demonstrably finds correct matches. [Optional stronger control: confirm afdb
  Rossmann-classified domains hit the Rossmann cell s5-0090-0010 — not needed for
  the verdict.]

## Phase-0 VERDICT: (c) is ruled out.
Build essentially lossless (0.1), sweep complete with sane runtimes (0.2),
searchmatrix streams the whole DB with no entry cap (0.3), and the search
demonstrably finds correct matches (0.4). **The 4.9M → 12,511 reduction is a real
biological/sampling result, not a pipeline artifact.** The apparent gap is a mix
of (a) AFDB S5-sparsity and (b) PDB redundancy+bias — quantified in Phases 1-3.

## Phase 1 — Input characterization (all three) — ✅ DONE
Parser `scratchpad/phase1_profile.py` (streams metamatricesDB; robust to the two
format quirks: non-A chain letters `EP`/`HX`, and residue numbers gluing to the
type on wide values `EA1083`). N matches DB entry counts exactly.

| dataset | N | median SSE | median span (res) | all-α | mixed α/β | all-β |
|---|---:|---:|---:|---:|---:|---:|
| manual_pdb | 19,015 | 12 | 129 | 11.6% | 85.0% | 3.5% |
| pdb_exp | 496,359 | 12 | 139 | 10.6% | 86.4% | 3.0% |
| **afdb** | 4,921,931 | **5** | **67** | **28.3%** | **64.9%** | 6.9% |
*(composition among ≥5-SSE capable; afdb header-count mismatch=0, clean parse)*

**Two large structural differences drive the hit-rate gap:**
1. **Size (dominant).** AFDB reps median 5 SSE / 67 res vs PDB 12 SSE / ~135 res;
   half of afdb has ≤5 SSE. Combinatorics: a 5-SSE structure has ~1 candidate
   5-SSE window; a 12-SSE domain has up to C(12,5)≈800. Hit propensity scales
   super-linearly with size, so a tens-of-× per-structure gap follows from a
   2.4× size difference alone.
2. **Composition (~2.7×).** AFDB capable is 28.3% all-α vs ~11% (PDB), and less
   α/β (64.9% vs 86%). All-α → generic grab-bag cells; α/β → S5-rich Rossmann
   cells. So AFDB under-reaches exactly the cells PDB is enriched in — this IS
   the α/β-99×-vs-α-10× cell-type differential, visible directly in the input.

**Reframing:** the low AFDB hit rate is not incompleteness (Phase 0) and not
"AFDB lacks folds" — AFDB non-singleton cluster reps are *small, more-helical
modules* (median 67 aa; conserved fragments / single small domains), while PDB
ECOD domains are large curated α/β folding units. Next: Phase 3 size-binned hit
rate to quantify the size vs composition split; Phase 2 redundancy for the census
effect.

## Phase 2 — Redundancy structure (crux of a vs b)
- pdb_exp = redundant census; afdb = dereplicated cluster reps. Quantify distinct
  ECOD T-groups among hitting structures + copies-per-T-group.
- Dereplicate pdb_exp to one-rep-per-T-group (or ~40% id) and recompute hit rate.
  If it collapses toward afdb → (b) redundancy is the driver.

## Phase 3 — Like-for-like hit rate
- Control simultaneously for size (≥5 SSE) AND redundancy (dereplicated);
  recompute matched hit rates, resolved by cell composition (α / α-β / β).
- Residual after matching = the genuine (a) component (AFDB diverse but S5-sparse).

## Phase 4 — Replace the rarefaction
- One artifact: completeness funnel + matched hit-rate table with the
  redundancy/bias decomposition explicit. Leads with "the sweep is complete."

## Data pointers
`afdb_db/`, `ecod_db_pdb_exp/`, `s5_full_afdb/summary.tsv`,
`s5_full_pdb_exp/summary.tsv`; SSE count via `grep -oP '\.ssd\s+\K[0-9]+'`;
searchmatrix `searchMatrix/src/searchControl.h`.
