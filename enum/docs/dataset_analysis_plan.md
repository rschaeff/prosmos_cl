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

**Two structural differences (NOTE: Phase 3 showed size is minor, ~1.7×, not
dominant — redundancy is the real driver; see below):**
1. **Size.** AFDB reps median 5 SSE / 67 res vs PDB 12 SSE / ~135 res; half of
   afdb has ≤5 SSE. Hit propensity rises with size (pdb_exp 5-SSE 2.5% → 25+
   24.6%), but standardizing only takes 16.5%→9.82% (1.7×). Not the main effect.
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

## Phase 2 + 3 — Redundancy & like-for-like — ✅ DONE (answers the puzzle)

**Size standardization (Phase 3a).** Applying pdb_exp's size-specific hit rates to
afdb's size distribution: expected afdb rate = **9.82%** vs pdb_exp overall 16.5%
— so size (afdb smaller) is only a **1.7×** effect. afdb ACTUAL rate = **0.55%**
(15,343 distinct hitting entries / 2,797,352 capable). → **18× residual after
size**, so size is NOT the dominant driver (corrects the Phase-1 hunch).
pdb_exp size-binned rate (the mechanism): 5-SSE 2.5% → 25+ SSE 24.6%, monotone.

**Redundancy / fold level (Phase 2 — the real answer).** Dereplicating pdb_exp by
ECOD T-group:
| | pdb_exp | afdb |
|---|---:|---:|
| capable entries (T-assigned) | 382,720 | 2.8M (mostly dark/unassigned) |
| per-structure hit rate | 18.8% | 0.55% |
| **distinct hitting FOLDS (T-groups)** | **1,360** | **1,086** |
| copies per hitting fold | 52.9 | 14.1 |

**At the fold level they are comparable — afdb recovers 1,086 of the 1,360 S5-
hitting folds that ALL of experimental PDB does (80%), from dereplicated reps.**

**Decomposition of the ~34× per-structure gap** (not missing folds, not an
incomplete search):
1. **Redundancy** (dominant): PDB is a crystallization census (~53 entries/fold);
   afdb is dereplicated (~14/fold) → ~3.8×.
2. **Size**: afdb median 5 SSE vs 12 → ~1.7×.
3. **Denominator composition**: afdb's 2.8M capable = the whole dereplicated
   proteome, mostly small/dark non-fold clusters; pdb_exp's 383k = curated folds.
   (Plus afdb 2.7× more all-α → under-reaches α/β cells.)

## FORK RESOLVED — corrected (dereplicated) rarefaction
The original rarefaction (cells vs hitting *structures*) conflated diversity with
redundancy: PDB accumulates cells slowly per structure *because* it is redundant
(copies of one fold add no cells), so "AFDB faster per structure" was really "PDB
is redundant." Corrected by dereplicating the x-axis → **cells vs distinct FOLDS
(T-groups) sampled**:
- pdb_exp: 1,360 folds → 2,466 cells (**1.81 cells/fold**)
- afdb:    1,086 folds → 2,074 cells (**1.91 cells/fold**)
- matched fold counts track within ~15% (200f: 1,279 vs 1,087; 500f: 1,775 vs 1,490)

Verdict: **AFDB folds are not cell-poorer** — per fold they reach S5 cells as
efficiently as PDB folds. AFDB's lower total cell count = fewer folds sampled,
not weaker folds. This is the defensible form of the parked rarefaction's shape
claim. Retired the per-structure curve from the deck (superseded by the fold-bar
+ diff-heatmap); the corrected fold-curve is a backup if a trajectory view is
ever wanted (`plot_foldrarefaction` derivable from the promiscuity fold→cell map).

## Phase 4 — Replace the rarefaction — the defensible statement
"Per-structure hit rate is not comparable across a redundant crystallization
census (PDB) and a dereplicated, smaller-on-average diversity sample (AFDB). At
the fold level, AFDB recovers 80% (1,086/1,360) of the distinct S5-hitting folds
that all of experimental PDB does. The 34× per-structure gap is redundancy (3.8×)
× size (1.7×) × denominator composition — the sweep is complete and no folds are
missing." → Figure: fold-level bar (1,360 vs 1,086) + the decomposition, NOT the
per-structure rarefaction curve.

## The "size 1.7×" factor is a rep-selection artifact (deep dive, 2026-07-09)

The Phase-3 size effect turned out **not** to be biology or DPAM partitioning —
it is that the AFDB search substrate is grey's **non-singleton cluster
representatives**, and those reps are systematically truncated relative to the
domains they represent.

**Chain of elimination.**
- afdb searched entries: median **5 SSE / 67-res SSE-span** vs pdb_exp ECOD
  domains 12 SSE / 139 res.
- NOT DPAM over-splitting: afdb_200m's DPAM over *all* AFDB = median **132 res**,
  ≈ ECOD 135; grey's DPAM ≈ afdb_200m's DPAM for the same proteins (coverage 0.98
  when the searched `_D<n>` bulk is joined to `ecod_domain_range` by index).
- NOT the singleton/non-singleton split: within-cluster, member domains =
  population (median ~130); only reps are short. Split is minor + wrong-direction
  (non-singleton 64 > singleton 47).
- **It is the clustering rep-selection.** Within the same clusters (roster p0,
  equal footing): **REP median 56 res vs non-rep MEMBER median 127.** Reps are
  ~half the length of the members they represent.

**Conditioned on T-group — worst for globular α/β, not repeats** (rep/member
median length ratio):
| fold | ratio | | fold | ratio |
|---|---:|---|---|---:|
| TIM barrel 2002.1.1 | **0.26** | | ARM repeat 109.4.1 | 0.77 |
| P-loop 2004.1.1 | **0.39** | | Ig 11.1.1 | 0.90 |
| kinase 206.1.1 | **0.39** | | SH3 4.1.1 | 0.87 |
| Rossmann 2003.1.1 | **0.48** | | winged-helix 101.1.2 | 0.99 |
| RNaseH 2484.1.1 | **0.55** | | OB-fold 2.1.1 | 1.00 |

Strongest truncation in the compact, S5-rich workhorse folds; least in repeats
(a short ARM rep = fewer units, mild); none in naturally-small folds (Ig/OB/SH3
already minimal — this is the eukaryotic-tilt control, and it comes out clean).
Visuals: `enum/docs/figures/s5_rep_truncation.png` (rep domain on full model,
globular reps are fragments, Ig reps whole).

**Nuance (heterogeneous clusters).** The clean 0.5 ratio is an *aggregate*
(rep-median 81 vs member-median 168 across all 2003.1.1 clusters), partly a
*between-cluster* effect. Within a single cluster the rep sits at the **short
mode** of a *bimodal* size/SSE distribution — clusters lump partial fragments +
full-length domains of the same fold, and the rep tends to be a short member. E.g.
cluster A0A2D6GC17_D1: rep 82 aa / 8 SSE, but the cluster contains full 179 aa /
**17 SSE** Rossmanns it omits. (`rep_truncation/renders/cluster_violin.png`,
`rep_vs_member.png`.)

**Partial-domain flag exists but only partly explains it.** Grey has
`summary/partial_domains.parquet` (query = fragment contained in larger target;
qcov high, tcov low) + `dpam_prob_max`. Partial-flag rate among reps is enriched
in globular folds (TIM 8.3%, P-loop 8.0%, RNaseH 7.0%, kinase 6.5% vs 4.9%
baseline; Ig 1.1%) — corroborates the direction — but absolute rates are modest,
so most short reps are *unflagged* short members, not DPAM-called partials. (Note
`assign_status` = great/good/acceptable/ambiguous is *assignment confidence*, not
completeness.)

**Provenance.** The rep is grey's — `cluster_name` IS the mmseqs `358M→20M`
representative. Our pipeline searched grey's reps; we never chose which member
represents a cluster. The fix is a rep *re-selection* (longest / highest
dpam_prob / non-partial member per cluster) upstream — grey's lever. Flagged to
grey (heads-up handled by RDS).

**Implication.** The Phase-3 "size 1.7×" factor is *not* a property of AFDB or of
PDB-vs-AFDB — it is which structures grey's pipeline chose as reps. Searching
full-length domains/members would roughly double SSE content for the affected
globular folds and shrink the size factor. The **redundancy 3.8×** effect (Phase
2) is unaffected and remains the real driver.

Durable materials: `~/work/prosmos_2026/rep_truncation/` (renders/, models/,
data/, README.md, render.pml, phase1_profile.py).

## Data pointers
`afdb_db/`, `ecod_db_pdb_exp/`, `s5_full_afdb/summary.tsv`,
`s5_full_pdb_exp/summary.tsv`; SSE count via `grep -oP '\.ssd\s+\K[0-9]+'`;
searchmatrix `searchMatrix/src/searchControl.h`.
afdb_200m schema (lotta/ecod_protein): `ecod_domain_range` (domain ranges, join
via `protein.uniprot_acc`→`protein.id`=`protein_id`), `cluster_summary_full_v4`,
`ecod_cluster_roster_p*`. grey parquets: `/home/grey/afdb.200m/summary/` +
`assignments/`.
