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

**At the fold level they hit a comparable NUMBER of T-groups (1,086 vs 1,360) —
but this "80%" is a COUNT ratio, NOT overlap; see the geometry-vs-fold section
below (real T-group overlap is only 606, Jaccard 0.33). What actually agrees is
the CELL (geometry) level.**

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

## Rep "truncation" — investigated and RETRACTED (Simpson's paradox); 2026-07-09

**Retraction.** An earlier version of this section (commit 41c9323) claimed the
AFDB searched reps are systematically truncated fragments (~half the length of
their cluster members; T-group ratios 0.26–0.55, worst for globular folds). **That
claim was a weighting artifact and is wrong.** Attempting to *quantify the fix*
(re-search with full-length reps) exposed it.

**What the wrong claim compared.** "REP median 56 vs MEMBER median 127" and the
T-group ratios compared **cluster-weighted reps** (one per cluster → dominated by
the *many small-domain clusters*) against **member-weighted members** (roster
rows → dominated by the *few huge clusters* with long domains). Different
populations → Simpson's paradox. Reps looked short only because small-domain
folds form many clusters while long-domain folds form few (large) clusters.

**Proper paired test (the correction).** For each cluster, rep length vs its own
members' median (n=2,649 sampled clusters, ≥2 members):
- **rep/member-median ratio = 0.99** (p25 0.96, p75 1.00).
- Clusters where rep < member-median: 52% (≈ chance). Clusters are length-
  *homogeneous*; the rep faithfully represents its cluster. **No within-cluster
  truncation.**

**The "fix" barely moves anything.** Re-selecting each rep to its cluster's
LONGEST member (random cluster sample, cluster-weighted like the real search set):
| cluster size | n | rep med | longest med | recovery |
|---|---:|---:|---:|---:|
| 2–3 | 1201 | 58 | 63 | 1.09× |
| 4–9 | 710 | 66 | 74 | 1.12× |
| 10–49 | 506 | 73 | 84 | 1.15× |
| 50–499 | 189 | 86 | 102 | 1.19× |
| 500+ | 45 | 114 | 137 | 1.20× |

Even the largest clusters recover only ~1.2×. So the A/B search was **not run** —
the input changes ~1.1× and the hit-rate recovery would be negligible. The
per-cluster-render examples (rep 82 aa vs a 179-aa / 17-SSE member) were **real
but unrepresentative** — cherry-picked from large *bimodal* clusters; typical
clusters are homogeneous. (Figures `s5_rep_truncation.png`, `cluster_violin.png`,
`rep_vs_member.png` are kept as the illustrative-but-atypical record.)

**Correct attribution (unchanged conclusion).** The AFDB searched reps ARE shorter
than the full domain population (median 66-res SSE-span / 5 SSE vs population
130 res). But that is the **dereplication effect**, not truncation: the search set
is one rep per cluster (cluster-weighted), and cluster-space is dominated by
small conserved-domain families. This is the same phenomenon as the redundancy
finding — reads through the cluster-weighting, not a rep defect. So:
- The Phase-3 **size 1.7×** stands as a real property of the dereplicated search
  set (small-domain-heavy cluster space), NOT a fixable rep artifact.
- The **redundancy 3.8×** (Phase 2) remains the dominant driver.
- **Nothing to flag to grey** — the reps are representative; the DPAM/clustering
  is fine. (The partial-domain flag `partial_domains.parquet` exists and is mildly
  enriched in globular reps ~8% vs 5%, but does not indicate a systematic problem.)

Durable materials: `~/work/prosmos_2026/rep_truncation/` (renders/, models/,
data/{fix_pairs,fix_lengths,fix_candidates}.tsv, README.md). Method note: any
rep-vs-member comparison MUST be paired within-cluster — never cluster-weighted
reps vs member-weighted members.

## Geometry vs fold: PDB and AFDB reach the same CELLS, diverge only in the rare
## fold tail (2026-07-09) — corrects the "80% recovery / subset" framing

**The correction.** "AFDB recovers 80% of PDB's folds" was a *count* ratio
(1,086/1,360), not overlap. At the **T-group** level the sets are NOT a subset:
| | count |
|---|---:|
| PDB hitting T-groups | 1,360 |
| AFDB hitting T-groups | 1,086 |
| **shared** | **606** (Jaccard **0.33**) |
| PDB-only | 754 · AFDB-only | 480 |

So *neither* is a subset of the other at the fold level. **What agrees is the
CELL (geometry) level** (binary occupancy: 1,832 both, 634 PDB-only, 279
AFDB-only — AFDB ~88% subset of PDB cells).

**The divergence is the rare-fold tail, ~orthogonal to ECOD class.** PDB-only
T-groups are rare (median 24 PDB hits vs 114 shared) and spread across 547
X-groups; AFDB-only rarer (median 10) across 352 X-groups (+ the X109 solenoid
surplus). Not a specific ECOD architecture — the marginal/rare folds each dataset
samples differently, plus AFDB's repeat families.

**Worked example — Ig (X-group 11), and 11.9.1.** The *core Ig fold* 11.1.1 hits
in BOTH, massively (PDB 7,666 / AFDB 732); 9/15 X11 T-groups are shared and carry
~all the hits. The 6 "PDB-only" X11 T-groups are rare sub-families (5 of them ≤23
PDB hits). So **Ig is NOT a PDB-only fold** — "X11 40% PDB-only" was T-group
*counting* of the rare tail. The one non-trivial PDB-only Ig sub-family, **11.9.1
(404 PDB, 0 AFDB)**, was probed directly: it occupies **46 S5 cells in PDB, and
all 46/46 are AFDB-occupied** (via other folds), each a promiscuous β-rich cell
shared with **7–257 other T-groups** (TIM 2002.1.1, kinase 206.1.1, P-loop
2004.1.1, core Ig 11.1.1…). **So 11.9.1 contributes ZERO PDB-only cells.**

**The definitive statement.** Even the folds PDB *uniquely* hits add **no new S5
cells** — they occupy geometry AFDB already reaches through other folds. So:
> Geometry (cells) agrees between PDB and AFDB; T-group membership diverges only
> in the rare/sub-family tail (a classification + sampling effect), and that
> divergence carries **no unique geometry** in either direction. Both datasets
> paint the same bounded 5-SSE geometric picture; AFDB just labels it with fewer,
> dereplicated folds.

Method note: compare fold repertoires at the **cell (geometry)** level, or coarser
than T-group (X/H-group), not raw T-group — T-group divergence is dominated by
rare-sub-family sampling and carries no geometry.

## Darkness is compact geometry S5 can't template (Case 2, not Case 1)

Why do most ≥5-SSE structures NOT hit? Checked the SSE contact graph of every
≥5-SSE pdb_exp structure (contact = interaction-matrix letter code), split by hit
status:
| pdb_exp ≥5-SSE | has compact 5-core (≥5 SSEs, each ≥2 contacts) | mean SSE degree |
|---|---:|---:|
| hitting (72,180) | 99.9% | 2.69 |
| **dark (373,658)** | **96.2%** | 2.39 |

Dark structures are **NOT sparse/extended** — 96% are compact globular units with
a 5-SSE core, barely less compact than hitters. They simply don't match the 198
2D-hex skeletons. **So S5 is blind to a large class of real *compact* 5-SSE
geometry** (not just non-motif structure). Consequence for the thesis: "saturated
/ no new topology" is scoped to the *narrow slice of compact 5-SSE space S5 can
template* — a minority of even the compact-domain space. This is the empirical,
quantitative form of the Lesk / "local, fixed-cardinality descriptor" caveat, and
it reinforces the two deck caveats (S5-is-local; AF-is-PDB-trained).

### Why a compact domain is dark: TOO connected, not disconnected (contact ceiling)
The hard question — "is there a matrix that *would* hit these that's just missing
from our set, or could no S5 ever hit them?" — has a precise, mechanistic answer.
These dark compact domains **do** form a connected 5-SSE motif — often a *complete*
one (every SSE contacts every other, Kₙ). The failure is the opposite of
disconnection: they are **denser than any skeleton can be**.
- **Contact ceiling.** Across all 198 S5 skeletons the densest has **7 of 10**
  possible 5-SSE contacts (dist: 4→12, 5→98, 6→28, 7→60 skeletons). **None is a
  5-clique** — and none can be: the 2D-hex lattice's max clique is 3 (same wall as
  the oracle-gap K₁,₄/C₅ non-realizability, [[project_s5_oracle_gap]]). So the
  198 are complete *for that model*; a denser template isn't "missing," it is
  **un-enumerable without leaving the planar-hex model**.
- **The topological wall is REAL but RARE (correction).** A *complete-graph* (Kₙ)
  domain matches no skeleton: every EEEEE skeleton has ≥1 require-non-contact among
  the mapped pairs, which a complete graph violates; the hand-picked gallery (all
  Kₙ) are genuine topological walls. **But this does NOT generalize to all dense β
  domains.** Earlier framing ("≥8-edge all-β ⇒ matches nothing, induced wall") was
  **too strong**: **58 of 198 EEEEE skeletons are multi-sheet and code their
  cross-sheet non-edges `X` (wildcard), not `-`.** So an 8–9-edge β graph *can*
  embed topologically into a multi-sheet skeleton — it is then dark on **geometry**
  (contact-type, strand register, handedness), not topology. The empirical **0 hits
  among 11,922 ≥8-edge 5-SSE all-β domains** (`figures/s5_contact_ceiling.png`)
  stands, but the *mechanism* is dominantly geometric; only the true complete-graph
  tail is a topological wall. See the population sieve below for the split.
- **α = monomorphism (geometry wall).** For H-typed pairs non-edges are coded `X`
  = wildcard, so extra contacts are tolerated; a compact all-α bundle's *contact
  pattern* is satisfiable. Dark for the same **geometry** reason — contact-type
  (v/u vs angle class), the 9 per-skeleton **handedness** triples, distance gates.

Bottom line for the thesis: darkness is **not** "no connected 5-SSE motif" and
**not** a fixable enumeration hole. It is **overwhelmingly geometric** (2D-hex
handedness/angle/contact-type rejects a topology that *is* representable) — see the
sieve; the pure-topology wall (Kₙ, denser than any skeleton) is a rare β-enriched
tail (~0.6% of the compact residual), not the main story.

### Full dark-sieve of all 4.92M searched reps (Tier-1 exact + Tier-2 proxy)
`sieve_dark.py` (one exact pass) + `proxy_d1d2.py` (50k topological-match sample).
Reconciled hit set via **15-char name truncation** (searchmatrix truncates output
names, dropping `_nD<n>` for 10-char UniProt accs; 939 ambiguous truncations →
2,272 sibling reps conservatively counted as hit). Mutually-exclusive buckets
(`figures/s5_dark_sieve.png`):

| bucket | count | % | meaning |
|---|---:|---:|---|
| A · <5 SSEs | 2,114,879 | 43.0% | too few SSEs to fill a 5-motif |
| B · ≥5 SSE, <5 pass length | 645,490 | 13.1% | SSEs too short (H≥8/E≥5 filter) |
| C · ≥5 pass, sparse | 54,528 | 1.1% | no connected 5-motif (Case-1 extended) |
| **D2 · compact, geometry-rejected** | **2,046,098** | **41.6%** | topology fits a skeleton; 2D-hex geometry doesn't |
| **D1 · compact, topology-absent** | **12,972** | **0.26%** | contact pattern fits NO skeleton (Kₙ / β tail) |
| HIT | 15,593 | 0.32% | — |

- **56% (A+B) is a size/length story**; only **1.1%** is genuinely extended (C).
- The compact residual (D, 42%) splits **D2 99.4% / D1 0.6%** (proxy, ±0.1% CI):
  darkness of compact domains is **geometric**, not topological. D1 is
  **β-enriched** (55% all-E vs 18% base) — the induced-non-contact/complete-graph
  tail, as predicted; D2 is composition-broad.
- Method note: the topological proxy embeds each typed skeleton's
  contact/`-`/`X` pattern (ignoring contact-type identity, handedness, distance —
  those are the "geometry" it isolates); validated on Kₙ all-β (→ 0 matches).

### Same sieve on PDB (pdb_exp) — is the gap AFDB-specific? No.
Ran the identical sieve on pdb_exp (`sieve_pdb.py`, 496,359 ECOD domains; note the
SSE regex must tolerate real chain ids `EP/HX/E3/EB` + glued/lettered residues —
AFDB was all chain `A`). Clean: buckets sum exactly, all 72,180 hitters in D, 0 in
A/B/C. Figure `figures/s5_sieve_pdb_vs_afdb.png`.

| bucket | PDB % | AFDB % |
|---|---:|---:|
| A <5 SSE | 10.2 | 43.0 |
| B SSEs too short | 6.9 | 13.1 |
| C sparse/extended | 0.9 | 1.1 |
| D2 compact · geometry-rejected | ~67.4 | ~41.6 |
| D1 compact · topology-absent | ~0.14 | ~0.26 |
| HIT | 14.5 | 0.32 |

**Three findings:**
1. **The mechanism is identical.** PDB compact-dark splits **D2 99.8% / D1 0.2%**
   (proxy, ±0.04) — the same as AFDB (99.35 / 0.65). Compact-domain darkness is
   geometry-rejection in *both* databases; **not an AFDB artifact.** The S5 matrix
   is geometrically narrow for real experimental PDB structures too.
2. **Composition differs by domain size.** AFDB is **56% too-small/short (A+B)**
   (dereplication keeps many small single-domain clusters); PDB only **17%** (big
   multi-SSE ECOD domains, median ~12 SSE). So PDB is dominated by the eligible
   compact set: **D-dark is 67.6% of all searched PDB vs 41.8% of AFDB** — PDB is
   proportionally *more* compact-but-dark.
3. **The residual gap is redundancy.** Within the eligible compact-motif set
   (bucket D), PDB hits **17.7%** vs AFDB **0.71%** (~25×). Since D controls for
   size/length/sparsity, this residual is the redundant-census vs dereplicated-
   sample effect: PDB crystallizes the hitting geometries many times over.

Bottom line: PDB shows the **same** stratification gap and the **same** geometric
darkness mechanism; it differs only in entry-size composition (fewer sub-floor
domains) and per-structure redundancy (25× within-eligible). The "S5 is a narrow
geometric descriptor, and most real compact domains have representable topology but
non-lattice geometry" conclusion is database-independent.

**Length-filter control (important).** The S5 queries hard-filter every position
by SSE length (`length` line: **H ≥ 8, E ≥ 5** residues, uniform across all 32
typings). So darkness could trivially be "too-few long-enough SSEs." Controlled
for: the worked set was re-derived counting **only** length-passing SSEs
(`find_dark_compact_v2.py` → **1,913,192** length-clean dark compact ≥5-SSE
domains, vs 2.46M before the filter — ~22% of v1 were length-confounded). The 20
gallery examples each have **≥5 SSEs that clear the length filter**, all in a
complete contact graph — so their darkness is geometry/topology, not sub-threshold
SSEs. Set: `~/work/prosmos_2026/dark_gallery/` (`dark_compact_montage_v2.png`,
`dark20b_annotated.tsv` with per-SSE lengths); assigned folds among them: winged
helix (101.1.2), LysM (101.15.1), HET-s left-handed β-helix (208.5.1), pectin-lyase
right-handed β-helix (207.2.1). (v1 `dark_compact_montage.png` had 8/20 failing the
length filter — superseded.)

## Data pointers
`afdb_db/`, `ecod_db_pdb_exp/`, `s5_full_afdb/summary.tsv`,
`s5_full_pdb_exp/summary.tsv`; SSE count via `grep -oP '\.ssd\s+\K[0-9]+'`;
searchmatrix `searchMatrix/src/searchControl.h`.
afdb_200m schema (lotta/ecod_protein): `ecod_domain_range` (domain ranges, join
via `protein.uniprot_acc`→`protein.id`=`protein_id`), `cluster_summary_full_v4`,
`ecod_cluster_roster_p*`. grey parquets: `/home/grey/afdb.200m/summary/` +
`assignments/`.
