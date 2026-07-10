---
title: "Where the S5 negative space lands — the full landscape, the novelty gauntlet, and one methodological trap"
author: "R. Dustin Schaeffer"
date: 2026-07-09
marp: true
paginate: true
---

# Where the S5 negative space lands

## Follow-up: from an 800-query hole to a saturation measurement

*The first deck ended with the 800-query whole-skeleton negspace holding at zero across 5.4M structures, and the full 6,336-cell sweeps in flight. They are now complete. We then ran every fold-novelty probe we could reach — AFDB dark matter, TED novel folds, and a multi-filter annotation funnel — through the same S5 instrument. All land on known geometry. One probe looked like it broke the pattern; it turned out to be a methodological trap worth stating carefully, because it is the objection a skeptic will raise.*

<!-- Notes: This is the closing arc. Deck 1 = enum + negspace + four bugs. Deck 2 = full landscape done, the negative reframed as a positive claim, the novelty gauntlet, and the one controlled experiment that seals it. Audience: lab + Jimin, who now has the interactive inspector. -->

---

# Recap — where deck 1 left off

- **Enumeration**: 198 S5 skeletons × 32 H/E typings = **6,336 queries** from the Chitturi 2016 2D-hex-lattice model.
- **Whole-skeleton negspace**: 25 skeletons (800 queries) return zero across *every* typing.
- **Scale-up test**: 800 queries × (19k curated + 496k experimental PDB + 4.9M AFDB) = **2,400 sweeps, 0 lit-up**. Not one negspace motif is realized at 5.4M structures.
- **Cost**: four latent `searchmatrix` buffer/validation bugs, all fixed.
- **Then in flight**: the full 6,336-cell sweeps + the interactive inspector.

> Deck 2 picks up exactly here: the sweeps finished, and the question became *how far does the negative generalize?*

<!-- Notes: One-slide bridge so this deck stands alone. The 800 = whole-skeleton-zero (strict). The full landscape introduces the per-CELL negspace, a superset. Keep those two numbers distinct on the next slide. -->

---

# The full landscape is in — four databases, one matrix

| Dataset | Searched | Occupied cells | Kind |
|---|---:|---:|---|
| ECOD manual reps | 19,015 | 2,328 | experimental |
| Full experimental PDB (ECOD domains) | 496,359 | 2,466 | experimental |
| AFDB non-singleton reps — ECOD-assigned | 4,921,931 | 2,074 | predicted |
| AFDB non-singleton reps — unassigned | 787,509 | 786 | predicted |

- **3,380 of 6,336 cells stay empty** across all four (the per-cell negspace; a superset of the 800 whole-skeleton holes).
- All browsable in the inspector: heatmaps, per-cell exemplars rolled up by ECOD T-group, real per-pair contact geometry, all-pairs compositional enrichment.

<!-- Notes: 3,380 is the per-cell union-of-absence; 800 is the strict whole-skeleton version. Occupied union ≈ 2,956 cells. The inspector is the deliverable Jimin gets — leda:3001/docs. -->

---

# The reframe — the negative *is* the result

- S5 is an **assignment-free instrument**: pure SSE geometry, fully enumerable, so it comes with a **denominator** (6,336). That makes "what's missing" a measurable quantity, not an anecdote.
- A novelty test built on CATH/ECOD begs the question (novel = "we didn't classify it"). S5 does not.
- **The claim**: the observable 5-SSE fold repertoire is effectively *closed*. AlphaFold multiplies the **membership** of known folds by orders of magnitude but adds essentially **no new local topology**.

> Sequence space is exploding; it maps many-to-one onto a bounded, nearly-filled fold codomain that stopped growing orders of magnitude of sequence ago. AF is the lens that lets you see the collapse.

<!-- Notes: This is the thesis. "Saturation measurement" not "we failed to find novelty." The many-to-one framing is the paper's spine and reframes ECOD's value: AF maps divergent members back to a reference — good news for classification, not a threat. -->

---

# AFDB coverage: per-structure misleads, fold level agrees

![w:900](../docs/figures/s5_foldlevel.png)

- **The trap**: per-structure hit rate looks like AFDB found almost nothing — PDB-exp 18.8% vs AFDB 0.55% among ≥5-SSE structures (**34×**). But that compares a *redundant crystallization census* to a *dereplicated, smaller* diversity sample — not a fair comparison.
- **The fair view (fold level)**: distinct S5-hitting ECOD T-groups — PDB **1,360** vs AFDB **1,086**. **AFDB recovers 80% of the fold diversity all of experimental PDB does**, from dereplicated reps.
- **The 34× decomposes**: redundancy **3.8×** (PDB ~53 entries/fold vs AFDB ~14) × size **1.7×** (AFDB search set is dereplicated *cluster reps*, and cluster-space is dominated by small-domain families → median 5 vs 12 SSE) × denominator composition (AFDB "capable" = the whole dereplicated proteome, mostly small/dark non-fold clusters). **Not** incompleteness — the sweep is complete (6,336 queries, 4.92M entries, verified).

<!-- Notes: This replaces the rarefaction curve, which confused people. Audit chain that kills the "did we finish the search" read: afdb DB has 4,921,931 entries (grep -c .ssd), sweep 6,336 queries all rc=0, searchmatrix streams the DB (no entry cap), a traced hit is correct. Redundancy is the driver (3.8x): dereplicating pdb_exp by T-group, 1,360 folds carry 52.9 entries each. SIZE (1.7x) IS A DEREPLICATION EFFECT, NOT REP TRUNCATION: reps faithfully represent their clusters (paired within-cluster rep/member-median length = 0.99); the search set is short only because it is cluster-weighted and cluster-space is small-domain-heavy. DO NOT say reps are truncated fragments or that grey's rep-selection is buggy -- an earlier truncation claim was a Simpson's-paradox artifact (cluster-weighted reps vs member-weighted members), retracted; re-selecting full-length reps recovers only ~1.1x. Interactive: inspect fold-nonredundant heatmaps pdb_exp_nr / afdb_assigned_nr + fold-vs-fold compare (1,815 both-occupied, only 30/2,725 significant). Full trail incl. retraction: enum/docs/dataset_analysis_plan.md. -->

---

# Cell by cell, PDB and AFDB agree at the fold level

![w:520](../docs/figures/s5_foldcompare.png)

- Redundancy-controlled diff: **log₂(fold-share PDB / AFDB)** per S5 cell, same skeleton×typing grid. **Pale = agree.**
- 1,815 cells occupied in both; **only 30 of 2,725 differentially enriched** (q<0.05, circled). A faint PDB lean (it carries ~25% more total folds) — but **no cell class is AFDB's alone, and none PDB's alone at scale**.
- Two independent samples of protein space — a redundant crystallographic census and a dereplicated predicted proteome — occupy the **same motif cells in the same proportions**. That convergence is the real evidence for a bounded repertoire.
- The **AFDB-enriched cells are all HHHHH** (all-α) — and it's a confound, not new diversity: the *same* α-solenoid/repeat families (ARM, Ankyrin, HEAT — X-group 109) + TM-helical (GPCR, MFS) **tile overlapping 5-helix windows across many skeletons at once**, so one family inflates every all-α cell. Composition (eukaryote repeats/TM) × S5-blind-to-helical-identity × solenoid tiling — not reachable topology.

<!-- Notes: The interactive twin is /compare/pdb_exp_nr/afdb_assigned_nr in the inspector. Grey = neither occupied. Slight red skew = PDB's 1,360 vs AFDB's 1,086 folds (mild, compositional). The 30 circled cells are the only significant differences after BH correction — trivially few for 2,725 cells. All 27 afdb-enriched cells are HHHHH (ty0); the SAME T-groups (109.4.1 ARM, 109.3.1 Ankyrin, 5001 GPCR, 5050 MFS) recur across all of them -> it's one phenomenon (α-solenoids/TM helices tiling across skeletons + eukaryotic repeat abundance), double-counted across cells, NOT 27 independent diversity gains. So even AFDB's only fold-surplus is an artifact stack, reinforcing "no new reachable topology." This is the diff heatmap the per-structure compare could never show cleanly (that one was dominated by the ~53x vs ~14x redundancy offset). -->

---

# Why the fold repertoire looks closed — and the honest caveat

- The fold-level agreement (AFDB ≈ 80% of PDB's S5-hitting folds, from a totally
  independent, dereplicated sample) is the positive statement: **the two views
  of protein space converge on the same bounded set of local motifs.**
- AFDB's *lower* count (1,086 vs 1,360) is shallower sampling of the rarer folds,
  **not** a different repertoire — it adds essentially no cells PDB lacks.
- **Caveat kept in view**: "fold" here = ECOD T-group as the dereplication unit;
  the per-structure rate is genuinely uninformative across datasets of different
  redundancy + size, so all cross-dataset claims are made at the fold level with
  the completeness audit attached.

<!-- Notes: Optional depth slide; can cut for time. The point: once you dereplicate, PDB and AFDB agree, which is the real evidence for a bounded repertoire — stronger than the rarefaction ever was. Keep the T-group-as-unit caveat so nobody reads "fold" as sequence-cluster. -->

---

---

# The novelty gauntlet — every probe we could reach

Each pushes a differently-defined "new fold" set through the *same* S5 instrument.

| Probe | Set (how "novel" is defined) | Result at S5 |
|---|---|---|
| Whole-skeleton negspace | 800 lattice motifs, never seen | 0 / 2,400 sweeps at 5.4M |
| AFDB dark matter | unassigned AFDB reps (no ECOD) | under-hits; 786 cells, 96% PDB-known |
| TED novel folds | 7,427 curated, **no CATH match** | 89% PDB-occupied; **1** robust novel |
| Annotation funnel (v3) | 4,201 dark clusters, multi-filter | whole-model 40 → **domain-chopped 0** |

- Four independent definitions of novelty — geometric, dark-sequence, no-CATH, and annotation-filtered.
- **All land on known S5 geometry.** The one apparent exception (v3) is a methodological artifact — next slides.

<!-- Notes: This table IS the deliverable the user asked to enumerate. Note the definitions are genuinely different and non-circular w.r.t. S5 (which is classification-independent). TED is the strongest single point: most-curated fold-novelty set in existence, still 89% geometrically known. -->

---

# Probes 1–2 — dark matter and TED novel folds

**AFDB dark matter** (sequence-dark ≠ fold-dark):
- ~16% of the DB but only ~4% of hits — it *under*-hits S5.
- 786 occupied cells, 96% already PDB-known; ~16 novel-vs-all-known, all singletons.

**TED novel folds** (defined by *no CATH match* — a classification-based novelty claim):
- 89% land on PDB-occupied cells.
- Exactly **1** robust (multi-structure) novel cell: skeleton 118 / EEHHE.
- Non-circular: S5 knows nothing about CATH, yet agrees the geometry is known.

> Even the most-curated fold-novelty set on the planet is, at the 5-SSE level, geometry we have already seen.

<!-- Notes: The spine finding. Sequence divergence and classification-novelty both dissociate from geometric novelty. The TED "1 robust" is the honest residual — reported, not hidden. -->

---

# Probe 3 — the annotation funnel, and the trap

- `~/work/afdb_200M` GENUINE funnel: dark clusters → no-Pfam → no-CATH → no-InterPro → pLDDT≥70 → low-decay → not-prophage → ≥10 members. **4,201** "genuine novelty" candidate domains.
- Searched as **whole multi-domain models**, they appear to hit **40 robust never-seen cells** (top: HEEHE, 6 independent proteins) — the *only* probe to break the pattern.
- **The catch**: those motifs draw their five SSEs from *scattered positions across a whole protein* (median span 234 res vs 153 for known-cell hits) — i.e. **across domain boundaries**. An arrangement a domain-chopped reference search *structurally cannot produce*.

> The reference databases were searched as **chopped domains**; v3 was searched as **whole models**. That granularity mismatch — not novelty — manufactures never-seen cells.

<!-- Notes: This is the honest "we almost fooled ourselves" slide. The 40 robust cells were real matches, but cross-domain. The right comparison must match input granularity. This is exactly the objection a ProSMoS partisan raises, so we ran the controlled experiment. -->

---

# The sealing test — a controlled A/B

Chop each candidate to **its own committed `_D<n>` domain boundary** (the pipeline's own decomposition — no arbitrary cutoff), rebuild with the **identical PALSSE + generateMatrix pipeline**. Only one variable changes: chopping.

| | whole-model v3 | domain-chopped v3 |
|---|---:|---:|
| occupied cells | 2,063 | 982 |
| total hits | 96,191 | 16,981 |
| never-seen cells | 150 | 25 |
| **never-seen w/ ≥2 proteins** | **40** | **0** |

- The 40 robust "novel" cells → **0**. The 25 residual never-seen cells are all singletons (noise floor).
- **Lesson for the paper**: whole-model vs domain-chopped input must be matched, or cross-domain SSE assembly masquerades as fold novelty.

<!-- Notes: Unimpeachable because (a) boundaries are the pipeline's own, not a cutoff we chose, and (b) it's a within-v3 A/B holding SSE-method fixed. This retires the "you cherry-picked the compact filter" objection. Builds at ~/work/prosmos_2026/{v3_dom,s5_v3_dom}. -->

---

# When to trust S5 — the instrument's own limits

![w:620](../docs/figures/s5_geom_determinism.png)

- **Composition-graded fidelity**: within-cell 5-SSE geometry tightens monotonically with strand count (~−1.9 Å/strand); a **1-strand penalty** — all-α / single-strand cells are grab-bags, β-registration makes a cell determinate.
- **It's baked into the query, not just the structures**: an `EEEEE` query carries ~700 hard `-` (non-contact) constraints + ~1,000 sheet-pairing codes (`c`/`t`) across the 198 skeletons; `HHHHH` is ~854 `X` wildcards with **zero** hard constraints. Same-sheet non-adjacent strands *cannot* touch (register-separated) → `-`; non-adjacent helices are unconstrained → `X`. So all-α cells are intrinsically *looser* queries and match promiscuously — the fidelity gradient is half structural, half encoding.
- **Blind to transmembrane**: GPCR (5001.1.1) and MFS (5050.1.1) — maximally different TM folds — both land in generic `HHHHH` cells.
- **Tiling inflates abundance**: ARM repeats (109.4.1) — one domain hits up to 123 cells. Abundant ≠ novel; controlled by compositional (share-of-hits) normalization.

<!-- Notes: Method introspection = credibility. We state where the descriptor is strong (β-rich) and where it's a grab-bag (all-α). The encoding bullet: the -/X asymmetry is NOT our invention. Vocabulary is paper §1.1.1 verbatim ("Non adjacent SSEs can either interact optionally (X) or not interact at all (-)"); the assignment rule (same-sheet non-adj strand -> '-', else -> 'X') exactly reproduces the CG-2012 oracle -- verified 0 violations across 2,807 IA-S5.txt records (2,766 '-', 13,132 'X'). It's physically forced: same-sheet non-register-adjacent strands are separated by intervening strands, can't contact -> hard '-'; helices/cross-sheet unconstrained -> wildcard 'X'. So all-helix cells are looser QUERIES, reinforcing the determinism/promiscuity gradient from the query side. TM-invisibility and ARM-tiling are the two cautionary examples; compare view is compositional, not per-structure rate. -->

---

# Where it lands

- **"ECOD is completed" — bounded to three regimes**: (1) the common core is discovery-complete; (2) the frontier moved from *discovery* to *membership* (AF = remote-homology-at-scale); (3) the tail is curation-hard, not discovery-rich (repeats, symmetry, ARM).
- **Bottom line**: AFDB does **not** robustly reach S5 cells PDB can't — it's a ~80% *subset* of PDB's folds sampled far more densely, and every novelty probe lands on known geometry. Two caveats bound this:
- **Caveat 1 — S5 is a *local, fixed-cardinality* descriptor (Lesk)**: it measures 5-SSE *local motif* geometry — blind to global topology, transmembrane arrangement, domain-scale novelty. "No new S5 cells" ≠ "no novel global folds"; an empty cell = "never observed here," not "unfoldable."
- **Caveat 2 — AlphaFold was trained on the PDB**: AFDB agreeing with PDB on what's reachable is *partly a training prior*. The data can't cleanly separate "the real proteome is saturated" from "AF only renders PDB-like geometry." Same practical result (broader prediction surfaces no new cells), genuinely ambiguous mechanism.
- **Still open**: are the 3,380 empty cells sterically reachable-but-untaken, or lattice over-generation? Resolution needs constructive 3D backbone modeling (RFdiffusion / Rosetta) on the negspace.

<!-- Notes: The careful version of "we're done." Three regimes so nobody hears "protein folding is solved." The two caveats are load-bearing for referees: (1) S5 is local so we only claim local-motif saturation, not global-fold closure; (2) AF is PDB-trained so the AFDB-doesn't-reach-new result is partly a prior, not independent evidence — state both explicitly and don't overclaim "fold space is closed." The reachability question is the same open item as deck 1 — constructive-modeling territory, not more searching. -->

---

# Recap

- **Full landscape done**: 6,336 cells × 4 databases; 3,380 stay empty; all in the interactive inspector.
- **The negative reframed**: a saturation measurement with an assignment-free, enumerable instrument — the observable 5-SSE repertoire is effectively closed.
- **PDB ≈ AFDB at the fold level**: per-structure hit rate misleads (34×, redundancy×size); dereplicated by fold, AFDB recovers **80%** of PDB's S5-hitting folds (1,086 vs 1,360) — coverage is redundancy, not incompleteness.
- **Novelty gauntlet**: four independent novelty definitions (lattice, dark-sequence, no-CATH, annotation-funnel) → all land on known geometry.
- **The one trap**: whole-model search fabricated 40 "novel" cells; a controlled domain-chopped re-search collapsed them to 0. Input granularity, not topology.
- **Deliverables**: inspector (`leda:3001/docs`), sealed capstone, full session log.
- **Open**: geometric reachability of the 3,380 — constructive modeling, not more search.

<!-- Notes: One-slide close. Order: what finished, the claim, the gauntlet, the trap+seal, deliverables, the single remaining open question. -->

---

# Backup — apples-to-apples rarefaction (per T-group)

![w:720](../docs/figures/s5_foldrarefaction.png)

- The honest replacement for the retired per-structure rarefaction: **x-axis is distinct ECOD T-groups, not structures** → removes the redundancy confound (PDB census vs AFDB dereplicated).
- **PDB and AFDB track within ~15%** — AFDB T-groups are about as cell-productive as PDB's (1.91 vs 1.81 cells/T-group). No dataset samples matrix-space dramatically better per T-group.
- Band = IQR over 40 random T-group orderings (tight → order-robust). Residual caveats *partially cancel*: PDB samples each T-group ~6–33× deeper (union → more cells); solenoid tiling inflates cells/T-group in both.

<!-- Notes: Backup for the reachability/coverage Q&A. Corrected form of the parked rarefaction fork. Line = mean over 40 random T-group-add orders; cloud = 25-75th-percentile band across those orders (sampling-order variability), NOT a confidence interval on a fit. Used "T-group" not "fold" — more precise for this ECOD-literate audience. manual_pdb deliberately excluded: different T-group universe (ecod_rep all classes 2,176 vs af2_pdb subset 1,360) + ~1 struct/T-group, not on the same footing. Data: pdb_exp_build / afdb_build s5_promiscuity.json; script ~/work/prosmos_2026/plot_foldrarefaction.py. -->
