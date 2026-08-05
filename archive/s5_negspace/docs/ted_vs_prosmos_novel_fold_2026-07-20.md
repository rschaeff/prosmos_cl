# TED's "candidate novel fold" detection vs. the ProSMoS S5 grid

Local sources (verified against the PDF text, not the summary):
`~/work/ted/science.adq4946.pdf` (Lau et al., *Science* 386, eadq4946, 2024;
TED) and `~/work/ted/gkae1087.pdf` (Waman et al., *NAR* 53, D348, 2025; CATH
v4.4). All TED numbers below are quoted verbatim from the Science methods.

The short answer: **the two methods do not measure the same thing, and that is
exactly why they disagree.** TED asks *"is this whole domain globally unlike
anything in the structure libraries?"* — a similarity-exclusion question answered
by superposition. The S5 grid asks *"does this domain contain a 5-SSE local
arrangement, and is that arrangement one the PDB never realises?"* — a
vocabulary-membership question answered by discrete pattern matching. A domain
can be a resounding "yes" to the first and "no" to the second, and TED's own
novel folds are precisely those domains.

---

## 1. TED's novel-fold funnel (verbatim from the Science methods)

TED classifies 324M domains, then treats "novel fold" as the **residue that
survives an exclusion cascade** applied to the *unlabeled* clusters:

| Stage | Filter | Survivors |
|---|---|---|
| Unlabeled after Foldseek + Merizo-search vs CATH SSG5 | (no CATH H/T label) | **41,879,858 clusters** |
| Globularity | normalized Rg < 0.356 **and** packing density > 10.333 (5th-percentile cutoffs from H-labeled domains) | — |
| + SSE count | **≥ 6 SSEs** (helices+strands, STRIDE) | 13,820,550 |
| + pLDDT80 | mean pLDDT of top-80% residues **≥ 90** | 8,612,318 clusters (19,816,697 domains) |
| Exhaustive structure search | Foldseek easy-search TM-align mode vs **PDB, CATH, ECOD, SCOPe**; eliminate any match with TM-align > 0.56 & cov > 60% | 240,674 candidate clusters |
| High-symmetry split | SymD Z > 9 removed as *repeat* set (6,433 clusters) | low-symmetry remainder |
| Chopping-quality gate | Foldclass variant trained on bad choppings | — |
| Final TM-align pass | direct exhaustive TM-align, TM > 0.5 & cov > 60% | **24,653 → 7,427 clusters** |

TED's own definition of the 7,427: *"well-folded domains with no significant
structural similarity to current PDB chains (as of Dec 2023) or current
structural domain libraries (as of Feb 2024)."* Novelty is then **ranked**, not
thresholded — a density-based anomaly score (mean Euclidean distance to the 50
nearest CATH/ECOD/SCOPe Foldclass embeddings), because "there is no exact
boundary between novel and just highly divergent examples of known folds."

Key architectural point: TED's unit is the **whole domain**, its comparison is
**global superposition**, and its reference is **the union of all known
structures**. Novel = "superposes onto nothing known."

---

## 2. The S5 grid's detection process

| Axis | S5 grid |
|---|---|
| Unit of novelty | a **5-SSE local topology** (one grid cell = skeleton × H/E typing) |
| Reference | the **enumerated** hex-lattice skeleton set (198 graph / 140 geom at n=5) × which cells PDB domains occupy |
| Comparison | discrete **SSE-adjacency pattern match** (ProSMoS BFS + geometric constraints), not superposition |
| "Novel" = | a cell **lit by an AFDB domain but by no PDB domain** (an "AFDB-only cell") |
| Detection logic | **positive / enumerative** — every possible motif exists a priori; we ask which are realised |
| Quality gating | none — no globularity, pLDDT, symmetry, or SSE-count filter (n is fixed at 5) |
| Output | a 198×32 grid of lit/dark cells, per database |

Novel here is **inclusion in an a-priori alphabet cell that one database occupies
and the other does not**. It never asks whether a whole domain is unlike known
structures; it asks whether a *local motif* falls in an unoccupied slot of a
finite, pre-enumerated table.

---

## 3. Side-by-side

| | **TED** | **ProSMoS S5 grid** |
|---|---|---|
| Question | Is this **domain** globally unlike everything known? | Does this domain contain a **5-SSE motif** the PDB never realises? |
| Scale of novelty | whole fold (global architecture) | local topology (5 SSEs) |
| Logic | **exclusion** (survives a similarity cascade) | **enumeration** (fills an a-priori cell) |
| Comparator | Foldseek + TM-align superposition | discrete graph/geometric pattern match |
| Reference universe | all PDB + CATH + ECOD + SCOPe structures | the enumerated hex-grid + PDB's occupancy of it |
| Sensitive to | global arrangement, symmetry, permutation, size | 5-SSE adjacency only; blind to everything above n=5 |
| Quality filters | heavy (Rg, packing, pLDDT80≥90, SSE≥6, SymD, chopping) | none |
| Handles symmetry | explicitly (SymD Z>9 → separate repeat set) | invisible — a 40-strand β-propeller and a 5-strand sheet can light the same cell |
| Novelty is | ranked continuously (anomaly score) | binary per cell (lit/dark in each DB) |
| Result | 7,427 putative novel-fold clusters | 619 raw AFDB-only cells → ~90 depth-robust |
| Failure mode | divergent-but-known folds leak in ("no exact boundary") | a globally novel fold made of common local motifs is **missed by construction** |

---

## 4. Why they must disagree — and the number that proves it

The two axes of novelty are almost orthogonal. A fold is globally novel (TED)
usually because of **symmetry, circular permutation, or a new *arrangement* of
familiar units** — a 13-strand barrel, an 11-bladed propeller, an extruded
solenoid. But those same units are locally ordinary: five strands from a
13-strand barrel look like five strands from anything else. So a TED-novel domain
almost always lights a *PDB-occupied* S5 cell.

We measured this directly (recorded in `project_s5_instrument_ceiling`):

- **99.66%** of a 28,902-domain CATH-unassigned (TED-novel-type) set land in
  **PDB-occupied** S5 cells. TED's whole discovery class is, at the 5-SSE scale,
  already lit by the PDB.
- The **β-flower** — a *certified* TED/Pereira novel fold — run through the exact
  S5 pipeline lights cell `sk132/ty31`, which already holds **366 PDB folds**.
  Its novelty is entirely global; locally it is saturated space.
- Inversely, **52.5%** of the AFDB domains *we* call lit are TED-novel. So
  "CATH-unassigned" describes half of predicted structure space — a coverage
  statement about a PDB-derived classification, not a discovery signal — and our
  "AFDB-only cells" are not the same objects as TED's novel folds either.

**Neither instrument is a superset of the other.** TED finds domains the grid
calls ordinary (global novelty in local-common motifs); the grid flags cells TED
would never surface (a rare 5-SSE typing realised only in predicted space, inside
domains that TM-align happily matches to something known). They are answering
different questions with different reference frames.

---

## 5. What TED does that the grid structurally cannot

1. **Global architecture.** TED superposes whole domains; the grid never sees
   past 5 SSEs. Every propeller, solenoid, and barrel novelty TED reports is
   above the grid's ceiling.
2. **Symmetry as novelty.** TED's SymD split *is* the point for a large share of
   its hits; the grid is symmetry-blind.
3. **Quality-gated confidence.** TED's pLDDT80/globularity/chopping filters mean
   its 7,427 are curated well-folded domains. The grid applies no such gate — an
   AFDB-only cell can be lit by a poorly-predicted or mis-chopped model, which is
   why depth-correction knocks 619 → ~90.
4. **A continuous novelty rank.** TED can say *how* novel; the grid gives a
   binary per cell.

## 6. What the grid does that TED cannot

1. **Enumerative completeness at its scale.** The grid knows the *entire* set of
   possible 5-SSE topologies a priori, so "dark" is a statement about an
   exhaustive alphabet, not about "nothing matched in the libraries I happened to
   search." TED can never say a fold is *impossible*; the grid's impossible-mask
   can.
2. **Provenance-symmetric comparison.** The grid runs PDB and AFDB through one
   pipeline and compares cell occupancy directly; TED's reference *is* the PDB
   libraries, so it cannot ask "what does the PDB have that AFDB doesn't."
3. **Vocabulary vs. sampling.** The grid separates "a topology the PDB never
   realises" from "a fold nobody has crystallised." TED's exclusion cascade
   conflates the two into "novel."

---

## 7. The one-line takeaway for the deck

> TED and the S5 grid define novelty at different scales, so they see different
> things. TED detects **globally dissimilar whole domains** by superposition
> against all known structures — and those domains are, 99.66% of the time, built
> from **locally common 5-SSE motifs the PDB already realises**. The S5 grid is
> not a fold-discovery instrument and was never going to reproduce TED's list:
> **a new fold is not a new 5-SSE topology.** The grid's contribution is
> orthogonal — an *exhaustive*, provenance-symmetric census of local topological
> vocabulary, where TED offers a *filtered, ranked* census of global architecture.
