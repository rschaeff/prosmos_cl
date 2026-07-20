# AFDB-only S5 cells: 619 naive → 90 after the depth control

**Deliverable premise (RS's people):** concrete examples of cells lit only by
predicted (AFDB) structures, rendered with adjacent-motif context. RS flagged the
drift where "novel cell" gets read as "novel topology." This is the depth control
that separates the two.

Data: 100%-complete sweeps. AFDB `grid_afdb_nT_rebuilt.npy`, PDB `grid_pdb_nT.npy`
(400/400, 295,360 lit), `impossible_mask.npy` (640 unsatisfiable queries excluded).

## The naive list

At true 100% PDB, fold-level, satisfiable cells only:

| | |
|---|---|
| AFDB lit | 3,624 |
| PDB lit | 3,173 |
| both | 3,005 |
| **AFDB-only** | **619** |
| PDB-only | 168 |

Trajectory 720 (68% PDB) → 630 (94.5%) → 619 (100%): flattened, so not merely PDB
undersampling. Robustness: 313 singletons (1 fold), 234 at 2-4 folds, 72 at 5-19.

## The decisive control: size-matched rarefaction

AFDB searched 4.92M records, PDB 496K — **10x depth** (4x in lit structures: 1.19M
vs 295K). Subsample AFDB structures down, recount folds/cell, recount AFDB-only vs
the FULL PDB grid. Sanity gate: full rebuild from the (cell,record) pairs
reproduced `grid_afdb_nT_rebuilt` at **100.0% pattern agreement, r=1.0000**.

| AFDB depth | AFDB-lit recs | AFDB-only | PDB-only |
|---|---|---|---|
| 100% | 1,191,270 | **619** | 168 |
| 75% | 893,452 | 543 | 217 |
| 50% | 595,635 | 426 | 302 |
| **25%** (≈PDB's 295K lit) | 297,817 | **284** | **501** |
| 10% (searched-matched) | 120,135 | 148 | 842 |

**At equal structure depth the asymmetry REVERSES.** PDB's lit depth (295K) ≈ AFDB
at 25%, where AFDB-only (284) < PDB-only (501). The crossover is between 25% and
50% depth — above PDB's actual depth. So the 619 excess is **overwhelmingly a depth
artifact**: AFDB looks richer only because it sampled ~4-10x deeper. Matched fairly,
predicted structure space is not topologically richer than experimental; if
anything the reverse.

Third independent confirmation of the same conclusion:
[[project_negspace_finding]] (Ig 0.01% dark to PDB), the unannotated eligibility
match (56.2% vs 57.1%), and now this.

## The defensible deliverable: 90 robust cells

Persistence over 20 matched-depth (25%) draws — a cell is ROBUST if it stays
AFDB-only in ≥80% of draws:

| | cells |
|---|---|
| **ROBUST (≥80%)** | **90** — survive the depth control; lit only by predicted structures at PDB-equal depth |
| middle | 407 |
| FRAGILE (≤20%) | 122 — depth artifacts |

The 90 are the honest basis for examples. Top-persistent (all 20/20, with AFDB fold
count): sk16/ty12 (17), sk91/ty28 (15), sk164/ty3 (13), sk171/ty7 (11), sk57/ty24
(10), sk16/ty19 (10), sk158/ty27 (10)… — multi-fold, not singletons.

Artifacts: `afdb_only_mask.npy`, `afdb_persist.pkl` (persist dict + robust/fragile
lists), `afdb_rare/*.pkl` (8.71M cell-record pairs), `afdb_rarefy.py`.

## Per-topology significance — and why Fisher alone misleads (added 2026-07-18)

Requested (Qian): Fisher's exact per topology, heatmap it, then focus on cells
where one DB is 0 and the other is significant. Done — `s5_grid/fisher_cells.py`,
figure `fisher_panels.png`. Fold counts, not structure counts (structure-level
would treat ~30x within-fold redundancy as independent evidence): per satisfiable
cell, 2x2 = [[a, NA-a], [p, NP-p]] with NA=3,656 / NP=3,935 folds, BH-corrected.

| | cells | Fisher q<0.05 | depth-corrected q<0.05 |
|---|---|---|---|
| all tested | 3,792 | **1,457 (38%)** | 579 |
| PDB=0, AFDB>0 | 619 | **42** | **5** (4 significant under both) |
| AFDB=0, PDB>0 | 168 | **0** | **0** |

**Two problems, the second fatal for the zero-cells.**
1. Fisher calls 38% of cells significant. The null is wrong: counts are not
   independent draws (one domain hits many cells; domains in a fold are near-identical).
2. **The reverse direction can never produce a hit.** PDB-only cells carry 1-5
   folds, so [[0,3656],[4,3931]] cannot reach significance after BH; AFDB-only
   cells reach 17 folds *because AFDB was sampled ~4-10x deeper*. "PDB=0 & AFDB
   significant" is findable; "AFDB=0 & PDB significant" is structurally impossible.
   The asymmetry is sampling depth wearing a p-value. So the empty reverse list is
   a statement about power, NOT evidence the PDB lacks unique topologies.

**The replacement, still a real p-value.** Thin AFDB to PDB's structure depth
(q = 295,360/1,191,270 = 0.248, lit-structure ratio). A fold with k member
structures in a cell survives with probability 1-(1-q)^k, so the thinned fold
count is **Poisson-binomial** over that cell's folds — exact, no resampling,
BH-correctable. Compare PDB's observed count against that distribution.

**The four cells significant under BOTH** (`fisher_zero_cells.json`):

| cell | AFDB folds | q(Fisher) | q(depth) |
|---|---|---|---|
| sk16 ty12 | 17 | 2.3e-05 | 7.5e-05 |
| sk171 ty7 | 11 | 1.3e-03 | 2.2e-02 |
| sk158 ty27 | 10 | 2.6e-03 | 6.2e-03 |
| sk34 ty0 | 8 | 9.5e-03 | 2.8e-02 |

Three were already in the 60-exemplar plates and sk16/ty12 is the worked
adjacency example — the rarefaction-robust set and the significance-surviving set
converge on the same topologies by independent routes.

## Framing for the deck (per RS)

Call them **"cells lit only by predicted structures (depth-controlled)"**, never
"novel topology." The one-sentence honest read: *of 619 cells lit only by predicted
structures, 90 survive downsampling AFDB to experimental depth — a real but small
predicted-enriched set, and at matched depth experimental space actually carries
MORE unique cells than predicted.*

## Selected examples (13 robust + 2 fragile contrast)

`afdb_selection.pkl`. Chosen greedily to span dominant-H-group, typing (0–4
strands), and skeleton. Diversity available: the 90 robust cells cover 56 dominant
H-groups / 71 skeletons / 27 typings.

Robust (sk,ty,typing,AFDB-folds,dominant fold): 16/12 HHEEH 17 uncharacterized
BT_1490; 91/28 HHEEE 15 UBC-like; 164/3 EEHHH 13 ARM repeat; 171/7 EEEHH 11
Hedgehog/DD-peptidase; 16/19 EEHHE 10 His-Me finger endonuclease; 57/24 HHHEE 10
Ankyrin; 158/27 EEHEE 10 RIFT; 134/0 HHHHH 9 HTH; 180/16 HHHHE 7 SET; 103/16 HHHHE
6 UDP-glycosyltransferase; 96/30 HEEEE 8 Nat acyltransferase; 44/8 HHHEH 6 DENND;
101/30 HEEEE 2 LRR. Fragile contrast: 128/3 EEHHH 2 His-Me finger (SAME fold as
robust 16/19 — vivid depth contrast); 34/10 HEHEH 1 α/β-hydrolase.

## Adjacent-motif context — the antidote to the "novel" drift

`afdb_context.txt`. For each selected cell, its 5 Hamming-1 typing neighbours (flip
one SSE H↔E), each labelled PDB-lit / AFDB-only / dark / impossible + dominant fold.

**13 of 15 selected AFDB-only cells sit one H↔E flip from a PDB-lit named fold**
(mean 2.6 PDB-lit neighbours of ≤5 possible flips). E.g. 16/12 (uncharacterised, 17
folds) is one flip from Restriction-endonuclease (PDB 5f), Metalloprotease/zincin
(7f), and Phosphotransferase (3f). So an AFDB-only cell is a **typing variant
adjacent to well-understood experimental folds**, not novel topology. The lone
exception (134/0 HHHHH) is all-helix, whose neighbours are dark/AFDB-only.

## Renders — DONE (`s5_grid/renders/`)

prosmos_inspect-style, one card per exemplar: **ProSMoS 5×5 query matrix** (code
colours from `interactions.ts`) | **domain by PALSSE element** (helix red / strand
yellow / loop grey) | **domain by matched-motif slot** (SLOT_COLORS blue→red N→C),
caption with the adjacent-motif context. Plus `contact_sheet.png` (all 15, blue
border = robust, red = depth-artifact).

Structure sourcing (answers "how did we have the AFDB structures"): ProSMoS never
touched atomic coords — PALSSE reduced each DPAM domain to SSEs (endpoints +
contact matrix = the metamatrix). For rendering, the **full AF models were fetched
fresh from EBI** (`AF-<acc>-F1-model_v6.pdb`, keyed by UniProt) into
`work/prosmos_2026/afdb_struct/` and chopped to the DPAM range
(`afdb_200m.ecod_domain_range`). Domains are DPAM (`source='dpam'`), NOT TED — the
`ted_domain_*` tables in that schema are a separate annotation layer.

**Caveats:** (1) search ran on AF **v4**, renders show **v6** — same confound
prosmos_inspect discloses ("coordinate-only, small in aggregate"). (2) One exemplar
(`A0A7K4A903`, HTH sk134/ty0) was retired from AFDB; swapped to `A0A286DT57_D5`
(same cell). (3) `render_pymol.py` needs `PYMOL_PATH=~/.pymol` or images watermark.

Scripts: `render_prep.py` (data) → `render_pymol.py` (PyMOL, per key) →
`render_compose.py` (cards + contact sheet). Data in `render_data.json`.
