# Are any S5 topologies serviced only by a structural context? No.

**Question (anticipated from the room):** are there cells lit only by cryo-EM
structures — i.e. topologies that only occur in large complexes — or only by
viral proteins? Expected answer "no", but it needs an actual number.

**Answer: no, and emphatically so at the unit that matters.** At FOLD level not a
single cell is significantly enriched for either EM or viral origin, and no cell
is viral-only across ≥3 distinct folds.

Script: `enum/scripts/strat_fold.py` → `s5_grid/strat_fold.json`.

## Result

Over the corrected 100% PDB sweep: 3,180 cells, 295,360 lit domains
(256,420 method-typed, 285,920 superkingdom-typed).

| | domain level | **fold level (T-groups)** |
|---|---|---|
| base rate | EM 9.7% of domains · viral 7.7% | EM 36.9% of T-groups · viral 19.0% |
| cells tested | 1,925 / 1,954 (n≥20 domains) | 1,725 (n≥5 T-groups) |
| **EM-enriched q<0.05** | **334** | **0** |
| **viral-enriched q<0.05** | **328** | **0** |
| EM-only cells | 0 at n≥20 (7 at n≤5) | 120, of which **1** spans ≥3 folds |
| viral-only cells | 60, three with >100 domains | 76, of which **0** span ≥3 folds |

The only cell worth a second look is `sk173/ty0` — EM-only across 3 distinct
T-groups. Everything else labelled "X-only" is a cell occupied by one or two
folds that happen to be X, which is a statement about what has been solved, not
about the topology.

## Why the domain-level version was wrong (worth keeping)

At domain level the viral cut looked spectacular: 60 viral-only cells, the largest
with 168 domains against a 7.7% baseline — p ≈ 1e-150. It is entirely redundancy:

| cell | domains | T-groups | organisms |
|---|---|---|---|
| sk178 ty16 | 168 | **1** (4967.1.1) | Hepatitis C virus |
| sk178 ty2 | 148 | **1** (4967.1.1) | Hepatitis C virus |
| sk130 ty24 | 130 | **1** (2484.1.1) | Enterobacteria phage |
| sk184 ty11 | 41 | **1** (304.48.1) | Poliovirus / Enterovirus |
| sk7 ty16 | 25 | **1** (4967.1.1) | Dengue / JEV |
| sk5 ty3 | 19 | **1** (1042.1.1) | Coronavirus spike |

Every one is a single T-group solved dozens-to-hundreds of times — HCV
protease/polymerase, poliovirus capsid, coronavirus spike: drug targets, not
topological classes. `n = 168` is one independent observation, not 168. Correcting
the unit takes 334 + 328 "significant" cells to **zero**: a ~100% artifact rate,
and the sharpest demonstration in this project of why fold-level normalisation is
not optional.

## Two data traps hit on the way (both silent)

1. **`afdb_200m.taxonomy` contains no Viruses at all** (Eukaryota 897,882 /
   Bacteria 96,647 / Archaea 4,221 / none 440). It is AlphaFold-derived and AFDB
   v4 excluded viral proteomes, so routing PDB domains through it returns
   *0 viral by construction*. The first viral run reported "global viral rate
   0.00%, 0 viral-only cells" — a guaranteed zero with no biology in it.
   **Correct source: `ecod_commons.domain_taxonomy` on dione** (PDB-native, keyed
   by `ecod_uid`, 3.1M rows, Viruses 159,444).
2. **Current UniProt resolves only 34% of our AFDB accessions**, and the survivors
   are 51% eukaryotic against a true 27% — TrEMBL purges since the AFDB v4
   snapshot skew toward environmental/redundant entries. Using it would have
   manufactured the expected eukaryote enrichment. Correct source:
   `afdb_200m.protein` (214.7M rows, `uniprot_acc` → `organism_tax_id`), 99%
   coverage, contemporaneous with the models searched.

## The structural argument that should accompany the number

The S5 motif is five SSEs **within one domain**. An inter-chain interface cannot
create a cell. So a genuinely "complex-only" topology would require the *fold
itself* to occur only in assembly-forming proteins — a far stronger claim than
"EM solves complexes", and one this data gives no support for.

Graded method preference does exist and is unsurprising (at domain level, large
all-helix cells run 14–17% EM vs a 9.7% baseline; helical assemblies are EM's
home turf) — but preference is not exclusivity, and at fold level even the
preference vanishes.
