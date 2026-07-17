# Does predicted Ig space contain topology the PDB lacks? — No.

**Question (RS):** build a covering set of templates from PDB Ig, apply it to AFDB
Ig — does it still leave AFDB Ig domains dark? I.e. is the Ig darkness evidence of
novel *predicted* topology, or only of a limited *vocabulary*?

**Answer: vocabulary. 0.01% of AFDB Ig is topologically dark to the PDB.**

Scripts: `s5_grid/igcross/subset_sig.py` → `igcross/sigs.pkl`. Data: the AFDB
(`s5_inv/chunks`, 4.92M records) and experimental-PDB (`s5_pdb_inv/chunks`,
496,359 records) metamatrix DBs; Ig = ECOD X-group 11.*.

## Method — signatures, not fitted templates

The template-fitting design was abandoned; three defects made it unable to answer
the question (see "Why not templates" below). Instead, since a ProSMoS motif
matches **any 5 SSEs in N→C order**, the natural unit is the 5-SSE *subset*:

For every Ig domain in each DB, enumerate every **connected, two-sheet 5-SSE
subset** and reduce it to a topology **signature** = (canonical sheet-partition,
10-bit adjacency). A domain is **dark to the other DB** iff *none* of its subset
signatures occurs in *any* domain there — exactly the condition under which no
template derived from that DB could ever hit it.

This has no fitting step, so the circularity ("derived from AFDB ⇒ covers AFDB")
cannot arise; it uses the whole population; and it permits size-matched controls.

## Result

```
AFDB: 131,665 Ig domains >=5 passing SSEs | 122,919 with a two-sheet 5-subset | 25,687 signatures
PDB :  41,587 Ig domains >=5 passing SSEs |  39,669 with a two-sheet 5-subset | 20,052 signatures
```

| | |
|---|---|
| **AFDB Ig dark to the PDB vocabulary** | **14 / 122,919 = 0.01%** |
| PDB Ig dark to the AFDB vocabulary (control) | 1 / 39,669 = 0.00% |

**Not a sample-size artifact.** Rarefying AFDB to PDB's domain count (5 draws)
gives 23,908–24,112 signatures and AFDB-dark of 0.00–0.02% every time.

**The signature-level asymmetry is real but thin:**

| | AFDB-only | shared | PDB-only |
|---|---|---|---|
| signatures | 5,719 (22.3%) | 19,968 | 84 (0.4%) |
| median domains carrying | 5 | 52 | 1 |
| singletons | 16% | 1% | 63% |
| share of all domain-signature observations | **1.01%** | 99% | — |

AFDB-only signatures are **not noise** (median 5 domains, max 137, only 16%
singletons) — they are reproducible variants. But 11.9% of AFDB Ig domains carry
≥1 AFDB-only signature while only **14** carry *nothing else*. Every AFDB Ig domain
showing a novel 5-subset also presents a PDB-known one. The PDB-only set (median 1,
63% singletons) is noise.

Predicted Ig space = experimental Ig space + a 1% decorative fringe. Same shape as
the negspace sweep's "thin fringe, not empty" ([[project_negspace_finding]]).

## What it settles

The **86.7% Ig darkness in the S5 grid is a property of the grid's vocabulary, not
of the structures.** AFDB holds no meaningfully novel Ig topology.

This kills the "the constraint buys specificity, darkness is the bill" story:
- Not protecting against promiscuity — the two-sheet templates are *more*
  discriminating than grid cells (median top-X purity 61.4% vs 39.6%; every one of
  17 is Ig-dominated, 39–89%, vs grid cells whose dominant group wanders across
  X109 repetitive-α-hairpins / X2003 Rossmann / X9 lipocalins / X243 cystatin).
- Not protecting against exotic predicted topology — there isn't any (0.01%).

**The constraint is representational, not statistical.** Ig alone — one X-group —
presents **25,687** distinct two-sheet 5-SSE topologies against a whole-grid
vocabulary of ~5,696 satisfiable queries. Chitturi's SCC-2 hex lattice is *planar*
and cannot express two stacked sheets at all. Ig is dark because the abstraction is
2D, not because of a tradeoff anyone was managing. Same root cause as the 640
physically-impossible degree-3 sheet queries (`s5_grid_defects_2026-07-16.md`): a
flat grid asked to describe 3D packing.

## Why not templates (the abandoned design)

1. **`merged_tmpl.py`'s core rule is `n_pass == 5`** — *exactly* five length-passing
   SSEs. A complete Ig β-sandwich has 7–9 strands, so the rule selects *incomplete*
   domains: it keeps 5.5% of AFDB Ig but only **0.8%** of PDB Ig. A covering set
   cannot be derived from a 0.8% tail, and would not be a covering set of PDB Ig
   but of unusually small PDB Ig.
2. **Size confound.** PDB Ig domains are more complete than AFDB Ig: median n_pass
   **9 vs 7**. Any AFDB-dark residual would be partly "AFDB Ig are smaller." The
   subset method allows conditioning on n_pass; template-fitting does not.
3. **Circularity.** The existing 17 templates were fitted to a 2,500-core sample of
   *AFDB grid-dark* Ig, so their 85.3% AFDB coverage is partly definitional.

## Trap: the chain-A regex

`merged_tmpl.py` (and my first cut of this analysis) uses `([HE])A\s+…`, hard-coding
chain **A**. AFDB metamatrix records are chain-A only, so it is harmless there —
but PDB records carry chains A/B/C/D/X/1/H, and the regex silently drops every
non-A SSE. My first run reported 75.8% of PDB Ig with `n_pass = 0` and a bogus 35×
core-yield gap (true gap ≈7×). Correct form, validated in `export_palsse.py`:

```python
SSEC = re.compile(r'([HE])([A-Za-z0-9])\s*(-?\d+[A-Za-z]?)\s*--\s*(-?\d+[A-Za-z]?)\s+(\d+)\s+-?\d+\.\d')
```

**`merged_tmpl.py` still carries this bug.** Fix before pointing it at any
multi-chain DB.

## Loose ends

- The **14 dark AFDB Ig domains** — worth an eyeball (pLDDT? DPAM mis-parse?), but
  at 0.011% they change nothing.
- The 5,719 AFDB-only signatures are a real 1% fringe. What *are* they — β-bulges,
  extra strand contacts, sheet-assignment edge cases? Possibly the same
  bulge/bifurcation biology as the 7 impossible-but-hitting queries.
