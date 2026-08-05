# Two defects in the fold-level S5 grid (found 2026-07-16)

Both were found by eye, from structure in the AFDB fold-level heatmap that had no
business being there: bars of depleted cells at the HHHHH edge covering several
adjacent typings. Neither is a searchmatrix bug; both are downstream of the sweep.

Corrected figure: `work/prosmos_2026/s5_grid/afdb_vs_pdb_fold_level_final.png`
(script `plot_nr_final.py`). Rebuilt array: `grid_afdb_nT_rebuilt.npy`.
Impossible-query mask: `impossible_mask.npy`.

---

## Defect 1 — 670 AFDB cells silently zeroed (a sharding race, aliased onto the all-helix columns)

**What happened.** `grid_afdb_nT.npy` was assembled from the 40 output files of
`nt_shard.py` while **the last 7 were still buffered**. `nt_shard.py` opens its
output with `open(...,"w")`, so the file exists immediately but stays *empty*
until `close()` flushes it. The assembler globbed `ntshards/*.tsv`, read seven
empty files, found no rows, and left those cells at **0**. No error, no warning.

Evidence: the 7 affected shards (24, 32, 34, 36, 37, 38, 39) are exactly the last
7 to finish (mtimes 18:29:17–18:30:26); `plot_nr.py` ran at 18:30:35 — **9 seconds**
after shard 32 closed. In all 7, *every* row is zero in the old array. The error is
strictly one-directional: 670 cells `old=0, new>0`, and **zero** cells the other way.

**Why it landed on the HHHHH edge.** This is the part worth remembering:

```
nt_shard.py:  qd = sorted(os.listdir(H))[shard::40]     # stride 40
grid row length = 32 typings;  index = sk*32 + ty;  shard = index % 40
gcd(32, 40) = 8   =>   shard ≡ ty (mod 8)
```

The sharding stride **aliases against the typing axis**. Each shard only ever
receives typings `ty ≡ shard (mod 8)` — shards 0/8/16/24/32 hold *only*
ty ∈ {0, 8, 16, 24}, the all-helix columns. Those are also the biggest cells
(sk=132 ty=0 alone has 421,065 hit files), hence the slowest to `scandir`, hence
the last to finish, hence the ones lost. Losing one worker erased a whole typing
column instead of scattering invisibly. With a stride of 41 (coprime to 32) the
same race would have produced a diffuse ~17% undercount that nobody would have
spotted.

**Damage.** 53,397 folds destroyed = **23% of the AFDB fold total**
(181,845 → 235,242, +29%). Concentrated at ty=0 (75 cells), ty=8 (63), ty=16 (58),
ty=24 (58).

**The shard TSVs were correct all along** — `s5-0132-0000` records
`421065 records / 1534 folds`. Only the assembly lost it, so the rebuild needed no
recomputation. The assembler itself was an unsaved inline command, which is why
there is no script to fix; `plot_nr_final.py` now rebuilds from the TSVs directly.

**PDB is unaffected.** `pdb_nt.py` is single-process with no sharding. Its 10
`structures>0, folds=0` cells are tiny (1–10 structures) and are genuine
unclassified domains. So the defect biased **AFDB only, downward**, while PDB was
correct — every AFDB-vs-PDB fold-level number was skewed the same way.

**Lessons.**
- A worker that writes `0` on failure is indistinguishable from a real zero.
- Never read a producer's output files without a completion barrier (a `.done`
  sentinel, or an explicit expected-row-count assert: 6,336 rows, not "whatever
  the glob returned").
- **Choose shard strides coprime to any axis of the output.** Aliasing turns a
  random loss into a structured one that reads as signal. Here it was a gift —
  the bars are what made it visible.

## Defect 2 — 640 of 6,336 queries are physically unsatisfiable

**What.** A query that declares `sheetS` over n strands and demands a contact
graph in which some strand has **≥3 lateral neighbours**. A β-sheet is
one-dimensional: strands sit side by side, so a strand has **at most 2** lateral
neighbours, and n strands in one sheet have exactly n−1 contacts.

The skeletons come from a **hex grid**, where a node may have up to 6 neighbours.
The enumeration types such a skeleton `EEEEE`, declares all five strands one
sheet, and never checks the result is realizable.

**Scale.** 640/6,336 (10%) grid-wide; **98.9% of them are dead** (633/640).
Sharpest in the 5-strand-sheet slice:

| max lateral neighbours demanded | lit | dead |
|---|---|---|
| 2 (possible) | 6 | 2 |
| 3 | 0 | 68 |
| 4 | 0 | 64 |

**Verified directly.** All 134 dead 5-strand-sheet queries were run against a
mini-DB of **8,155 real Rossmann domains** (ECOD T-group `2003.1.1`, extracted
from the PDB metamatrix): **0 hits**, rc=0 — with a positive control (the 6 lit
5-strand queries) hitting on the same harness. Build:
`rossmann.db`, manifests `five_dead.man` / `five_lit.man`.

**Rossmann is NOT excluded by the grid.** `sk=126` carries exactly the Rossmann
contact graph (the path 3-2-1-4-5, edges {(1,2),(2,3),(1,4),(4,5)}) and is lit.

**The rule is not exact:** 7 impossible queries do hit. Real sheets are not
perfectly 1-D — β-bulges and bifurcated sheets can give a strand a third
neighbour. So 640 is an **upper bound**, and "essentially unsatisfiable" is the
honest phrase.

**Two hypotheses the data refuted** (recorded so they are not re-proposed):
1. *"The generator forces `-` on non-adjacent strands it declared same-sheet,
   making the query self-contradictory."* — **No.** In real records, 84% of
   same-sheet non-adjacent pairs are coded `-` (4805/5700). The queries match
   reality. Deadness vs n_dash is non-monotonic (48/53/98/58/100/100/25%).
2. *"The generator only emits `t` (antiparallel), so parallel sheets can never
   match."* — **No.** Sheet-internal codes are `t` 5,168 / `c` 3,444 / `-` 2,919.
   Both orientations are emitted. (Real Rossmann sheets *do* code adjacent pairs
   `c` 1431 vs `t` 163 — parallel — but `sk=98` and others demand `c`.)

---

## Corrected numbers

| | before | after |
|---|---|---|
| AFDB folds total | 181,845 | **235,242** (+29%) |
| dead cells | 3,375 | **2,705** (670 were Defect 1) |
| …physically impossible | — | 633 (23% of dead) |
| …genuinely dark | — | **2,072** |
| grid occupancy | 47% (of 6,336) | **64%** (of 5,696 satisfiable) |
| median log₂(PDB/AFDB), all cells | +0.000 (1.00×) | −0.585 |
| …data-bearing, satisfiable | — | **−1.222 (0.43×)** |

**A third, separate defect in the old figure:** the log₂ panel let cells dead in
*both* datasets evaluate to `log2(NA/NP) = +0.0004` and render as neutral grey —
**pixel-identical to "the two agree perfectly."** At 42% of cells they were the
largest block in the population and **pinned the all-cell median at exactly 0.00**.
The reported *"1.00×, experimental and predicted structure space contain S5 motifs
at indistinguishable fold-level rates"* was measuring the pseudo-count, not
biology. **Retracted.**

**Do not replace it with 0.43×.** The PDB sweep is 68% complete and partial
coverage depresses PDB counts in exactly that direction. The comparison is not
interpretable until `624866` finishes. This is an independent reason to hold the
prosmos_inspect re-export.

## What survives

The **dark-side story stands**. Only 633/2,705 (23%) of the darkness is
definitional, leaving 2,072 genuinely dark cells. Neither defect touches the
**domain-level** dark fraction: an impossible query would never have hit anything,
and Defect 1 only ever destroyed fold *counts* in cells that were lit anyway. Both
change the cell denominator and the heatmap's interpretation, not which domains
are lit.

**Still to do:** propagate the corrected numbers into the deck and
`session_2026-07-16_handoff.md` (both still carry 1.00× and the 6,336
denominator); decide whether the paper's S5 vocabulary was ever meant to include
the 640 (if not, the enumeration count is ~5,696, which touches the enum work);
and look at the 7 impossible-but-hitting queries — if they are bulged/bifurcated
sheets, that is real biology the hex abstraction catches by accident.
