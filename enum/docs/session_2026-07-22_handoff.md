# Session handoff 2026-07-21/23 — structure-gap closure, searchmatrix bug, targeted extension, query-set correction

**Full artifact inventory: `~/work/prosmos_2026/INVENTORY.md`** (audited from disk).
This handoff is the narrative; the inventory is the ground truth.

## Arc of the session

1. Built the `/ruczinski` browse section in prosmos_inspect (list + detail + Mol*).
2. Found and fixed three data bugs in it (missing AFDB structures, ecod_ routing,
   segment chain).
3. Discovered the AFDB sweep is contaminated with experimental ECOD reference
   domains (7,735 `ecod_*` records, all `source_type='pdb'`).
4. RS generated the missing experimental structures → coverage 33.5% → **93.8%**
   (500,313 → 1,400,111 domains). Rebuilt the experimental DB with TRUE provenance.
5. Re-ran the Ruczinski census on 1.4M: **27 exp / 47 AFDB / 20 predicted-only /
   1 absent**. Only ONE motif crossed from zero → the predicted-only set is not a
   sampling artifact. Panel 28 still absent from both.
6. **Depth analysis** (the key conceptual result): the raw exp/predicted split is
   half selection-effect (the 48 were *chosen* for experimental absence in 2002)
   and half missing depth correction. Depth-corrected, only 10/20 predicted-only
   survive BH; pooled observed/expected ratio is 0.88 → no global depletion.
7. Built two targeted-extension enumerators (S6 near an S5 seed; S4→S5 in the
   Ruczinski encoding) — the route to S6 space near specific motifs.
8. **Found a searchmatrix parsing bug** that had corrupted two full sweeps.
9. **Found the production binary was stale AND unoptimised** (14× slow).
10. **Caught my own error**: re-swept the corrected DB with the wrong query set
    (adjonly, not the canonical typed graph198). Re-sweeping correctly now.

## The searchmatrix bug (commit 1e5546c) — READ THIS

`intMnumofele()` bounded the SSE loop by LINE LENGTH (`while beginposition<size`),
never consulting `numberElment`. Records whose coordinates have five integer
digits (large assemblies, e.g. `46198.078`) overflow their 8-char field, making
the line longer than `36 + 68*numberElment`. The loop parsed the surplus as
PHANTOM SSEs from misaligned bytes; where it ran past the end, an uncaught
`std::out_of_range` aborted the chunk. 75 of 1.4M records overflow (0.005%); 4
crashed a chunk apiece (~11,200 domains), the rest injected phantom SSEs.

`array_inverted.sbatch` never propagated the exit code, so SLURM logged
`COMPLETED 0:0` on the crashes — invisible through two sweeps. Now `exit "$RC"`.

The old mixed-PDB sweep had the same 4 crashes, so the **published 26/47/21/1
census was computed with chunks missing**. AFDB is unaffected (origin-centred
models never reach 5-digit coords). Verified: on chunk_0321, unpatched exits 134
after 1,760 hits; patched exits 0 with 5,063.

## The stale/unoptimised binary (commit fbd3d0f)

`build/searchmatrix` (the path sweeps invoke) was never refreshed after commit
8a853a7's 2.9× hot-loop optimisation — it stayed at the 07-11 build. Compounded
by the documented build command missing `-O` (a further 3.8×). Net ~14× slow,
which is why S5 chunks took 7h instead of 15min. Added `searchMatrix/Makefile`
(depends on `src/*.h`, pins `-O2`) and a `make stale` guard wired into
`submit_inverted.sh` that refuses to dispatch against a stale binary.

## The query-set error (commit d795100) — the reason for the census

After the binary fix I re-swept the 1.4M DB with **`queries_adjonly`**
(orientation-blind) instead of the canonical **`queries_graph198_alltypings`**
(typed). adjonly matches ~20× more, inflating sk016 ty12 from 0 to 229,424
domains and giving a depth ratio of 427 (should be <1; all depth-q came back
`nan`). The two share `s5-XXXX-XXXX` hit naming, which hid the swap.

Fixed: recompute defaults to `s5_exp_full_g198`, uses DISTINCT lit domains for Q
(not `Pd.sum()`, which double-counts), and asserts `0<Q<1` so the wrong query
set fails immediately. The adjonly sweep is quarantined
(`s5_exp_full_inv/WRONG_QUERY_SET.md`); the invalid grids are renamed
`*.INVALID_adjonly`.

**Currently re-sweeping** `s5_exp_full_g198` (job 636570, graph198, fixed binary,
64 concurrency). This is the sole blocker for a valid significance recompute.

## Targeted extension (commit 416589f)

Two enumerators in `enum/scripts/`:
- `generate_s6_targeted.py` — extends an S5 grid cell by one SSE via the lattice
  combine. **CAVEAT: defaults to `--matrix adjonly`**, so it extended the wrong
  baseline; must be redone against graph198.
- `extend_motif_query.py` — extends a typed motif (Ruczinski encoding) by one
  SSE; used for panels 3/10/23 (job 633378, DONE). This one is correct.

Verified insight: the ~15–60 extensions of one seed form a SIBLING SET (same
parent, same depth), so read occupancy relative to siblings — absolute S6
emptiness is uninformative (constraints multiply, some cell is always empty).

## The contamination / depth-significance thread (NOT yet recomputed correctly)

The published 5 depth-significant AFDB-only cells include TWO
(`sk016 ty12`, `sk033 ty07`) whose only "AFDB-only" support came from
experimental ECOD domains that had no structure file at the time — the
structure-gap generating false positives, not just capping coverage. Both now
have searchable structures. The recompute (pending the g198 sweep) will show
whether they survive. Expectation: they resolve to PDB-lit.

## Open, in dependency order

1. g198 sweep 636570 finishes → run `recompute_significance.py`.
2. Cluster-unit depth correction (per-structure Poisson over-calls given AFDB
   redundancy; matched 50%-cluster unit is the honest test).
3. Re-derive "panels 3/10 have no S5 parent" against graph198 — conclusion likely
   holds, but the stated reason ("grid is orientation-blind") is FALSE for the
   typed grid.
4. Regenerate `generate_s6_targeted` output against graph198; re-run s6_targeted.
5. Fold corrected numbers into `/ruczinski` and the significant-browse app.
6. Push prosmos_cl (main still has the `PGPASSWORD='***REMOVED***'` placeholder;
   the branch fixed it in f141456).

## Commits this session (prosmos_cl)

- `c67c0af` Ruczinski 48-panel transcription + corrected census
- `f141456` read ECOD password from ~/.pgpass
- `416589f` two targeted-extension enumerators
- `1e5546c` searchmatrix numberElment bound + exit-code propagation
- `fbd3d0f` searchMatrix Makefile + stale-binary guard
- `d795100` canonical significance recompute + query-set/depth-ratio guards

(prosmos_inspect: /ruczinski section + the three data-bug fixes + deep links,
committed separately in that repo.)
