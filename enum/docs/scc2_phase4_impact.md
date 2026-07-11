# Phase 4 — paper-faithful query generation + impact probe

Generates the paper-faithful S5 query set and measures how the basis change
(198→140) + the paper's three rules move the darkness result vs the old basis.
This is the plan's **go/no-go** for a full darkness/negspace re-run.

## Query set generated

708 paper-faithful `.query` files at
`work/prosmos_2026/queries_paper_faithful_140/s5/` (naming
`s5-<skel>-<typing>-<variant>`). Grid-e SSPs carry two variants (the "≥1
mandatory" disjunction). Companion sets for the decomposition (fresh, same code
path): `queries_geom140_alltypings` (140 skels, rules off),
`queries_graph198_alltypings` (198 skels, rules off).

## Harness

`work/prosmos_2026/run_probe.sh <qdir> <db> <out>` — every query vs a
metamatricesDB, isolated cwd per worker (searchmatrix writes a DB-named log +
`../sheetbug` into cwd, so parallel workers must not share it; **query paths must
be absolute** — a relative path breaks after the worker `cd`s). Positive control:
1tim record → 38 hit pairs / 16 skeletons (matches the historical result).

## Decisive per-structure result (1tim)

| query set (all from the same code path) | 1tim S5 hits |
|---|---|
| graph-198 all-typings (old basis) | 38 (16 skeletons) |
| **geom-140 all-typings** (basis change only) | **38 (16 skeletons)** |
| geom-140 paper-faithful (all rules) | **0** |

- **Geometric SCC-2 is hit-neutral on 1tim**: all 16 matching skeletons survive
  (none was in the dropped-58 bent variants). The 198→140 change does not lose
   1tim's motifs.
- **The typing rule (Rule 1) alone removes all 38 hits.** 1tim matches grid-e
  queries typed with a helix on the 1-2-3-4 collinear run (e.g. `HHEHE`), which
  Rule 1 ("four or more collinear points always form a beta sheet", Chitturi 2016
  SI §1) forbids. So a canonical TIM barrel goes **dark at S5** under paper-
  faithful — some paper-faithful darkness is a Rule-1 modeling artifact, not
  geometric absence.

## Population probe (2050 random AFDB reps, fresh same-code-path)

| query set | structures hit | hit % |
|---|---|---|
| graph-198 all-typings (old basis) | 478 | 23.3 |
| **geom-140 all-typings** (basis change only) | **478** | **23.3** |
| geom-140 paper-faithful (all 3 rules) | 241 | 11.8 |

**Basis change 198→140 is exactly hit-neutral:** the two hit *sets* are
identical (0 lost, 0 gained, 478 both). The dropped-58 bent variants contribute
**zero** hits on the sample — geometric SCC-2 removes only redundant geometry, no
matches.

**Rules on the identical 140 set** (all-typings → paper-faithful): 234 both,
**244 lost to the typing rule**, **7 recovered by minimal handedness/
interactions**. The paper rules **halve** the hits (478 → 241, −50%), and the
effect is **almost entirely the typing rule** (Rule 1); the minimality relaxation
(Rules 2/3) recovers only a handful.

### Caveat — these are RELATIVE, not population, rates

The fresh `_build_record` queries hit 23% of this stride-sample, whereas the full
S5 search hit ~0.55% of AFDB, and the qian-package 198 queries hit 4.5% of *this*
sample. So (a) the stride sample is enriched for hitters vs a uniform draw, and
(b) the fresh queries are more permissive than the qian set (478 vs 93 on the
same sample — an unreconciled generation difference). **Absolute dark% here is
not the population darkness.** Only the *relative* comparisons above — all within
the fresh same-code-path family — are valid; they isolate the basis-change and
rule effects cleanly.

## Verdict (go/no-go)

- **Adopt geometric SCC-2 (140) as the skeleton basis — GO, low risk.** It is
  provably **hit-neutral** (identical hit sets to the 198 on both 1tim and the
  population sample), so every darkness/negspace/D1-D2 result computed on the 198
  **transfers unchanged** to the 140. The basis change alone does *not* require a
  re-run. This is the clean win: the paper-faithful S5 count is 140 and the
  analysis stands.
- **The full paper-faithful query set is a separate, larger intervention.** It
  roughly **halves** hits, driven almost entirely by the **typing rule** — which
  makes even a canonical TIM barrel (1tim) dark at S5. Adopting it would
  materially change *which* structures hit, so it should **not** be a blind
  re-run. First decide, as a scientific/framing question, whether Rule 1's
  "collinear-4 ⇒ β-sheet" is the right lens given it darkens real folds.
- **New finding to foreground:** separate **geometric darkness** (dark with rules
  off) from **rule-induced darkness** (dark only because Rule 1 forbids a typing
  the structure actually shows). The 20 dark-gallery structures are *geometrically*
  dark (0 hits even rules-off); 1tim is *rule-induced* dark. That distinction is
  a cleaner, more defensible story than a single "dark fraction."
- **To reconcile before any absolute-number claim:** the fresh-vs-qian query
  permissiveness gap (478 vs 93) — the darkness analysis's real "old basis" is
  whatever query set produced the 0.55%, which must be pinned before quoting a
  paper-faithful population dark fraction.

## Reproduce

```
# generate query sets (enum/):
python -c "... skeletons_to_records(s5, paper_faithful=True) ..."   # 708
# probe:
work/prosmos_2026/run_probe.sh <qdir> <db> <out.txt>
# sample DB: work/prosmos_2026/sample_afdb_2000.db (2050 reps, stride 2400)
```
