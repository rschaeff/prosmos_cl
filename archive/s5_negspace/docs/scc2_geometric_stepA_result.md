# Step A result — geometric SCC-2, S5 (paper-faithful count `N'`)

Executes Step A of `scc2_geometric_plan.md`: Phase 0 (encode Fig-S3 grids
geometrically) + Phase 1 (geometric SCC-2 behind a flag) + Phase 2
(re-enumerate → `N'`, diff vs 198, clean-mapping proof). **No analysis re-run.**

## Headline

| SCC-2 backend | S5 skeleton count |
|---|---|
| graph (unlabeled adjacency, current default) | **198** |
| geometric (Fig-S3 hex-congruence, paper-faithful) | **140** = `N'` |

The paper-faithful S5 skeleton count is **140**. All 66 existing tests still
pass with the default (graph) backend unchanged.

## The fix is sound (Phase 2 diff)

`enumerate_skeletons(5)` under geometric SCC-2 is a **strict subset** of the
graph result: geom(140) ⊂ graph(198), symmetric difference = **58 dropped, 0
added**. The 58 dropped are *entirely* accounted for:

- **all 58** live in grid-e's 5-edge adjacency graph
  (`((0,1),(0,2),(0,3),(1,2),(1,4))`);
- **none** is hex-congruent to any Fig-S3 grid.

They are exactly the **bent 1-2-3-4 path** variants: graph-isomorphic to the
straight Fig-S3(e) grid but not the Fig-S3 *shape*. The graph test admitted
them; the geometric test rejects them. This is the root cause named in the plan.

Per-grid breakdown of the 140 survivors (every survivor maps to a Fig-S3 grid;
`None` count = 0):

| Fig-S3 grid | edges | survivors |
|---|---|---|
| (d) P5 path              | 4 | 12 |
| (e) collinear + corner   | 5 | 40 |
| (g/h) star               | 6 | 28 |
| (f) dense                | 7 | 60 |
| **total** | | **140** |

The graph 5-edge class had 98 members (40 congruent-e + 58 bent); geometric
SCC-2 keeps the 40, drops the 58.

## Clean-mapping proof (the 100% vs ~80% claim)

The prior blocker: under the graph whitelist, ~20-30% of skeletons in the
"grid e" class had a **paper handedness triple collapse to coplanar** (sign 0)
in our 3D realization, so the paper's per-grid handedness set could not be
applied — a *dirty* mapping.

With geometric SCC-2:

- **grid e (the only split class): 40/40 survivors** have *both* paper
  handedness triples `{(1,4,5),(2,3,5)}` (reference labeling) non-coplanar,
  under the congruence-induced labeling. **100% clean.**
- grids d / g-h / f: **0 misses** — every skeleton in those classes is
  hex-congruent to its grid (1:1), so their paper handedness applies uniformly
  by construction.
- The 58 dropped bent-e skeletons carry the coplanar collapses (e.g. 1/10,
  1/10, 4/10 coplanar triples in the first three) — the concrete reason they
  have no paper grid and no paper handedness.

Old clean fraction over the 198 (grids d+gh+f fully congruent, e-class 40/98) =
140/198 ≈ **70.7%**; new = **140/140 = 100%**. The "dirt" was precisely the 58.

## What landed (code)

- `enum/src/ssp_enum/geom_scc2.py` — Fig-S3 S5 grid point-sets, hex-isometry
  canonical form `geometric_canonical`, `grid_of`, `passes_geom_scc_2`.
- `enum/src/ssp_enum/compactness.py` — `_SCC2_MODE` toggle + `set_scc2_mode`;
  `passes_scc_2` dispatches graph|geometric. **Default stays `"graph"`** so the
  current 198-based darkness/negspace analysis is untouched.
- `enum/tests/test_geom_scc2.py` — locks 198→140, the per-grid counts, the
  strict-subset + 58-drop diagnosis, S3/S4 invariance, isometry-invariance.
- Repro: `work/prosmos_2026/geom_scc2_stepA.py`.

## Reproduce

```python
from ssp_enum.compactness import set_scc2_mode
from ssp_enum.enumerate import enumerate_skeletons
set_scc2_mode("graph");     len(enumerate_skeletons(5))   # 198
set_scc2_mode("geometric"); len(enumerate_skeletons(5))   # 140
```

## What Step A does NOT do (deferred)

- **Phase 3** — apply minimal handedness (incl. type/sheet-conditional f-triples)
  + minimal interactions (solid=required, broken=≥1-mandatory) behind
  `paper_faithful=True`. Needs f's 7-list and gh's 6-list read off Fig S3.
- **Phase 4** — regenerate the paper-faithful query set on the 140 basis + the
  cheap impact probe (does the darkness/negspace conclusion survive the basis
  change?). **This is the go/no-go for a full analysis re-run.**
- S3/S4 geometric whitelists (panels a-c of Fig S3) — S5 was the analysis target.

## Bottom line

`N' = 140`. The geometric fix is sound: the 58 dropped are exactly the
non-Fig-S3 bent variants, and every one of the 140 survivors maps cleanly to a
Fig-S3 grid with non-coplanar paper handedness. The 198→140 change is opt-in and
does not touch the current analysis until Phase 4's probe says it should.
