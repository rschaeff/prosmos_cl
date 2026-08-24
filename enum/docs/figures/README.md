# Dark-but-compact AFDB domains (Case-2 darkness examples)

20 AlphaFold/AFDB domains that are **compact, ≥5 SSEs, and hit ZERO S5 cells**.
Concrete instances of the "darkness is compact geometry S5 can't template"
finding (dataset_analysis_plan.md → *Darkness is compact geometry S5 can't
template (Case 2, not Case 1)*).

## How they were found
- Source DB: `afdb_db/metamatricesDB.clean` — the 4,921,931 AFDB domains searched
  against the full S5 matrix (`s5_full_afdb/`).
- **Dark** := domain name never appears as a hit file `s5_full_afdb/hits/*/<name>.txt`
  (all hits are `pdb(dpam|ecod)_…` prefixed; 15,327 distinct hitting domains total,
  ≈0.55%). Spot-checked candidates directly against the whole hits tree → 0 files.
- **≥5 SSEs** from the `.ssd` header count.
- **Compact** := **complete (or near-complete) SSE contact graph** in the ProSMoS
  interaction matrix — i.e. ≥5 SSEs each with degree ≥2; the 20 chosen all have
  mean degree = N−1 (**every SSE contacts every other, Kₙ**). This is ProSMoS's
  *own* contact metric — the same min-CA-distance criterion that defines a hit —
  so these are not sparse/extended in the sense that would trivially exclude them.
- Parser: `../find_dark_compact.py` → `../dark_compact_afdb.tsv`
  (**2,459,512** dark compact ≥5-SSE domains — the bulk of AFDB's compact ≥5-SSE
  space is dark). Selection/annotation reproducible from `../dark20_annotated.tsv`.

## The 20 (see `dark_compact_montage.png`)
Span α, β, α/β. 8 carry confident ECOD T-groups — they are **real named folds**,
not junk:
| fold (T-group) | examples | class |
|---|---|---|
| winged helix (101.1.2) | A0A0D1DCD2 | HTH α+β |
| LysM (101.15.1) | A0A2H0H7D5, A9AVL8 | HTH α+β |
| F-box (145.1.1) | A0A177CFU8 | α |
| EGF-like growth-factor (389.6.1) | T1IN09, A0A812UZP3 | small β |
| HET-s prion left-handed β-helix (208.5.1) | A0A239CTP7 | β-solenoid |
The other 12 are compact but carry no confident ECOD assignment.

## Honest caveat on "compact"
"Compact" here = **dense SSE-contact graph**, ProSMoS's hit-defining metric.
The all-β 8-strand examples (`EEEEEEEE`, incl. the HET-s β-helix) are dense in
contact yet **elongated in space** (β-sandwich / single-turn β-solenoid) — they
read as ribbons in the montage, not globular balls. The α and α/β examples
(winged-helix, LysM, 5-helix bundles) are the genuinely globular ones. Either
way the point holds: by ProSMoS's own contact criterion these have a fully
connected 5+ SSE core, yet match **none** of the 198 S5 2D-hex skeletons.

## Files
- `dark_compact_montage.png` — 4×5 labeled gallery
- `models/AF-*.pdb` — AF models (EBI; A0A0X7JK86 & A0A2G7E281 absent at EBI, swapped)
- `renders/*.png`, `render.py`, `render.pml`
- `../dark20_annotated.tsv` — name, unp, range, N, H/E string, T-group, fold, status

---

## v2 (length-filtered — supersedes the above 20)
The S5 queries hard-filter each position by SSE length (**H ≥ 8, E ≥ 5** residues).
The v1 selection ignored this and **8/20 failed** it (dark for a trivial
too-short-SSE reason, esp. the small disulfide-rich EGF cases). Re-derived counting
only length-passing SSEs (`../find_dark_compact_v2.py` → 1,913,192 length-clean
dark compact ≥5-SSE domains). New set: **`dark_compact_montage_v2.png`**,
`../dark20b_annotated.tsv` (per-SSE lengths + assignment), renders in `renders_b/`.
All 20 have ≥5 length-passing SSEs in a complete contact graph, all verified dark.
Assigned folds: winged helix (101.1.2), LysM (101.15.1), HET-s β-helix (208.5.1),
pectin-lyase β-helix (207.2.1). β-helix examples are dense-in-contact yet elongated
(honest caveat retained).
