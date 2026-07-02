# S5 matrix — temporal & ECOD-promiscuity analysis

Analyses layered on the S5 skeleton×typing hit matrix (198 skeletons × 32 H/E
typings = 6,336 cells) searched against the **ECOD manual-representative DB**
(the 19,085-domain experimental set behind `figures/s5_hit_grid.png`).

Two questions:
1. **Did the PDB *change* the sampled matrix over time, or just *fill it in*?**
2. **Is each cell unitypical (one ECOD group) or promiscuous (many)?**

## Data sources
- **Hits**: `~/work/prosmos_2026/ecod_search_v4/hits/s5-SK-TY/pdbNNNNNNNNN.txt`
  (+ `ecod_search_v4_retry/`). Each file = one matching rep domain; the name is
  the zero-padded **ecod_uid**; the first `MOTIF` block gives the 5 matched SSE
  residue ranges (PDB author numbering). 196,928 s5 hits.
- **Classification / PDB code**: `ecod_rep.domain` on **dione:45000** (db
  `ecod_protein`). `t_id` = `X.H.T` → H-group = first two components, T-group =
  full `t_id`. PDB code parsed from `ecod_domain_id` (`e2gp4A1` → `2gp4`);
  domain chain + range from `chain_id` / `pdb_range`.
- **Deposition dates**: `/usr2/pdb/derived_data/index/entries.idx`
  (IDCODE, ACCESSION DATE MM/DD/YY). 98.9% of hits dateable.
- **Structures** (for exemplar renders): `/usr2/pdb` mmCIF mirror.

## Pipeline (scripts in `../scripts/`)
Each script has an `SP` constant pointing at a scratch dir holding the
(regenerable) intermediates: `manual_hits.txt` (`find hits -name '*.txt'`),
`ecod_uid_class.tsv` / `ecodrep_uid_map.tsv` (dumped from `ecod_rep.domain`),
`pdb2date.tsv` (from `entries.idx`). Set `SP` before running.

| step | script | output |
|------|--------|--------|
| time cube + viz | `build_s5_timeseries.py` | `figures/s5_timeseries.{json,html}` |
| time static figs | `plot_s5_timeseries_static.py` | `s5_timeseries_fill_vs_change.png`, `s5_timeseries_frames.png` |
| promiscuity | `build_s5_promiscuity.py` | `figures/s5_promiscuity.json` |
| promiscuity heatmap | `plot_s5_promiscuity.py` | `s5_promiscuity.png` |
| promiscuity explorer | `build_s5_promiscuity_explorer.py` | `figures/s5_promiscuity_explorer.html` |
| curated exemplar select | `build_s5_curated_exemplars.py` | `s5_exemplar_render_spec.json` |
| PyMOL render | `render_s5_exemplars.py` (`pymol -cq`) | `exemplar_renders/*.png` |
| montage assembly | `assemble_s5_montage.py` | `s5_exemplars_unitypical.png`, `s5_montage_sk*.png` |

> PyMOL needs `PYMOL_PATH=~/.pymol` (license) and `OMP_PROC_BIND` unset.

## Findings
**Time (`s5_timeseries*`).** Distinct occupied cells reach 90% of their final
extent by ~2009 and plateau at ~2,326, while cumulative hits keep climbing
(107k in 2005 → 195k in 2019). New cells opened per 1,000 new hits crashes from
~1,000 (1980s) to near-zero by ~2000. → the reachable region was **mapped by the
mid-2000s; later structures densify it** (filled in, didn't keep changing).
Note: this manual-rep DB is a ~2019 ECOD snapshot, so the axis ends at 2019
(saturation happens well inside the window).

**Promiscuity (`s5_promiscuity*`).** 27% of occupied cells are unitypical
(nT==1 — the geometry uniquely tags one topology); the rest span up to 721
T-groups (median 4). Promiscuity tracks commonness: simple/all-α arrangements
(e.g. `HHHHH` 5-helix bundles) recur across hundreds of folds; rare geometries
are fold-specific. `s5_promiscuity_explorer.html` gives click-to-inspect
group breakdowns with exemplar domains linked to ECOD/RCSB.

**Curated exemplars.** `s5_exemplars_unitypical.png` (4 fold-specific
geometries) and `s5_montage_sk{0132_ty00,0087_ty24,0126_ty28}.png` (one geometry,
6 folds each; three distinct skeletons/typings). Each montage carries a banner
with the skeleton schematic (hex-lattice nodes, rainbow N→C, ○ helix / □ strand)
and the ProSMoS query matrix (interaction codes, via `plot_skeleton_schematic`) — cartoon with the 5 matched SSEs colored rainbow N→C, embedded in the
grey domain.
