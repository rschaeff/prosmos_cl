# Single-violation S5 unobserved SSPs

The 14 S5 super-secondary structure patterns (SSPs) that are absent from the
PDB but contain only one "unfavorable" folding feature — flagged in Chitturi
et al. 2016 as candidate de novo protein-design targets.

Source: http://prodata.swmed.edu/ssps/S5/Resort_unob.html
Paper:  Chitturi B, Shi S, Kinch LN, Grishin NV. *Compact Structure Patterns
        in Proteins.* J Mol Biol 428(21):4392–4412 (2016).
        DOI 10.1016/j.jmb.2016.07.022, PMID 27498165.

These are *not* SSPs found in the example/query.txt — they are diagrams of
patterns nature has not produced (as of PDB Dec 2011, 75,196 structures, the
snapshot the paper analyzed). To run ProSMoS against any of them you would
need to translate each diagram into a ProSMoS query matrix; the BMPs alone
are not searchable input.

Files are named `5-<skeleton_id>-<X>-<Y>.bmp` from the upstream website. The
extension is misleading: the files are actually PNG.

## The 14 design targets

| Panel | File | Group | Single feature |
|---:|---|---|---|
| 1 | 5-269-0-0.bmp | (i) five-stranded β-sheet | jump |
| 2 | 5-311-0-0.bmp | (i) | jump |
| 3 | 5-289-0-0.bmp | (i) | jump |
| 4 | 5-288-0-0.bmp | (i) | jump |
| 5 | 5-280-0-0.bmp | (i) | jump |
| 6 | 5-282-0-0.bmp | (i) | jump |
| 7 | 5-306-0-0.bmp | (i) | jump |
| 8 | 5-309-0-0.bmp | (i) | jump |
| 36 | 5-283-1-2.bmp | (iii) 4-strand β-sheet + 1 α-helix | jump |
| 40 | 5-307-1-2.bmp | (iii) | jump |
| 41 | 5-243-1-2.bmp | (iii) | jump |
| 56 | 5-265-7-7.bmp | (vi) 2-strand β-sheet + 3 α-helices | crossing loop |
| 57 | 5-234-7-7.bmp | (vi) | crossing loop |
| 58 | 5-141-7-7.bmp | (vi) | crossing loop |

Panel numbers match Fig. 6 of the paper and the order on the source page.

## Group key (Fig. 6, Chitturi 2016)

- (i)   panels 1–8:   five-stranded β-sheets containing psi-loops
- (ii)  panels 9–34:  β-sandwich (3+2 strands across two sheets) — all have ≥2 features, none single-violation
- (iii) panels 35–43: four-stranded β-sheet flanked by one α-helix
- (iv)  panels 44–47: three-stranded β-sheet with α-helices on both sides — all ≥2 features
- (v)   panels 48–55: three-stranded β-sheet with two α-helices on one side — all ≥2 features
- (vi)  panels 56–63: two-stranded β-sheet with three α-helices on one side

## Feature definitions (Chitturi 2016, Discussion)

- **jump** — a β-strand pair adjacent in the β-sheet but non-adjacent in sequence (split connection)
- **crossing loop** — two inter-SSE connections cross in the plane of the SSP
- **left-handed βxβ** — left-handed chirality of a β-strand–X–β-strand triple (X = α-helix, β-strand, or longer connection); the right-handed form is strongly preferred in nature
