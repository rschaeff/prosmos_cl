# Interim hand-translated queries for the 14 design targets

**Status:** interim, superseded by the upcoming `enum/` package output.

These 11 `.query` files were hand-translated from the 14 unobserved-S5 BMP
diagrams on http://prodata.swmed.edu/ssps/S5/ . They were produced by
reading the BMP renderings and applying the ProSMoS query-matrix format
conventions, validated for internal consistency against Chalam Chitturi's
CG-2012 IA.txt for the one byte-identical case (panel 14, 5-141-7-7).

## Why these exist at all

Before `~/chalam/CG-2012/` was restored from archive (May 2026), the only
source of canonical SSP query matrices was the website's BMP diagrams.
These translations were the best-effort attempt at recovering the design
targets from those BMPs.

## Why they should not be used for publishable work

For 13 of 14 design-target BMPs, the website's BMPs are visually distinct
from the byte-identical CG-2012 BMPs (different file sizes, different
renderings), and the underlying CG-2012 matrices imply topologies
inconsistent with my BMP-based reading. The 2014 website used a later
enumeration than CG-2012 — the third- and fourth-component sub-codes in
the website filenames (`-1-2`, `-7-7`) are not present in CG-2012, and
the bytes of `5-269-0-0.bmp` differ between CG-2012 and the website.

Specifically:
- For 5 of the 14 (`5-283-1-2`, `5-307-1-2`, `5-243-1-2`, `5-265-7-7`,
  `5-234-7-7`), the website's third-component value doesn't exist in any
  CG-2012 block.
- For 8 of the 14 (panels 1-8), CG-2012 has a same-ID block but with a
  different physical topology than the website BMP renders.
- Only `5-141-7-7` (panel 14) is byte-identical between CG-2012 and the
  website.

## What's next

`~/dev/prosmos_cl/enum/` is implementing a fresh paper-traceable SSP
enumerator (Chitturi 2016) that will output ProSMoS queries directly.
Once the enumeration is complete and validated against CG-2012 counts at
S3-S5 plus the post-2012 refinements that bring 1472 → 1239, the 14
design targets should be regenerated from that pipeline and these files
deleted.

Files retained for now because:
- they document the BMP-reading convention that turned out to match
  CG-2012's encoding (length-on-every-SSE, `sheetS` keyword, etc.);
- they were used to validate the patched `searchmatrix` parser and the
  sheetbug-cwd workaround;
- the SLURM array runner in `scripts/slurm_search/` was developed
  against this query set.
