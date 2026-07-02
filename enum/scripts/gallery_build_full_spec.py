#!/usr/bin/env python3
"""Phase A: build the full per-cell render spec for ALL occupied S5 cells.
  unitypical (nT==1) -> 1 exemplar ; promiscuous -> top-K experimental exemplars.
Emits full_spec.json (+ the distinct PDB list to stage).
"""
import json, re
from pathlib import Path
import psycopg2

GAL = Path("/home/rschaeff/work/prosmos_2026/s5_gallery")
PROM = Path("/home/rschaeff/dev/prosmos_cl/enum/docs/figures/s5_promiscuity.json")
K = 3

D = json.load(PROM.open())
row = D["rowSkeleton"]


def tystr(ty):
    return "".join("E" if (ty >> b) & 1 else "H" for b in range(4, -1, -1))


# collect needed exemplar domains
cells = []
need = set()
for k, v in D["cellGroups"].items():
    r, ty = map(int, k.split(","))
    sk = row[r]
    nT = len(v["T"])
    if nT == 1:
        view = "unitypical"
        ex = v["T"][:1]
    else:
        view = "promiscuous"
        ex = [g for g in v["T"] if g[2].startswith("e")][:K]
    ex = [g for g in ex if g[2].startswith("e")]
    if not ex:
        continue
    for g in ex:
        need.add(g[2])
    cells.append({"sk": sk, "ty": ty, "typing": tystr(ty), "view": view,
                  "nT": nT, "nH": len(v["H"]), "nhit": v["nhit"], "ex": ex})

# pdb_range/chain for exemplar domains
conn = psycopg2.connect(host="dione", port=45000, dbname="ecod_protein", user="ecod")
cur = conn.cursor()
cur.execute("select ecod_domain_id,chain_id,pdb_range from ecod_rep.domain "
            "where ecod_domain_id = any(%s)", (list(need),))
info = {d: (c, pr) for d, c, pr in cur.fetchall()}

pdbs = set()
out_cells = []
for c in cells:
    exs = []
    for gid, cnt, did, euid in c["ex"]:
        chain, pr = info.get(did, (None, None))
        pdb = did[1:5]
        pdbs.add(pdb)
        exs.append({"did": did, "euid": euid, "pdb": pdb, "chain": chain,
                    "group": gid, "count": cnt})
    out_cells.append({**{k: c[k] for k in ("sk", "ty", "typing", "view", "nT", "nH", "nhit")},
                      "exemplars": exs})

spec = {"K": K, "cells": out_cells, "pdbs": sorted(pdbs)}
(GAL / "full_spec.json").write_text(json.dumps(spec))
(GAL / "pdbs.txt").write_text("\n".join(sorted(pdbs)) + "\n")
n_uni = sum(1 for c in out_cells if c["view"] == "unitypical")
n_pro = sum(1 for c in out_cells if c["view"] == "promiscuous")
n_rend = sum(len(c["exemplars"]) for c in out_cells)
print(f"cells: {len(out_cells)}  unitypical {n_uni}  promiscuous {n_pro}")
print(f"render jobs: {n_rend}  distinct PDBs to stage: {len(pdbs)}")
