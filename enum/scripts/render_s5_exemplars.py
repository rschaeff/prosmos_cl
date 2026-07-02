#!/usr/bin/env python3
"""PyMOL driver: render each exemplar domain cartoon (grey) with its 5 matched
SSEs highlighted rainbow N->C.  Run with:  pymol -cq render_pymol.py -- [only_uni]

Structures pulled from the local /usr2/pdb mmCIF mirror.
"""
import json, gzip, os, sys, shutil, tempfile
from pathlib import Path
from pymol import cmd

SP = Path("/tmp/claude-1219/-home-rschaeff-dev-prosmos-cl/02f20625-920b-46f1-a2bf-bc06d84727af/scratchpad")
OUT = SP / "renders"; OUT.mkdir(exist_ok=True)
MIRROR = Path("/usr2/pdb/data/structures/divided/mmCIF")
SPEC = json.load((SP / "render_spec.json").open())
SEG_COLORS = ["marine", "green", "yellow", "orange", "red"]  # SSE 1..5, N->C


def structure_path(pdb):
    return MIRROR / pdb[1:3] / f"{pdb}.cif.gz"


def load_pdb(pdb, obj):
    src = structure_path(pdb)
    if not src.exists():
        return False
    tmp = tempfile.NamedTemporaryFile(suffix=".cif", delete=False)
    with gzip.open(src, "rb") as fi:
        shutil.copyfileobj(fi, tmp)
    tmp.close()
    cmd.load(tmp.name, obj)
    os.unlink(tmp.name)
    return True


def render(spec, outpng, side=520):
    """spec: {pdb,chain,pdb_range,segments}. Returns True on success."""
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("cartoon_transparency", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    obj = "m"
    if not load_pdb(spec["pdb"], obj):
        print("  MISSING structure", spec["pdb"]); return False
    cmd.hide("everything")
    cmd.remove("solvent")
    chain = spec["chain"] or (spec["segments"][0][3] if spec["segments"] else "A")
    # domain selection = matched-segment span within the chain (robust vs range parse)
    segs = spec["segments"]
    lo = min(s[1] for s in segs); hi = max(s[2] for s in segs)
    dsel = f"{obj} and chain {chain} and resi {lo-8}-{hi+8}"
    cmd.show("cartoon", dsel)
    cmd.color("grey70", dsel)
    cmd.set("cartoon_side_chain_helper", 1)
    # highlight matched SSEs rainbow N->C
    for i, (typ, a, b, ch) in enumerate(segs):
        sel = f"{obj} and chain {ch} and resi {a}-{b}"
        cmd.color(SEG_COLORS[i % 5], sel)
    cmd.orient(dsel)
    cmd.turn("y", 15)
    cmd.ray(side, side)
    cmd.png(str(outpng), dpi=150)
    return True


def main():
    only_uni = "only_uni" in sys.argv
    manifest = []
    for j in SPEC["unitypical"]:
        e = j["exemplar"]
        name = f"uni_sk{j['sk']:04d}_ty{j['ty']:02d}_{e['did']}.png"
        ok = render(e, OUT / name)
        manifest.append({"kind": "uni", "png": name, "ok": ok, **j})
        print(f"uni sk{j['sk']:04d} ty{j['ty']:02d} {j['typing']} {e['did']} -> {ok}")
    if not only_uni:
        for j in SPEC["promiscuous"]:
            for e in j["exemplars"]:
                name = f"pro_sk{j['sk']:04d}_ty{j['ty']:02d}_{e['group']}_{e['did']}.png"
                ok = render(e, OUT / name)
                e["png"] = name; e["ok"] = ok
                print(f"pro sk{j['sk']:04d} ty{j['ty']:02d} {e['group']:>10} {e['did']} -> {ok}")
            manifest.append({"kind": "pro", **j})
    (SP / "render_manifest.json").write_text(json.dumps(SPEC if not only_uni else {"unitypical": SPEC["unitypical"]}, indent=1))
    print("done")


main()
