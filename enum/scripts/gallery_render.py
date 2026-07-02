#!/usr/bin/env python3
"""Phase B (SLURM array task): render this task's slice of exemplar structures.
Run:  NTASKS=$N SLURM_ARRAY_TASK_ID=$i pymol -cq render_gallery.py
Reads staged mmCIF from s5_gallery/struct/, matched SSE ranges from the hit files.
Idempotent: skips renders that already exist.
"""
import json, os, re
from pathlib import Path
from pymol import cmd

GAL = Path("/home/rschaeff/work/prosmos_2026/s5_gallery")
STR = GAL / "struct"
OUT = GAL / "renders"; OUT.mkdir(exist_ok=True)
HITROOTS = [Path.home() / "work/prosmos_2026/ecod_search_v4/hits",
            Path.home() / "work/prosmos_2026/ecod_search_v4_retry/hits"]
SEG_COLORS = ["marine", "green", "yellow", "orange", "red"]
SIDE = 400

seg_re = re.compile(r"segment-Type:\s*(\S+)\s+Position:\s*\d+\s+Range:\s*(\d+)\s*--\s*(\d+)\s+(\S+)")


def read_segs(sk, ty, euid):
    fn = f"pdb{euid:09d}.txt"
    for rd in HITROOTS:
        p = rd / f"s5-{sk:04d}-{ty:04d}" / fn
        if p.exists():
            segs = []
            for line in p.read_text().splitlines():
                m = seg_re.match(line)
                if m:
                    segs.append((m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)))
                elif line.strip() == "END" and segs:
                    break
            return segs[:5]
    return None


def main():
    N = int(os.environ.get("NTASKS", "1"))
    T = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    spec = json.load((GAL / "full_spec.json").open())
    jobs = []
    for c in spec["cells"]:
        for e in c["exemplars"]:
            jobs.append((c["sk"], c["ty"], e["euid"], e["pdb"], e["chain"], e["did"]))
    jobs.sort(key=lambda j: (j[3], j[0], j[1]))  # by pdb -> reuse loaded structure
    sz = (len(jobs) + N - 1) // N
    mine = jobs[T * sz:(T + 1) * sz]

    cmd.set("ray_opaque_background", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 1)
    cmd.bg_color("white")

    loaded = None
    done = miss = 0
    for sk, ty, euid, pdb, chain, did in mine:
        outpng = OUT / f"{sk:04d}_{ty:02d}_{did}.png"
        if outpng.exists():
            done += 1; continue
        segs = read_segs(sk, ty, euid)
        if not segs:
            miss += 1; continue
        cif = STR / f"{pdb}.cif"
        if not cif.exists():
            miss += 1; continue
        if loaded != pdb:
            cmd.delete("all")
            try:
                cmd.load(str(cif), "m")
            except Exception:
                miss += 1; loaded = None; continue
            cmd.remove("solvent")
            loaded = pdb
        ch = chain or segs[0][3]
        lo = min(s[1] for s in segs); hi = max(s[2] for s in segs)
        dsel = f"m and chain {ch} and resi {lo-8}-{hi+8}"
        cmd.hide("everything")
        if cmd.count_atoms(dsel) == 0:
            # fall back to whole chain if the range/chain didn't resolve
            dsel = f"m and chain {ch}"
            if cmd.count_atoms(dsel) == 0:
                miss += 1; continue
        cmd.show("cartoon", dsel)
        cmd.color("grey70", dsel)
        for i, (typ, a, b, sch) in enumerate(segs):
            cmd.color(SEG_COLORS[i % 5], f"m and chain {sch} and resi {a}-{b}")
        cmd.orient(dsel)
        cmd.turn("y", 12)
        cmd.ray(SIDE, SIDE)
        cmd.png(str(outpng), dpi=110)
        done += 1
    print(f"task {T}/{N}: {len(mine)} jobs, done {done}, miss {miss}")


main()
