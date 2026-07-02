#!/usr/bin/env python3
"""Build a time-series of the S5 (198 skeleton x 32 typing) ProSMoS hit matrix
against the ECOD manual-representative DB, binned by PDB *deposition* year.

Question it answers: did the region of skeleton x typing space that the PDB
samples *change* (new cells lighting up in later years) over time, or did it
basically *fill in* (the same cells getting deeper while the occupied set
saturates early)?

Inputs (all in scratchpad):
  manual_hits.txt      lines: hits/s5-SK-TY/pdbNNNNNNNNN.txt   (one per hit domain)
  ecodrep_uid_map.tsv  uid \t ecod_uid \t pdb \t ecod_domain_id
  pdb2date.tsv         pdb \t MM/DD/YY   (RCSB accession/deposition date)

Output:
  s5_timeseries.json   the year cube + first-appearance grid (record/debug)
  s5_timeseries.html   self-contained interactive viz (slider + first-appearance)
"""
import json
import re
from pathlib import Path
from collections import defaultdict

SP = Path("/tmp/claude-1219/-home-rschaeff-dev-prosmos-cl/02f20625-920b-46f1-a2bf-bc06d84727af/scratchpad")
HITS = SP / "manual_hits.txt"
UIDMAP = SP / "ecodrep_uid_map.tsv"
PDBDATE = SP / "pdb2date.tsv"
OUT_JSON = SP / "s5_timeseries.json"
OUT_HTML = SP / "s5_timeseries.html"

NSK, NTY = 198, 32


def parse_year(mmddyy: str) -> int | None:
    m = re.match(r"(\d{2})/(\d{2})/(\d{2})", mmddyy.strip())
    if not m:
        return None
    yy = int(m.group(3))
    # PDB spans 1972..present; pivot at 30 (nothing deposited 2031+ yet).
    return 1900 + yy if yy > 30 else 2000 + yy


def main():
    # ---- pdb -> year ----
    pdb_year = {}
    for line in PDBDATE.open():
        pdb, date = line.rstrip("\n").split("\t")
        y = parse_year(date)
        if y is not None:
            pdb_year[pdb.lower()] = y
    print(f"pdb_year: {len(pdb_year)}")

    # ---- uid/ecod_uid -> pdb ----
    uid2pdb, euid2pdb = {}, {}
    for line in UIDMAP.open():
        uid, euid, pdb, _did = line.rstrip("\n").split("\t")
        if pdb:
            if uid:
                uid2pdb[int(uid)] = pdb.lower()
            if euid:
                euid2pdb[int(euid)] = pdb.lower()
    print(f"uid2pdb: {len(uid2pdb)}  euid2pdb: {len(euid2pdb)}")

    # ---- decide which key the hit files use (coverage probe) ----
    sample, tot = 0, 0
    hit_re = re.compile(r"s5-(\d{4})-(\d{4})/pdb0*(\d+)\.txt$")
    with HITS.open() as f:
        for line in f:
            m = hit_re.search(line)
            if not m:
                continue
            tot += 1
            if tot > 20000:
                break
            uidn = int(m.group(3))
            if uidn in uid2pdb:
                sample |= 1
            if uidn in euid2pdb:
                sample |= 2
    hit_uid, hit_euid = 0, 0
    with HITS.open() as f:
        for line in f:
            m = hit_re.search(line)
            if not m:
                continue
            uidn = int(m.group(3))
            hit_uid += uidn in uid2pdb
            hit_euid += uidn in euid2pdb
            if hit_uid + hit_euid > 40000:
                break
    key = "uid" if hit_uid >= hit_euid else "ecod_uid"
    keymap = uid2pdb if key == "uid" else euid2pdb
    print(f"key={key}  (uid hits {hit_uid} vs ecod_uid hits {hit_euid})")

    # ---- accumulate per-cell year increments ----
    # cell = (sk, ty).  inc[cell][year] = count of hit domains deposited that year
    inc = defaultdict(lambda: defaultdict(int))
    n_hits = n_mapped = n_dated = 0
    unmapped_uid = set()
    with HITS.open() as f:
        for line in f:
            m = hit_re.search(line)
            if not m:
                continue
            n_hits += 1
            sk, ty, uidn = int(m.group(1)), int(m.group(2)), int(m.group(3))
            pdb = keymap.get(uidn)
            if pdb is None:
                unmapped_uid.add(uidn)
                continue
            n_mapped += 1
            y = pdb_year.get(pdb)
            if y is None:
                continue
            n_dated += 1
            inc[(sk, ty)][y] += 1
    print(f"hits {n_hits}  mapped {n_mapped}  dated {n_dated}  "
          f"unmapped_uids {len(unmapped_uid)}")

    # ---- year axis ----
    all_years = sorted({y for cell in inc.values() for y in cell})
    y0, y1 = all_years[0], all_years[-1]
    years = list(range(y0, y1 + 1))
    yidx = {y: i for i, y in enumerate(years)}
    print(f"years {y0}..{y1} ({len(years)})")

    # ---- per-cell final total + first-appearance year ----
    final_total = {}
    first_year = {}
    for cell, ymap in inc.items():
        final_total[cell] = sum(ymap.values())
        first_year[cell] = min(ymap)

    # ---- fixed row order: skeletons by final total desc (stable across frames) ----
    sk_total = defaultdict(int)
    for (sk, ty), t in final_total.items():
        sk_total[sk] += t
    row_order = sorted(range(NSK), key=lambda s: (-sk_total[s], s))
    row_pos = {sk: i for i, sk in enumerate(row_order)}

    # ---- sparse series for JS: list of [rowPos, ty, yearIdx, inc] ----
    series = []
    for (sk, ty), ymap in inc.items():
        r = row_pos[sk]
        for y, c in sorted(ymap.items()):
            series.append([r, ty, yidx[y], c])

    # first-appearance grid in display order: firstGrid[rowPos*32+ty] = yearIdx or -1
    first_grid = [-1] * (NSK * NTY)
    for (sk, ty), y in first_year.items():
        first_grid[row_pos[sk] * NTY + ty] = yidx[y]

    # ---- per-year aggregate curves: occupied cells & cumulative total hits ----
    occ_by_year = [0] * len(years)   # cells with >=1 hit by year Y (cumulative)
    tot_by_year = [0] * len(years)   # total hits deposited by year Y (cumulative)
    # occupancy: a cell counts from its first_year onward
    for cell, fy in first_year.items():
        for i in range(yidx[fy], len(years)):
            occ_by_year[i] += 1
    for cell, ymap in inc.items():
        run = 0
        for i, y in enumerate(years):
            run += ymap.get(y, 0)
            tot_by_year[i] += run

    n_occupied = len(final_total)
    data = {
        "nsk": NSK, "nty": NTY,
        "years": years,
        "series": series,               # [rowPos, ty, yearIdx, inc]
        "firstGrid": first_grid,        # rowPos*32+ty -> yearIdx | -1
        "rowSkeleton": row_order,       # display row -> original skeleton id
        "occByYear": occ_by_year,
        "totByYear": tot_by_year,
        "nOccupied": n_occupied,
        "stats": {
            "hits": n_hits, "mapped": n_mapped, "dated": n_dated,
            "occupiedCells": n_occupied, "totalCells": NSK * NTY,
            "key": key,
        },
    }
    OUT_JSON.write_text(json.dumps(data))
    print(f"wrote {OUT_JSON}  ({OUT_JSON.stat().st_size} bytes)")
    print(f"occupied cells (final): {n_occupied}/{NSK*NTY}")
    print(f"occupancy first year {years[0]}: {occ_by_year[0]}  "
          f"last year {years[-1]}: {occ_by_year[-1]}")

    write_html(data)
    print(f"wrote {OUT_HTML}")


def write_html(data):
    payload = json.dumps(data, separators=(",", ":"))
    OUT_HTML.write_text(HTML_TEMPLATE.replace("__DATA__", payload))


HTML_TEMPLATE = r"""<!doctype html>
<html><head><meta charset="utf-8"><title>S5 PDB-deposition time series</title>
<style>
  body{font:13px/1.4 system-ui,sans-serif;margin:0;background:#0d1117;color:#e6edf3}
  header{padding:14px 20px;border-bottom:1px solid #30363d}
  h1{font-size:16px;margin:0 0 4px}
  .sub{color:#8b949e;font-size:12px}
  main{display:flex;gap:26px;padding:18px 20px;flex-wrap:wrap}
  .panel{background:#161b22;border:1px solid #30363d;border-radius:8px;padding:12px}
  .panel h2{font-size:12px;margin:0 0 8px;color:#c9d1d9;font-weight:600}
  canvas{image-rendering:pixelated;display:block}
  .ctrl{display:flex;align-items:center;gap:12px;padding:10px 20px;border-top:1px solid #30363d;border-bottom:1px solid #30363d;flex-wrap:wrap}
  input[type=range]{width:420px}
  button{background:#21262d;color:#e6edf3;border:1px solid #30363d;border-radius:6px;padding:5px 12px;cursor:pointer}
  button:hover{background:#30363d}
  #yr{font-size:22px;font-weight:700;min-width:70px;text-align:center}
  .legend{display:flex;gap:2px;align-items:center;font-size:11px;color:#8b949e}
  .lg{width:14px;height:12px;display:inline-block}
  .stat{color:#8b949e;font-size:12px}
  .stat b{color:#e6edf3}
  .axis{font-size:10px;color:#8b949e}
  .tip{position:fixed;pointer-events:none;background:#000c;border:1px solid #30363d;padding:4px 7px;border-radius:5px;font-size:11px;display:none;z-index:9}
</style></head>
<body>
<header>
  <h1>S5 skeleton &times; H/E typing matrix &mdash; sampled by PDB over deposition time</h1>
  <div class="sub">ECOD manual-representative DB &middot; 198 skeletons &times; 32 typings = 6,336 cells &middot; rows sorted by final hit total (fixed across all years)</div>
</header>
<div class="ctrl">
  <button id="play">&#9654; Play</button>
  <input type="range" id="slider" min="0" max="0" value="0" step="1">
  <div id="yr">&mdash;</div>
  <div class="stat">occupied cells &le; year: <b id="s_occ">0</b>/<span id="s_tot">0</span>
    &nbsp;&middot;&nbsp; cumulative hits: <b id="s_hits">0</b></div>
</div>
<main>
  <div class="panel">
    <h2>Cumulative hit count &le; selected year <span id="cm_leg" class="legend"></span></h2>
    <canvas id="cum" width="640" height="990"></canvas>
    <div class="axis">x: typing 0 (HHHHH) &rarr; 31 (EEEEE) &nbsp;|&nbsp; y: skeleton rank</div>
  </div>
  <div class="panel">
    <h2>First-appearance year (when each cell first got a hit) <span id="fa_leg" class="legend"></span></h2>
    <canvas id="first" width="640" height="990"></canvas>
    <div class="axis">grey = never occupied &middot; same row order as left</div>
  </div>
  <div class="panel">
    <h2>Occupancy vs. total-hits over time</h2>
    <canvas id="curve" width="360" height="300"></canvas>
    <div class="axis" id="curve_note"></div>
    <p class="stat" style="max-width:340px;margin-top:10px" id="verdict"></p>
  </div>
</main>
<div class="tip" id="tip"></div>
<script>
const D=__DATA__;
const NSK=D.nsk,NTY=D.nty,YEARS=D.years,NY=YEARS.length;
// ---- viridis-ish colormap (t in 0..1) ----
function viridis(t){t=Math.max(0,Math.min(1,t));
 const s=[[68,1,84],[72,40,120],[62,74,137],[49,104,142],[38,130,142],
 [31,158,137],[53,183,121],[110,206,88],[181,222,43],[253,231,37]];
 const x=t*(s.length-1),i=Math.floor(x),f=x-i,a=s[i],b=s[Math.min(i+1,s.length-1)];
 return`rgb(${a[0]+(b[0]-a[0])*f|0},${a[1]+(b[1]-a[1])*f|0},${a[2]+(b[2]-a[2])*f|0})`;}
// plasma-ish for years
function plasma(t){t=Math.max(0,Math.min(1,t));
 const s=[[13,8,135],[84,2,163],[139,10,165],[185,50,137],[219,92,104],
 [244,136,73],[254,188,43],[240,249,33]];
 const x=t*(s.length-1),i=Math.floor(x),f=x-i,a=s[i],b=s[Math.min(i+1,s.length-1)];
 return`rgb(${a[0]+(b[0]-a[0])*f|0},${a[1]+(b[1]-a[1])*f|0},${a[2]+(b[2]-a[2])*f|0})`;}

// cell geometry
const CW=Math.floor(640/NTY), CH=Math.floor(990/NSK);
const W=CW*NTY, H=CH*NSK;
for(const id of["cum","first"]){const c=document.getElementById(id);c.width=W;c.height=H;}

// per-cell cumulative counts, rebuilt as slider moves
const cum=new Float64Array(NSK*NTY);
// group series by yearIdx for incremental update
const byYear=Array.from({length:NY},()=>[]);
for(const[r,ty,yi,inc]of D.series) byYear[yi].push([r*NTY+ty,inc]);
let curYear=-1, maxCum=1;
for(const[r,ty,yi,inc]of D.series){} // maxCum computed after full accumulate below
// precompute global max cumulative (at final year) for stable color scale
{const tmp=new Float64Array(NSK*NTY);
 for(const[r,ty,yi,inc]of D.series)tmp[r*NTY+ty]+=inc;
 for(const v of tmp)if(v>maxCum)maxCum=v;}
const logMax=Math.log10(maxCum+1);

const ctxC=document.getElementById("cum").getContext("2d");
const ctxF=document.getElementById("first").getContext("2d");

function drawCum(){
 ctxC.fillStyle="#0d1117";ctxC.fillRect(0,0,W,H);
 for(let i=0;i<cum.length;i++){const v=cum[i];if(v<=0)continue;
  const r=(i/NTY)|0,ty=i%NTY;
  ctxC.fillStyle=viridis(Math.log10(v+1)/logMax);
  ctxC.fillRect(ty*CW,r*CH,CW,CH);}
}
function setYear(yi){
 if(yi>curYear){for(let k=curYear+1;k<=yi;k++)for(const[idx,inc]of byYear[k])cum[idx]+=inc;}
 else{cum.fill(0);for(let k=0;k<=yi;k++)for(const[idx,inc]of byYear[k])cum[idx]+=inc;}
 curYear=yi;drawCum();
 let occ=0,hits=0;for(const v of cum){if(v>0)occ++;hits+=v;}
 document.getElementById("yr").textContent=YEARS[yi];
 document.getElementById("s_occ").textContent=occ;
 document.getElementById("s_hits").textContent=hits.toLocaleString();
 slider.value=yi;drawCurveMarker(yi);
}

// first-appearance grid (static)
function drawFirst(){
 ctxF.fillStyle="#0d1117";ctxF.fillRect(0,0,W,H);
 for(let i=0;i<D.firstGrid.length;i++){const yi=D.firstGrid[i];
  const r=(i/NTY)|0,ty=i%NTY;
  if(yi<0){ctxF.fillStyle="#21262d";}
  else ctxF.fillStyle=plasma(yi/(NY-1));
  ctxF.fillRect(ty*CW,r*CH,CW,CH);}
}
drawFirst();

// legends
function legend(el,fn,lab0,lab1){el.innerHTML=
 `<span>${lab0}</span>`+[...Array(10)].map((_,k)=>
 `<span class="lg" style="background:${fn(k/9)}"></span>`).join("")+`<span>${lab1}</span>`;}
legend(document.getElementById("cm_leg"),viridis,"1",maxCum.toLocaleString());
legend(document.getElementById("fa_leg"),plasma,YEARS[0],YEARS[NY-1]);

// occupancy / total curves
const curve=document.getElementById("curve"),cx=curve.getContext("2d");
const CWd=curve.width,CHd=curve.height,PAD=34;
const occ=D.occByYear,tot=D.totByYear;
const occMax=Math.max(...occ),totMax=Math.max(...tot);
function px(i){return PAD+(CWd-PAD-8)*i/(NY-1);}
function pyOcc(v){return CHd-PAD-(CHd-PAD-8)*v/occMax;}
function pyTot(v){return CHd-PAD-(CHd-PAD-8)*v/totMax;}
let markX=null;
function drawCurve(){
 cx.fillStyle="#161b22";cx.fillRect(0,0,CWd,CHd);
 cx.strokeStyle="#30363d";cx.lineWidth=1;
 cx.beginPath();cx.moveTo(PAD,8);cx.lineTo(PAD,CHd-PAD);cx.lineTo(CWd-8,CHd-PAD);cx.stroke();
 // total hits (grey, right-implied)
 cx.strokeStyle="#8b949e";cx.lineWidth=1.5;cx.beginPath();
 tot.forEach((v,i)=>{const X=px(i),Y=pyTot(v);i?cx.lineTo(X,Y):cx.moveTo(X,Y);});cx.stroke();
 // occupancy (green)
 cx.strokeStyle="#3fb950";cx.lineWidth=2;cx.beginPath();
 occ.forEach((v,i)=>{const X=px(i),Y=pyOcc(v);i?cx.lineTo(X,Y):cx.moveTo(X,Y);});cx.stroke();
 // labels
 cx.fillStyle="#3fb950";cx.font="10px sans-serif";cx.fillText("occupied cells",PAD+4,16);
 cx.fillStyle="#8b949e";cx.fillText("cumulative hits",PAD+4,30);
 cx.fillText(YEARS[0],PAD-6,CHd-PAD+12);cx.fillText(YEARS[NY-1],CWd-30,CHd-PAD+12);
 if(markX!==null){cx.strokeStyle="#f0883e";cx.lineWidth=1;cx.beginPath();
  cx.moveTo(markX,8);cx.lineTo(markX,CHd-PAD);cx.stroke();}
}
function drawCurveMarker(yi){markX=px(yi);drawCurve();}
drawCurve();

// verdict text
{const occEnd=occ[NY-1],occHalf=occ[occ.findIndex((_,i)=>YEARS[i]>= (YEARS[0]+YEARS[NY-1])/2)]||0;
 const half=Math.round((YEARS[0]+YEARS[NY-1])/2);
 const idxHalf=YEARS.indexOf(half)>=0?YEARS.indexOf(half):Math.floor(NY/2);
 const fracHalf=(occ[idxHalf]/occEnd*100).toFixed(0);
 document.getElementById("verdict").innerHTML=
  `By <b>${YEARS[idxHalf]}</b> the PDB already occupied <b>${fracHalf}%</b> of the cells it ever reaches (${occEnd} of ${D.stats.totalCells}). `+
  `If the green curve saturates early while the grey (hit count) keeps climbing, the matrix space <b>filled in</b> rather than <b>changed</b>.`;}
document.getElementById("s_tot").textContent=D.stats.totalCells.toLocaleString();

// controls
const slider=document.getElementById("slider");slider.max=NY-1;slider.value=NY-1;
slider.oninput=e=>setYear(+e.target.value);
let playing=false,timer=null;
document.getElementById("play").onclick=function(){
 playing=!playing;this.innerHTML=playing?"&#10073;&#10073; Pause":"&#9654; Play";
 if(playing){let yi=(+slider.value>=NY-1)?0:+slider.value;
  timer=setInterval(()=>{setYear(yi);if(yi>=NY-1){clearInterval(timer);playing=false;
   document.getElementById("play").innerHTML="&#9654; Play";}yi++;},360);}
 else clearInterval(timer);};

// tooltip on cumulative canvas
const tip=document.getElementById("tip");
document.getElementById("cum").onmousemove=e=>{
 const c=e.target.getBoundingClientRect();
 const ty=Math.floor((e.clientX-c.left)/(c.width/NTY));
 const r=Math.floor((e.clientY-c.top)/(c.height/NSK));
 if(r<0||r>=NSK||ty<0||ty>=NTY){tip.style.display="none";return;}
 const sk=D.rowSkeleton[r],v=cum[r*NTY+ty],fa=D.firstGrid[r*NTY+ty];
 tip.style.display="block";tip.style.left=(e.clientX+12)+"px";tip.style.top=(e.clientY+12)+"px";
 tip.innerHTML=`skeleton ${String(sk).padStart(4,"0")} &middot; typing ${ty}<br>`+
  `hits &le; ${YEARS[curYear]}: <b>${v}</b>`+(fa>=0?`<br>first seen: ${YEARS[fa]}`:`<br>never occupied`);
};
document.getElementById("cum").onmouseleave=()=>tip.style.display="none";

setYear(NY-1);
</script>
</body></html>
"""


if __name__ == "__main__":
    main()
