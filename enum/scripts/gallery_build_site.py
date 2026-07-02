#!/usr/bin/env python3
"""Phase D: build the lightweight gallery website.
  index.html : promiscuity matrix (toggle T/H) -> click a cell to load its
  composite figure from cells/.  Embeds only the small grid data.
"""
import json
from pathlib import Path

GAL = Path("/home/rschaeff/work/prosmos_2026/s5_gallery")
SITE = GAL / "site"; SITE.mkdir(exist_ok=True)
PROM = Path("/home/rschaeff/dev/prosmos_cl/enum/docs/figures/s5_promiscuity.json")


def main():
    D = json.load(PROM.open())
    have = {p.stem for p in (GAL / "cells").glob("*.png")}  # "SKKK_TT"
    row = D["rowSkeleton"]
    present = [[0] * D["nty"] for _ in range(D["nsk"])]
    for r in range(D["nsk"]):
        sk = row[r]
        for ty in range(D["nty"]):
            if f"{sk:04d}_{ty:02d}" in have:
                present[r][ty] = 1
    embed = {
        "nsk": D["nsk"], "nty": D["nty"], "rowSkeleton": row,
        "nH": D["nH"], "nT": D["nT"], "nHits": D["nHits"],
        "present": present, "stats": D["stats"],
    }
    payload = json.dumps(embed, separators=(",", ":"))
    (SITE / "index.html").write_text(TPL.replace("__DATA__", payload))
    print(f"wrote site/index.html ({len(payload)} bytes), cells present: "
          f"{sum(sum(r) for r in present)}")


TPL = r"""<!doctype html>
<html><head><meta charset="utf-8"><title>S5 matrix gallery — ECOD exemplars</title>
<style>
 body{font:13px/1.45 system-ui,sans-serif;margin:0;background:#0d1117;color:#e6edf3}
 header{padding:12px 20px;border-bottom:1px solid #30363d}
 h1{font-size:15px;margin:0 0 3px}.sub{color:#8b949e;font-size:12px}
 main{display:flex;gap:20px;padding:16px 20px;align-items:flex-start}
 .panel{background:#161b22;border:1px solid #30363d;border-radius:8px;padding:12px}
 canvas{image-rendering:pixelated;display:block;cursor:crosshair}
 .toolbar{display:flex;gap:12px;align-items:center;margin-bottom:10px;flex-wrap:wrap}
 button{background:#21262d;color:#e6edf3;border:1px solid #30363d;border-radius:6px;padding:5px 11px;cursor:pointer}
 button.on{background:#1f6feb;border-color:#1f6feb}
 .legend{display:flex;gap:2px;align-items:center;font-size:11px;color:#8b949e}
 .lg{width:13px;height:12px;display:inline-block}
 #side{width:720px;min-height:560px}
 #side img{max-width:100%;border-radius:6px;background:#fff}
 .axis{font-size:10px;color:#8b949e;margin-top:4px}
 .hint{color:#6e7681;font-size:12px;margin-top:40px;text-align:center}
 .tip{position:fixed;pointer-events:none;background:#000d;border:1px solid #30363d;padding:4px 7px;border-radius:5px;font-size:11px;display:none;z-index:9}
 kbd{background:#21262d;border:1px solid #30363d;border-radius:4px;padding:0 5px}
</style></head>
<body>
<header>
 <h1>S5 skeleton×typing matrix — ECOD exemplar gallery</h1>
 <div class="sub">198 skeletons × 32 typings · ECOD manual reps · click any lit cell to see its <b>unitypical</b> exemplar or <b>promiscuous</b> montage (skeleton schematic + ProSMoS query + folds)</div>
</header>
<main>
 <div class="panel">
  <div class="toolbar">
   <span>color by:</span>
   <button id="bT" class="on">T-groups</button>
   <button id="bH">H-groups</button>
   <span class="legend" id="leg"></span>
  </div>
  <canvas id="cv" width="512" height="990"></canvas>
  <div class="axis">x: typing 0 (HHHHH) → 31 (EEEEE) · y: skeleton rank · grey = unoccupied</div>
 </div>
 <div class="panel" id="side"><div class="hint">◀ click a lit cell</div></div>
</main>
<div class="tip" id="tip"></div>
<script>
const D=__DATA__;const NSK=D.nsk,NTY=D.nty;let mode="T";
function grid(){return mode==="T"?D.nT:D.nH;}
const nHits=D.nHits,present=D.present;
function turbo(t){t=Math.max(0,Math.min(1,t));
 const s=[[48,18,59],[65,69,171],[57,118,226],[42,161,213],[36,197,158],[93,220,96],
 [176,221,47],[238,190,42],[251,133,29],[228,66,12],[153,17,7]];
 const x=t*(s.length-1),i=Math.floor(x),f=x-i,a=s[i],b=s[Math.min(i+1,s.length-1)];
 return`rgb(${a[0]+(b[0]-a[0])*f|0},${a[1]+(b[1]-a[1])*f|0},${a[2]+(b[2]-a[2])*f|0})`;}
let vmax=1;for(const r of D.nT)for(const v of r)if(v>vmax)vmax=v;const lv=Math.log(vmax);
const cv=document.getElementById("cv"),CW=Math.floor(512/NTY),CH=Math.floor(990/NSK),W=CW*NTY,H=CH*NSK;
cv.width=W;cv.height=H;const ctx=cv.getContext("2d");let sel=null;
function draw(){ctx.fillStyle="#0d1117";ctx.fillRect(0,0,W,H);const g=grid();
 for(let r=0;r<NSK;r++)for(let ty=0;ty<NTY;ty++){
  if(nHits[r][ty]<=0){ctx.fillStyle="#21262d";}
  else{const v=g[r][ty];ctx.fillStyle=v>0?turbo(Math.log(v)/lv):"#21262d";}
  ctx.fillRect(ty*CW,r*CH,CW,CH);}
 if(sel){ctx.strokeStyle="#fff";ctx.lineWidth=1.5;ctx.strokeRect(sel.ty*CW,sel.r*CH,CW,CH);}}
function legend(){document.getElementById("leg").innerHTML=`<span>1</span>`+
 [...Array(10)].map((_,k)=>`<span class="lg" style="background:${turbo(k/9)}"></span>`).join("")+`<span>${vmax}</span>`;}
legend();
function pad(n,w){return String(n).padStart(w,"0");}
function typingStr(ty){let s="";for(let b=4;b>=0;b--)s+=((ty>>b)&1)?"E":"H";return s;}
function pick(r,ty){
 const side=document.getElementById("side");
 if(nHits[r][ty]<=0){side.innerHTML='<div class="hint">unoccupied cell — no hits</div>';sel=null;draw();return;}
 sel={r,ty};draw();
 const sk=D.rowSkeleton[r];
 if(!present[r][ty]){side.innerHTML=`<div class="hint">skeleton ${pad(sk,4)} · typing ${ty} (${typingStr(ty)})<br>no exemplar figure (non-experimental only)</div>`;return;}
 const img=`cells/${pad(sk,4)}_${pad(ty,2)}.png`;
 side.innerHTML=`<img src="${img}" alt="cell ${pad(sk,4)}_${pad(ty,2)}">`;
}
cv.onclick=e=>{const b=cv.getBoundingClientRect();
 const ty=Math.floor((e.clientX-b.left)/(b.width/NTY)),r=Math.floor((e.clientY-b.top)/(b.height/NSK));
 if(r>=0&&r<NSK&&ty>=0&&ty<NTY)pick(r,ty);};
const tip=document.getElementById("tip");
cv.onmousemove=e=>{const b=cv.getBoundingClientRect();
 const ty=Math.floor((e.clientX-b.left)/(b.width/NTY)),r=Math.floor((e.clientY-b.top)/(b.height/NSK));
 if(r<0||r>=NSK||ty<0||ty>=NTY||nHits[r][ty]<=0){tip.style.display="none";return;}
 tip.style.display="block";tip.style.left=(e.clientX+12)+"px";tip.style.top=(e.clientY+12)+"px";
 tip.innerHTML=`sk ${pad(D.rowSkeleton[r],4)} · ty ${ty} (${typingStr(ty)})<br>T=${D.nT[r][ty]} H=${D.nH[r][ty]} · ${nHits[r][ty]} hits`;};
cv.onmouseleave=()=>tip.style.display="none";
document.getElementById("bT").onclick=()=>{mode="T";bT.classList.add("on");bH.classList.remove("on");draw();};
document.getElementById("bH").onclick=()=>{mode="H";bH.classList.add("on");bT.classList.remove("on");draw();};
draw();
</script></body></html>
"""

if __name__ == "__main__":
    main()
