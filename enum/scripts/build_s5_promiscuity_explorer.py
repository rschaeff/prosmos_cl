#!/usr/bin/env python3
"""Interactive promiscuity explorer: click a cell -> ECOD-group breakdown +
exemplar domains linked to ECOD/RCSB.  Embeds a trimmed (top-K groups/cell)
version of s5_promiscuity.json into a self-contained HTML.
"""
import json
from pathlib import Path

SP = Path("/tmp/claude-1219/-home-rschaeff-dev-prosmos-cl/02f20625-920b-46f1-a2bf-bc06d84727af/scratchpad")
TOPK = 20


def main():
    D = json.load((SP / "s5_promiscuity.json").open())
    cg = {}
    for k, v in D["cellGroups"].items():
        cg[k] = {
            "nhit": v["nhit"],
            "H": v["H"][:TOPK],
            "T": v["T"][:TOPK],
            "nH": len(v["H"]),
            "nT": len(v["T"]),
        }
    embed = {
        "nsk": D["nsk"], "nty": D["nty"],
        "rowSkeleton": D["rowSkeleton"],
        "nH": D["nH"], "nT": D["nT"], "nHits": D["nHits"],
        "cellGroups": cg,
        "stats": D["stats"],
        "topk": TOPK,
    }
    payload = json.dumps(embed, separators=(",", ":"))
    (SP / "s5_promiscuity_explorer.html").write_text(TPL.replace("__DATA__", payload))
    print(f"wrote explorer ({len(payload)} bytes payload)")


TPL = r"""<!doctype html>
<html><head><meta charset="utf-8"><title>S5 ECOD-promiscuity explorer</title>
<style>
 body{font:13px/1.45 system-ui,sans-serif;margin:0;background:#0d1117;color:#e6edf3}
 header{padding:12px 20px;border-bottom:1px solid #30363d}
 h1{font-size:15px;margin:0 0 3px}.sub{color:#8b949e;font-size:12px}
 main{display:flex;gap:20px;padding:16px 20px;align-items:flex-start}
 .panel{background:#161b22;border:1px solid #30363d;border-radius:8px;padding:12px}
 canvas{image-rendering:pixelated;display:block;cursor:crosshair}
 .toolbar{display:flex;gap:14px;align-items:center;margin-bottom:10px;flex-wrap:wrap}
 button{background:#21262d;color:#e6edf3;border:1px solid #30363d;border-radius:6px;padding:5px 11px;cursor:pointer}
 button.on{background:#1f6feb;border-color:#1f6feb}
 .legend{display:flex;gap:2px;align-items:center;font-size:11px;color:#8b949e}
 .lg{width:13px;height:12px;display:inline-block}
 #side{width:430px;min-height:600px}
 #side h2{font-size:14px;margin:0 0 4px}
 .badge{display:inline-block;padding:2px 8px;border-radius:10px;font-size:11px;font-weight:600}
 .uni{background:#0d3320;color:#3fb950;border:1px solid #238636}
 .promisc{background:#3d1c00;color:#f0883e;border:1px solid #bd561d}
 .kv{color:#8b949e}.kv b{color:#e6edf3}
 table{border-collapse:collapse;width:100%;margin-top:10px;font-size:12px}
 th,td{text-align:left;padding:3px 6px;border-bottom:1px solid #21262d}
 th{color:#8b949e;font-weight:600}
 a{color:#58a6ff;text-decoration:none}a:hover{text-decoration:underline}
 .bar{height:9px;background:#1f6feb;border-radius:2px;display:inline-block;vertical-align:middle}
 .axis{font-size:10px;color:#8b949e;margin-top:4px}
 .tip{position:fixed;pointer-events:none;background:#000d;border:1px solid #30363d;padding:4px 7px;border-radius:5px;font-size:11px;display:none;z-index:9}
 .hint{color:#6e7681;font-size:12px;margin-top:40px;text-align:center}
</style></head>
<body>
<header>
 <h1>S5 matrix &mdash; ECOD-group promiscuity explorer</h1>
 <div class="sub">198 skeletons &times; 32 typings &middot; ECOD manual reps &middot; click a cell for its ECOD-group breakdown &amp; exemplar domains</div>
</header>
<main>
 <div class="panel">
  <div class="toolbar">
   <span>color by:</span>
   <button id="bT" class="on">T-groups (topology)</button>
   <button id="bH">H-groups (homology)</button>
   <span class="legend" id="leg"></span>
  </div>
  <canvas id="cv" width="640" height="990"></canvas>
  <div class="axis">x: typing 0 (HHHHH) &rarr; 31 (EEEEE) &nbsp;|&nbsp; y: skeleton rank (shared order) &nbsp;|&nbsp; grey = unoccupied</div>
 </div>
 <div class="panel" id="side">
  <div class="hint">&#9664; click a cell to inspect its ECOD groups</div>
 </div>
</main>
<div class="tip" id="tip"></div>
<script>
const D=__DATA__;
const NSK=D.nsk,NTY=D.nty;
let mode="T";
function grid(){return mode==="T"?D.nT:D.nH;}
const nHits=D.nHits;
// turbo-ish colormap
function turbo(t){t=Math.max(0,Math.min(1,t));
 const s=[[48,18,59],[65,69,171],[57,118,226],[42,161,213],[36,197,158],[93,220,96],
 [176,221,47],[238,190,42],[251,133,29],[228,66,12],[153,17,7]];
 const x=t*(s.length-1),i=Math.floor(x),f=x-i,a=s[i],b=s[Math.min(i+1,s.length-1)];
 return`rgb(${a[0]+(b[0]-a[0])*f|0},${a[1]+(b[1]-a[1])*f|0},${a[2]+(b[2]-a[2])*f|0})`;}
let vmax=1;for(const row of D.nT)for(const v of row)if(v>vmax)vmax=v;
const lvmax=Math.log(vmax);
function col(v){return v<=0?null:turbo(Math.log(v)/lvmax);}

const cv=document.getElementById("cv");
const CW=Math.floor(640/NTY),CH=Math.floor(990/NSK),W=CW*NTY,H=CH*NSK;
cv.width=W;cv.height=H;const ctx=cv.getContext("2d");
function draw(){ctx.fillStyle="#0d1117";ctx.fillRect(0,0,W,H);const g=grid();
 for(let r=0;r<NSK;r++)for(let ty=0;ty<NTY;ty++){
  if(nHits[r][ty]<=0){ctx.fillStyle="#21262d";}
  else{const c=col(g[r][ty]);if(!c)continue;ctx.fillStyle=c;}
  ctx.fillRect(ty*CW,r*CH,CW,CH);}
 if(sel){ctx.strokeStyle="#fff";ctx.lineWidth=1.5;ctx.strokeRect(sel.ty*CW,sel.r*CH,CW,CH);}}
function legend(){document.getElementById("leg").innerHTML=
 `<span>1 (unitypical)</span>`+[...Array(10)].map((_,k)=>
 `<span class="lg" style="background:${turbo(k/9)}"></span>`).join("")+`<span>${vmax} (promiscuous)</span>`;}
legend();

function typingStr(ty){let s="";for(let b=4;b>=0;b--)s+=((ty>>b)&1)?"E":"H";return s;}
function ecodLink(did){return`http://prodata.swmed.edu/ecod/complete/domain/${did}`;}
function rcsb(did){return did[0]==="e"?`https://www.rcsb.org/structure/${did.slice(1,5).toUpperCase()}`:null;}

let sel=null;
function select(r,ty){
 sel={r,ty};draw();
 const sk=D.rowSkeleton[r], key=r+","+ty, cg=D.cellGroups[key];
 const side=document.getElementById("side");
 if(!cg){side.innerHTML=`<div class="hint">skeleton ${String(sk).padStart(4,"0")} &middot; typing ${ty} (${typingStr(ty)})<br>no hits &mdash; unoccupied cell</div>`;return;}
 const nT=cg.nT,nH=cg.nH;
 const uni = (mode==="T"?nT:nH)===1;
 const groups = mode==="T"?cg.T:cg.H;
 const ng = mode==="T"?nT:nH;
 const shown=groups.length;
 let rows="";const maxc=groups.length?groups[0][1]:1;
 for(const[g,c,did] of groups){
  const bw=Math.max(3,Math.round(70*c/maxc));
  const rl=rcsb(did);
  rows+=`<tr><td><code>${g}</code></td><td>${c}</td>`+
   `<td><a href="${ecodLink(did)}" target="_blank">${did}</a>`+
   (rl?` &middot; <a href="${rl}" target="_blank">${did.slice(1,5).toUpperCase()}</a>`:"")+`</td>`+
   `<td><span class="bar" style="width:${bw}px"></span></td></tr>`;}
 side.innerHTML=
  `<h2>skeleton ${String(sk).padStart(4,"0")} &middot; typing ${ty}</h2>`+
  `<div class="kv">H/E pattern (SSE1&rarr;5): <b>${typingStr(ty)}</b></div>`+
  `<div class="kv">hits (distinct rep domains): <b>${cg.nhit}</b> &middot; T-groups: <b>${nT}</b> &middot; H-groups: <b>${nH}</b></div>`+
  `<div style="margin:8px 0">`+
   (uni?`<span class="badge uni">UNITYPICAL &mdash; one ${mode}-group</span>`
       :`<span class="badge promisc">PROMISCUOUS &mdash; ${ng} ${mode}-groups</span>`)+`</div>`+
  `<div class="kv" style="font-size:11px">${mode==="T"?"topology":"homology"} groups${shown<ng?` (top ${shown} of ${ng})`:""}, by # rep domains:</div>`+
  `<table><tr><th>ECOD ${mode}-group</th><th>#dom</th><th>exemplar (ECOD &middot; PDB)</th><th></th></tr>${rows}</table>`+
  `<div class="axis" style="margin-top:8px">exemplar = a representative domain from that group matching this geometry; links open ECOD domain page / RCSB entry.</div>`;
}

cv.onclick=e=>{const b=cv.getBoundingClientRect();
 const ty=Math.floor((e.clientX-b.left)/(b.width/NTY)),r=Math.floor((e.clientY-b.top)/(b.height/NSK));
 if(r>=0&&r<NSK&&ty>=0&&ty<NTY)select(r,ty);};
const tip=document.getElementById("tip");
cv.onmousemove=e=>{const b=cv.getBoundingClientRect();
 const ty=Math.floor((e.clientX-b.left)/(b.width/NTY)),r=Math.floor((e.clientY-b.top)/(b.height/NSK));
 if(r<0||r>=NSK||ty<0||ty>=NTY||nHits[r][ty]<=0){tip.style.display="none";return;}
 tip.style.display="block";tip.style.left=(e.clientX+12)+"px";tip.style.top=(e.clientY+12)+"px";
 tip.innerHTML=`sk ${String(D.rowSkeleton[r]).padStart(4,"0")} &middot; ty ${ty} (${typingStr(ty)})<br>`+
  `${nHits[r][ty]} hits &middot; T=${D.nT[r][ty]} H=${D.nH[r][ty]}`;};
cv.onmouseleave=()=>tip.style.display="none";

document.getElementById("bT").onclick=()=>{mode="T";bT.classList.add("on");bH.classList.remove("on");draw();if(sel)select(sel.r,sel.ty);};
document.getElementById("bH").onclick=()=>{mode="H";bH.classList.add("on");bT.classList.remove("on");draw();if(sel)select(sel.r,sel.ty);};
draw();
</script></body></html>
"""

if __name__ == "__main__":
    main()
