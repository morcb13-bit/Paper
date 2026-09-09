#  gen_html.py ── well_rt.html を書き出す
import json

DATA = open("rt_data.json").read()

HTML = r"""<!DOCTYPE html>
<meta charset="utf-8">
<title>井戸のかたち ── 黄金のひし形30面体の上を走る整数のオートマトン</title>
<style>
  :root{
    --paper:#eceef2; --ink:#232830; --hair:#5b5b6a;
    --wall:#1b2540; --plus:#a8452e; --minus:#26417f;
  }
  *{box-sizing:border-box}
  body{
    margin:0; background:var(--paper); color:var(--ink);
    font-family:"Hiragino Mincho ProN","Yu Mincho",YuMincho,"Noto Serif JP",serif;
    font-feature-settings:"palt"; -webkit-font-smoothing:antialiased;
  }
  .wrap{max-width:780px; margin:0 auto; padding:28px 20px 56px}
  h1{font-size:19px; font-weight:600; letter-spacing:.04em; margin:0 0 2px}
  .lede{font-size:13.5px; line-height:1.75; color:#4b515c; margin:0 0 20px; max-width:60ch}
  canvas{display:block; width:100%; height:auto; cursor:grab; touch-action:none; background:#0d1018}
  canvas:active{cursor:grabbing}
  .readout{
    display:flex; gap:26px; flex-wrap:wrap; align-items:baseline;
    border-top:1px solid var(--hair); border-bottom:1px solid var(--hair);
    padding:9px 2px; margin:6px 0 16px; font-size:13px;
    font-variant-numeric:tabular-nums;
  }
  .readout b{font-weight:600; font-size:15px}
  .readout span{color:#5f6672}
  .controls{display:flex; gap:8px; flex-wrap:wrap; align-items:center; font-size:13px}
  button{
    font:inherit; font-size:13px; color:var(--ink); background:transparent;
    border:1px solid var(--hair); border-radius:0; padding:5px 13px; cursor:pointer;
  }
  button:hover{background:#e0e3ea}
  button[aria-pressed="true"]{background:var(--ink); color:var(--paper); border-color:var(--ink)}
  button:focus-visible{outline:2px solid var(--plus); outline-offset:2px}
  .sep{width:1px; height:20px; background:var(--hair); margin:0 4px}
  label{display:inline-flex; align-items:center; gap:7px; color:#4b515c}
  input[type=range]{width:104px; accent-color:var(--ink)}
  .note{font-size:12.5px; line-height:1.8; color:#5f6672; margin:18px 0 0; max-width:62ch}
  .note code{font-family:ui-monospace,monospace; font-size:12px}
  @media (prefers-reduced-motion:reduce){ }
</style>

<div class="wrap">
  <h1>井戸のかたち ── ひし形30面体の上を走る整数のオートマトン</h1>
  <p class="lede">番地は面の辺ごとに一つ。隣から来た量を足して、自分の来た分を引く。これだけを繰り返す。
  壁のある担体（笠）では量が出られず、閉じた担体（30面全体）では戻ってくる。</p>

  <canvas id="cv" width="1560" height="1000"></canvas>

  <div class="readout">
    <div><span>歩数</span> <b id="rStep">0</b></div>
    <div><span>二乗和（毎歩 t² 倍）</span> <b id="rNorm">1</b></div>
    <div><span>量のある番地</span> <b id="rOcc">1</b> / <b id="rTot">10</b></div>
    <div><span>担体</span> <b id="rCar">笠 10面・壁5枚</b></div>
  </div>

  <div class="controls">
    <button id="bCap"  aria-pressed="true">笠（壁あり）</button>
    <button id="bRt"   aria-pressed="false">30面全体</button>
    <button id="bBelt" aria-pressed="false">帯</button>
    <div class="sep"></div>
    <button id="bPlay">再生</button>
    <button id="bStep">一歩</button>
    <button id="bReset">はじめに戻す</button>
    <div class="sep"></div>
    <label>速さ <input id="sSpeed" type="range" min="1" max="12" value="5"></label>
    <label><input id="cSpin" type="checkbox" checked> 回す</label>
  </div>

  <p class="note">番地は面の辺ごとに一つ、30面全体で120本。担体は稜線だけで描き、玉の飛び出す高さがその番地の量。
  色は<span style="color:#3d6fbf">外へ出る</span>（光る）か<span style="color:#5b6f9c">内へ入る</span>（沈んで暗い）かの二つだけ。
  太い稜線が壁（その先に番地が無い辺）。面をクリックするとそこに量を置き直す。
  計算は整数のまま行い、表示のときだけ毎歩 t で割り戻している（割り戻すと二乗和はちょうど1）。</p>
</div>

<script>
const DATA = __DATA__;
const PHI = (1+Math.sqrt(5))/2;

/* ---------- 幾何 ---------- */
const XYZ = DATA.xyz, FACES = DATA.faces;
const CENT = FACES.map(f=>{
  const c=[0,0,0]; f.forEach(i=>{for(let k=0;k<3;k++)c[k]+=XYZ[i][k]/4;}); return c;
});
const NORM = FACES.map((f,fi)=>{
  const a=XYZ[f[0]],b=XYZ[f[1]],c=XYZ[f[2]];
  const u=[b[0]-a[0],b[1]-a[1],b[2]-a[2]], v=[c[0]-a[0],c[1]-a[1],c[2]-a[2]];
  let n=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];
  const L=Math.hypot(...n); n=n.map(x=>x/L);
  const d=n[0]*CENT[fi][0]+n[1]*CENT[fi][1]+n[2]*CENT[fi][2];
  return d<0? n.map(x=>-x): n;
});

/* ---------- 一歩（整数のまま） ---------- */
let mode="cap", C=DATA.carriers[mode];
let state=[], step=0, seedFace=null;
const T_OF = ()=>BigInt(C.t);

function outArcs(f){
  const out=[]; C.arcs.forEach((a,i)=>{ if(a[0]===f) out.push(i); }); return out;
}
function reset(face){
  C=DATA.carriers[mode]; buildArcAt();
  seedFace = (face!==undefined && C.faces.includes(face))? face : C.seed;
  state = new Array(C.arcs.length).fill(0n);
  const o = outArcs(seedFace); if(o.length) state[o[0]] = 1n;
  step = 0; vCur=arcValues(state,0); vPrev=vCur.slice(); mix=1;
  draw(); readout();
}
function advance(){
  const n=C.arcs.length, next=new Array(n).fill(0n);
  for(let i=0;i<n;i++){
    let s=0n;
    for(const [j,c] of C.rows[i]) s += BigInt(c)*state[j];
    next[i]=s;
  }
  vPrev=vCur.slice(); state=next; step++;
  vCur=arcValues(state,step); mix=0;
}

/* 量は有向辺ごとに持つ。t^歩数 で割ると二乗和がちょうど1になるので目盛りは固定 */
function arcValues(st,k){
  const den = T_OF()**BigInt(k);
  const out=new Array(st.length);
  for(let i=0;i<st.length;i++) out[i]=Number(st[i]*1000000n/den)/1000000;
  return out;
}
let vPrev=[], vCur=[], mix=1;

/* ---------- 描画 ---------- */
const cv=document.getElementById("cv"), ctx=cv.getContext("2d");
let ax=-0.35, ay=-2.5918, drag=null, spin=true;
const LIGHT=[-0.35,-0.62,0.70];

function rot(p){
  let [x,y,z]=p;
  let c=Math.cos(ay), s=Math.sin(ay); [x,z]=[x*c+z*s, -x*s+z*c];
  c=Math.cos(ax); s=Math.sin(ax);    [y,z]=[y*c-z*s, y*s+z*c];
  return [x,y,z];
}
function proj(p){
  const S=Math.min(cv.width,cv.height)*0.225;
  return [cv.width/2 + p[0]*S, cv.height/2 - p[1]*S];
}
let hit=[];
function poly(P){
  ctx.beginPath(); ctx.moveTo(P[0][0],P[0][1]);
  for(let k=1;k<P.length;k++) ctx.lineTo(P[k][0],P[k][1]);
  ctx.closePath();
}
/* 稜線の一覧（面の対つき） */
const EDGES=[];
{
  const m=new Map();
  FACES.forEach((f,fi)=>{ for(let k=0;k<4;k++){
    const a=f[k], b=f[(k+1)%4], key=a<b?a+"_"+b:b+"_"+a;
    if(!m.has(key)) m.set(key,{a,b,fs:[]});
    m.get(key).fs.push(fi);
  }});
  for(const e of m.values()) EDGES.push(e);
}
const ARC_AT={}, SITE={};
function buildArcAt(){
  for(const k in ARC_AT) delete ARC_AT[k];
  for(const k in SITE) delete SITE[k];
  C.arcs.forEach((a,i)=>{
    (ARC_AT[a[0]]=ARC_AT[a[0]]||{})[a[2]]=i;
    const f=FACES[a[0]], c=CENT[a[0]], n=NORM[a[0]];
    const m=[0,1,2].map(j=>(XYZ[f[a[2]]][j]+XYZ[f[(a[2]+1)%4]][j])/2);
    SITE[i]={b:[0,1,2].map(j=>c[j]+(m[j]-c[j])*0.56), n};
  });
}
function height(v){ const a=Math.min(1,Math.abs(v)); return Math.sign(v)*Math.pow(a,0.5)*1.30; }
const GLOW=[126,180,255], DIM=[54,86,150];
function draw(){
  ctx.setTransform(1,0,0,1,0,0);
  ctx.fillStyle="#0d1018"; ctx.fillRect(0,0,cv.width,cv.height);
  const inSet=new Set(C.faces);
  const S=Math.min(cv.width,cv.height)*0.225;
  const items=[];
  for(const e of EDGES){
    const A=rot(XYZ[e.a]), B=rot(XYZ[e.b]);
    const act=e.fs.some(f=>inSet.has(f));
    let wall=false;
    for(const fi of e.fs){
      if(!inSet.has(fi)) continue;
      const f=FACES[fi];
      for(let k=0;k<4;k++){
        const p=f[k], q=f[(k+1)%4];
        if(((p===e.a&&q===e.b)||(p===e.b&&q===e.a)) && (ARC_AT[fi]||{})[k]===undefined) wall=true;
      }
    }
    items.push({kind:"edge", A, B, act, wall, z:(A[2]+B[2])/2});
  }
  for(let i=0;i<C.arcs.length;i++){
    const v=vPrev[i]+(vCur[i]-vPrev[i])*mix, h=height(v), s=SITE[i];
    const b=rot(s.b), c=rot([0,1,2].map(j=>s.b[j]+s.n[j]*h));
    items.push({kind:"ball", v, h, b, c, z:c[2]});
  }
  items.sort((a,b)=>a.z-b.z);
  for(const it of items){
    if(it.kind==="edge"){
      const t=(it.z+2)/4;                       /* 奥ほど淡く */
      const A=proj(it.A), B=proj(it.B);
      ctx.beginPath(); ctx.moveTo(A[0],A[1]); ctx.lineTo(B[0],B[1]);
      if(it.wall){ ctx.lineWidth=4.5; ctx.strokeStyle=`rgba(190,214,255,${0.35+0.55*t})`; }
      else if(it.act){ ctx.lineWidth=1.8; ctx.strokeStyle=`rgba(150,178,215,${0.18+0.42*t})`; }
      else { ctx.lineWidth=1.2; ctx.strokeStyle=`rgba(120,138,170,${0.08+0.16*t})`; }
      ctx.stroke();
      continue;
    }
    const a=Math.min(1,Math.abs(it.v));
    if(a<0.012) continue;
    const out=it.h>=0, col=out?GLOW:DIM;
    const B=proj(it.b), Cp=proj(it.c);
    ctx.beginPath(); ctx.moveTo(B[0],B[1]); ctx.lineTo(Cp[0],Cp[1]);
    ctx.lineWidth=1.6; ctx.strokeStyle=`rgba(${col.join(",")},${out?0.45:0.30})`; ctx.stroke();
    const r=S*(0.045+0.115*Math.sqrt(a));
    const g=ctx.createRadialGradient(Cp[0],Cp[1],0,Cp[0],Cp[1],r*2.25);
    if(out){
      g.addColorStop(0,"rgba(238,246,255,.98)");
      g.addColorStop(0.30,`rgba(${col.join(",")},.92)`);
      g.addColorStop(0.62,`rgba(${col.join(",")},.28)`);
    }else{
      g.addColorStop(0,`rgba(120,150,205,.85)`);
      g.addColorStop(0.34,`rgba(${col.join(",")},.65)`);
      g.addColorStop(0.66,`rgba(${col.join(",")},.14)`);
    }
    g.addColorStop(1,"rgba(0,0,0,0)");
    ctx.save();
    if(out) ctx.globalCompositeOperation="lighter";
    ctx.beginPath(); ctx.arc(Cp[0],Cp[1],r*2.25,0,Math.PI*2);
    ctx.fillStyle=g; ctx.fill();
    ctx.restore();
  }
  hit=[];
  for(let fi=0; fi<FACES.length; fi++){
    if(!inSet.has(fi)) continue;
    const nr=rot(NORM[fi]); if(nr[2]<=0.02) continue;
    hit.push({fi, P:FACES[fi].map(i=>proj(rot(XYZ[i])))});
  }
}
function readout(){
  document.getElementById("rStep").textContent=step;
  let s=0n; for(const x of state) s+=x*x;
  document.getElementById("rNorm").textContent =
    (step<=24)? s.toString() : `${C.t*C.t}^${step}`;
  document.getElementById("rOcc").textContent=state.filter(x=>x!==0n).length;
  document.getElementById("rTot").textContent=C.arcs.length;
  const names={cap:"笠 10面・壁5枚", rt:"30面全体・閉じている", belt:"帯 10面・輪"};
  document.getElementById("rCar").textContent=names[mode];
}

/* ---------- 操作 ---------- */
let playing=false, acc=0, last=0;
function loop(t){
  const dt=Math.min(50,t-last); last=t;
  if(spin && !document.hidden){ ay+=dt*0.00016; }
  if(playing){
    acc+=dt;
    const per=1000/ (+document.getElementById("sSpeed").value);
    if(acc>per){ acc=0; advance(); readout(); }
    mix=Math.min(1, acc/per);
  }
  draw();
  requestAnimationFrame(loop);
}
cv.addEventListener("pointerdown",e=>{drag={x:e.clientX,y:e.clientY,moved:false}; cv.setPointerCapture(e.pointerId);});
cv.addEventListener("pointermove",e=>{
  if(!drag) return;
  ay+=(e.clientX-drag.x)*0.008; ax+=(e.clientY-drag.y)*0.008;
  ax=Math.max(-1.4,Math.min(1.4,ax));
  if(Math.hypot(e.clientX-drag.x,e.clientY-drag.y)>3) drag.moved=true;
  drag.x=e.clientX; drag.y=e.clientY;
});
cv.addEventListener("pointerup",e=>{
  if(drag && !drag.moved){
    const r=cv.getBoundingClientRect();
    const x=(e.clientX-r.left)*cv.width/r.width, y=(e.clientY-r.top)*cv.height/r.height;
    for(let i=hit.length-1;i>=0;i--){
      const P=hit[i].P;
      ctx.beginPath(); ctx.moveTo(P[0][0],P[0][1]);
      for(let k=1;k<4;k++) ctx.lineTo(P[k][0],P[k][1]);
      ctx.closePath();
      if(ctx.isPointInPath(x,y)){ reset(hit[i].fi); break; }
    }
  }
  drag=null;
});
function setMode(m){
  mode=m;
  for(const [id,mm] of [["bCap","cap"],["bRt","rt"],["bBelt","belt"]])
    document.getElementById(id).setAttribute("aria-pressed", String(mm===m));
  reset();
}
document.getElementById("bCap").onclick=()=>setMode("cap");
document.getElementById("bRt").onclick=()=>setMode("rt");
document.getElementById("bBelt").onclick=()=>setMode("belt");
document.getElementById("bStep").onclick=()=>{advance();mix=1;readout();draw();};
document.getElementById("bReset").onclick=()=>reset(seedFace);
document.getElementById("bPlay").onclick=e=>{
  playing=!playing; e.target.textContent=playing?"止める":"再生";
};
document.getElementById("cSpin").onchange=e=>{spin=e.target.checked;};
if(window.matchMedia("(prefers-reduced-motion:reduce)").matches){
  spin=false; document.getElementById("cSpin").checked=false;
}
reset();
requestAnimationFrame(t=>{last=t; loop(t);});
</script>
"""

open("well_rt.html", "w").write(HTML.replace("__DATA__", DATA))
print("well_rt.html", len(HTML.replace("__DATA__", DATA)), "bytes")
