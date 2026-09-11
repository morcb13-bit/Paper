#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1〜13 の段を、ひし形30面体の鎖で組む。

  奇数の段   折れ角 108度のジグザグ（n個の鎖）
  偶数の段   折れ角 180度の直線（n個の鎖）

どの段も、各長さで**面外の広がりが最小**になるものを選んである（全部 0＝平面に収まる）。
重なりの判定は30本の面法線と支持関数で、通ったものだけを使っている。
"""

import json
import math

import rt_units as ru
from make_rt3d import faces_of_rt

RT = faces_of_rt()
ch = json.load(open("chains_fix.json"))

HUE = {}
for n in range(1, 14):
    if n == 1:
        HUE[n] = (210, 6)           # 1 は灰
    elif n in (5, 7, 11, 13):
        HUE[n] = (34, 46)           # 素数
    elif n % 2 == 0:
        HUE[n] = (212, 26)          # 2 の倍数
    else:
        HUE[n] = (142, 30)          # 3 の倍数

rows = []
gap = 0.6
y = 0.0
allpts = []
for n in range(1, 14):
    pts = [tuple(p) for p in ch[str(n)]]
    # 主軸を x に、面の法線を z に向ける
    g = [sum(p[i] for p in pts) / len(pts) for i in range(3)]
    Q = [[p[i] - g[i] for i in range(3)] for p in pts]
    if len(pts) >= 2:
        ax = ru.unit(tuple(Q[-1][i] - Q[0][i] for i in range(3)))
    else:
        ax = (1.0, 0.0, 0.0)
    nz = None
    if len(pts) >= 3:
        for q in Q[1:]:
            cr = ru.cross(ax, q)
            if math.dist((0.0, 0.0, 0.0), cr) > 1e-9:
                nz = ru.unit(cr)
                break
    if nz is None:
        tmp = (0.0, 0.0, 1.0) if abs(ax[2]) < 0.9 else (1.0, 0.0, 0.0)
        nz = ru.unit(ru.cross(ax, tmp))
    ny = ru.unit(ru.cross(nz, ax))
    R = [ax, ny, nz]
    P = [[sum(R[a][i] * q[i] for i in range(3)) for a in range(3)] for q in Q]
    rows.append({"n": n, "pts": P, "hue": HUE[n][0], "sat": HUE[n][1]})

# 段を縦に積む
H = 3.4
for k, r in enumerate(rows):
    dy = -k * H
    r["pts"] = [[p[0], p[1] + dy, p[2]] for p in r["pts"]]

centers, hues = [], []
for r in rows:
    for p in r["pts"]:
        centers.append(p)
        hues.append([r["hue"], r["sat"]])
rad = max(math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2) for p in centers) + 1.62

data = {"rt": RT, "centers": centers, "hues": hues, "r": rad,
        "rows": [{"n": r["n"], "k": len(r["pts"])} for r in rows]}

HTML = """<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>1〜13 を立体で積む</title>
<style>
  :root{--paper:#f7f4ee;--ink:#2b3138;--sub:#7c827d}
  *{box-sizing:border-box}
  body{margin:0;background:var(--paper);color:var(--ink);
    font-family:"Hiragino Mincho ProN","Yu Mincho",serif;
    display:flex;justify-content:center;padding:16px 12px 28px}
  .wrap{width:100%;max-width:820px}
  h1{font-size:19px;font-weight:400;letter-spacing:.04em;margin:0 0 3px}
  .sub{font-size:13px;color:var(--sub);margin:0 0 14px;line-height:1.75}
  canvas{width:100%;height:auto;display:block;background:#fffdf8;border:1px solid #ddd7cc;
    cursor:grab;touch-action:none}
  canvas:active{cursor:grabbing}
  .row{display:flex;flex-wrap:wrap;gap:14px;margin-top:12px;font-size:13px;color:var(--sub)}
  .k{display:flex;align-items:center;gap:6px}
  .sw{width:13px;height:13px;display:inline-block;border:1px solid #bbb}
  .hint{font-size:12px;color:#9a958a;margin-top:10px;line-height:1.8}
</style>
</head>
<body>
<div class="wrap">
  <h1>1〜13 を立体で積む</h1>
  <p class="sub">上から 1・2・3 …・13。奇数の段は折れ角108度で左右きっちり交互のジグザグ、偶数の段は折れ角180度の直線。ドラッグで回り、ホイールで寄れる。</p>
  <canvas id="cv" width="1100" height="1500"></canvas>
  <div class="row">
    <span class="k"><i class="sw" style="background:hsl(34 46% 72%)"></i>素数 5,7,11,13</span>
    <span class="k"><i class="sw" style="background:hsl(212 26% 72%)"></i>2 の倍数</span>
    <span class="k"><i class="sw" style="background:hsl(142 30% 72%)"></i>3 の倍数</span>
    <span class="k"><i class="sw" style="background:hsl(210 6% 72%)"></i>1</span>
  </div>
  <p class="hint">奇数は LRLRL… と向きが毎回入れ替わり、偶数は一直線。13段まで食い込みゼロ、面外の広がりもゼロ。回すと、平面で見ていた並びが厚みを持っていないことが見える。</p>
</div>
<script>
const D=__DATA__;
const cv=document.getElementById('cv'),g=cv.getContext('2d');
let rx=0, ry=0, zoom=1, drag=null;
function rot(p){
  let [x,y,z]=p;
  let c=Math.cos(ry),s=Math.sin(ry); let X=x*c-z*s, Z=x*s+z*c;
  c=Math.cos(rx); s=Math.sin(rx); let Y=y*c-Z*s; Z=y*s+Z*c;
  return [X,Y,Z];
}
function draw(){
  const S=Math.min(cv.width,cv.height)/(2.1*D.r)*zoom;
  g.fillStyle='#fffdf8'; g.fillRect(0,0,cv.width,cv.height);
  const L=(()=>{const v=[0.42,0.78,0.47];const m=Math.hypot(...v);return v.map(x=>x/m);})();
  const polys=[];
  for(let ci=0;ci<D.centers.length;ci++){
    const c=D.centers[ci], hu=D.hues[ci];
    for(const f of D.rt){
      const n=rot(f.n);
      if(n[2]<0) continue;
      const q=f.q.map(v=>rot([v[0]+c[0],v[1]+c[1],v[2]+c[2]]));
      polys.push({q,n,hu,z:(q[0][2]+q[1][2]+q[2][2]+q[3][2])/4});
    }
  }
  polys.sort((a,b)=>a.z-b.z);
  for(const p of polys){
    const d=p.n[0]*L[0]+p.n[1]*L[1]+p.n[2]*L[2];
    const t=Math.max(0,Math.min(1,(d+1)/2));
    g.beginPath();
    for(let i=0;i<4;i++){
      const X=cv.width/2+p.q[i][0]*S, Y=cv.height/2-p.q[i][1]*S;
      i?g.lineTo(X,Y):g.moveTo(X,Y);
    }
    g.closePath();
    g.fillStyle=`hsl(${p.hu[0]} ${p.hu[1]}% ${58+28*t}%)`;
    g.fill();
    g.strokeStyle='#5b5b6a'; g.lineWidth=0.7; g.stroke();
  }
}
cv.addEventListener('pointerdown',e=>{drag=[e.clientX,e.clientY];cv.setPointerCapture(e.pointerId);});
cv.addEventListener('pointermove',e=>{if(!drag)return;
  ry+=(e.clientX-drag[0])*0.008; rx+=(e.clientY-drag[1])*0.008;
  drag=[e.clientX,e.clientY]; draw();});
cv.addEventListener('pointerup',()=>drag=null);
cv.addEventListener('wheel',e=>{e.preventDefault();zoom*=e.deltaY<0?1.08:1/1.08;
  zoom=Math.max(.4,Math.min(5,zoom));draw();},{passive:false});
draw();
</script>
</body>
</html>"""

open("rt_tower_fix_3d.html", "w").write(HTML.replace("__DATA__", json.dumps(data, separators=(",", ":"))))
print("rt_tower_fix_3d.html  %.1f KB ／ RT %d 個"
      % (len(open("rt_tower_fix_3d.html").read()) / 1024, len(centers)))
