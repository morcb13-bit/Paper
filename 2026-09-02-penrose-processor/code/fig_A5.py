# 図A5：黄金整数の半径だけを残した二族の同心円。交点を結ぶと離心率が φ の冪の
#        楕円・双曲線が出る。
import json, math
from collections import defaultdict
import b13_chain_units as U
phi=(1+5**0.5)/2
cells={tuple(int(x) for x in k.split(",")):a for k,a in json.load(open("carrier_1245.json"))["cells"].items()}
POLY=[[U.xy(U.zadd(q,U.zt(a+2*i))) for i in range(5)] for q,a in cells.items()]
SG=[]
for area,cyc in U.gaps(cells):
    if abs(area-2.9389)<0.01:
        c=(0,0,0,0)
        for p in cyc: c=U.zadd(c,p)
        SG.append(((U.xy(c)[0]/10,U.xy(c)[1]/10),[U.xy(p) for p in cyc]))
S=[s for s,_ in SG]
best=None
for a in S:
    if not (10<math.hypot(*a)<20): continue
    for b in S:
        if abs(math.hypot(a[0]-b[0],a[1]-b[1])-phi**4)<1e-6:
            m=math.hypot((a[0]+b[0])/2,(a[1]+b[1])/2)
            if best is None or m<best[0]: best=(m,a,b)
A0,B0=best[1],best[2]; mx,my=(A0[0]+B0[0])/2,(A0[1]+B0[1])/2
th=-math.atan2(B0[1]-A0[1],B0[0]-A0[0]); ct,st=math.cos(th),math.sin(th)
def tr(p):
    x,y=p[0]-mx,p[1]-my; return (x*ct-y*st, x*st+y*ct)
POLY=[[tr(p) for p in q] for q in POLY]; SG=[(tr(s),[tr(p) for p in g]) for s,g in SG]
S=[s for s,_ in SG]; A,B=tr(A0),tr(B0); N=math.hypot(A[0]-B[0],A[1]-B[1])
LAD=[(2,3),(4,7),(5,8),(6,10),(7,11),(9,15),(11,18),(12,19),(16,26)]
val=lambda t:t[0]+t[1]*phi
RD=[val(t) for t in LAD]
# 両方の距離が黄金整数になる五芒星
def isg(v):
    for b in range(-90,91):
        a=v-b*phi
        if abs(a-round(a))<1e-6: return True
    return False
both=[p for p in S if isg(math.hypot(p[0]-A[0],p[1]-A[1])) and isg(math.hypot(p[0]-B[0],p[1]-B[1]))]
V=30.0
o=['<svg xmlns="http://www.w3.org/2000/svg" viewBox="%.1f %.1f %.1f %.1f" width="1000" height="1000">'%(-V,-V,2*V,2*V),
   '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#0e1014"/>'%(-V,-V,2*V,2*V),
   '<clipPath id="c"><rect x="%.1f" y="%.1f" width="%.1f" height="%.1f"/></clipPath>'%(-V,-V,2*V,2*V)]
o.append('<g stroke="#333b46" stroke-width="0.055" fill="none">')
for p in POLY:
    if all(abs(x)<V+2 and abs(y)<V+2 for x,y in p):
        o.append('<polygon points="%s"/>'%(" ".join("%.4f,%.4f"%v for v in p)))
o.append('</g>')
o.append('<g fill="#c8912f" fill-opacity="0.22">')
for s,g in SG:
    if abs(s[0])<V+2 and abs(s[1])<V+2:
        o.append('<polygon points="%s"/>'%(" ".join("%.4f,%.4f"%v for v in g)))
o.append('</g>')
o.append('<g clip-path="url(#c)" fill="none" stroke-width="0.10">')
for r in RD:
    o.append('<circle cx="%.5f" cy="%.5f" r="%.5f" stroke="#e0b050" stroke-opacity="0.55"/>'%(A[0],A[1],r))
    o.append('<circle cx="%.5f" cy="%.5f" r="%.5f" stroke="#6fa8dc" stroke-opacity="0.55"/>'%(B[0],B[1],r))
o.append('</g>')
cx=(A[0]+B[0])/2
o.append('<g clip-path="url(#c)" font-family="Menlo,Consolas,monospace">')
for k,lab in ((7,"e=φ⁻³"),(8,"e=φ⁻⁴"),(9,"e=φ⁻⁵")):
    j=phi**k; a=j/2; c=N/2; b=math.sqrt(a*a-c*c)
    o.append('<ellipse cx="%.6f" cy="0" rx="%.6f" ry="%.6f" fill="none" stroke="#7ee0a0" stroke-width="0.20"/>'%(cx,a,b))
    o.append('<text x="%.3f" y="%.3f" fill="#7ee0a0" font-size="1.3" text-anchor="middle">%s（j=φ^%d）</text>'%(cx,-b-0.6,lab,k))
for k,lab in ((2,"e=φ²"),(3,"e=φ"),(4,"e=1")):
    kk=phi**k; a=kk/2; c=N/2
    if a>c: continue
    b=math.sqrt(c*c-a*a)
    for s in (1,-1):
        pts=[]; t=-3.0
        while t<=3.0:
            x=cx+s*a*math.cosh(t); y=b*math.sinh(t)
            if abs(x)<V and abs(y)<V: pts.append((x,y))
            t+=0.03
        if len(pts)>1:
            o.append('<polyline points="%s" fill="none" stroke="#e08a8a" stroke-width="0.20"/>'%(" ".join("%.4f,%.4f"%p for p in pts)))
    o.append('<text x="%.3f" y="1.7" fill="#e08a8a" font-size="1.3" text-anchor="middle">%s</text>'%(cx+a,lab))
o.append('</g>')
o.append('<g fill="#ffffff">')
for m in RD:
    for n in RD:
        if m+n<N or abs(m-n)>N: continue
        x=(N*N+m*m-n*n)/(2*N); y2=m*m-x*x
        if y2<=0: continue
        y=math.sqrt(y2)
        for s in (1,-1):
            px,py=A[0]+x, s*y
            if abs(px)<V and abs(py)<V: o.append('<circle cx="%.4f" cy="%.4f" r="0.24"/>'%(px,py))
o.append('</g>')
o.append('<g fill="#ff6b6b"><circle cx="%.4f" cy="%.4f" r="0.42"/><circle cx="%.4f" cy="%.4f" r="0.42"/></g>'%(A[0],A[1],B[0],B[1]))
o.append('<g font-family="Menlo,Consolas,monospace" fill="#8fa0b4" font-size="1.35">')
o.append('<text x="%.1f" y="%.1f">五芒星どうしの距離のうち Z[φ] に入るものだけで同心円を張る（半径 %d 本）</text>'%(-V+1.2,-V+2.6,len(RD)))
o.append('<text x="%.1f" y="%.1f">φ⁴ 4+7φ φ⁶ 6+10φ 7+11φ 9+15φ 11+18φ 12+19φ 16+26φ　焦点距離 N=φ⁴</text>'%(-V+1.2,-V+4.5))
o.append('<text x="%.1f" y="%.1f" fill="#7ee0a0">和が φ の冪 → 楕円　e=φ⁴⁻ᵏ</text>'%(-V+1.2,-V+6.4))
o.append('<text x="%.1f" y="%.1f" fill="#e08a8a">差が φ の冪 → 双曲線　e=φ⁴⁻ᵏ　（k=4 は境界 e=1）</text>'%(-V+1.2,-V+8.3))
o.append('</g></svg>')
open("/mnt/user-data/outputs/fig_A5_golden_conics.svg","w").write("\n".join(o))
import os
print("両方の距離が黄金整数になる五芒星 %d 個（この担体の %d 個中）" % (len(both),len(S)))
print("size",os.path.getsize("/mnt/user-data/outputs/fig_A5_golden_conics.svg"))
