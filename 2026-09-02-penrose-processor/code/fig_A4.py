# 図A4：中心の異なる二つの同心円族をペンローズに重ねる
#   焦点 A, B は隣り合う五芒星の中心（距離 φ⁴）。各点のまわりに、そこから
#   他の五芒星までの距離（＝√(黄金整数)）で同心円を張る。
#   五芒星はすべて、二つの族の交点に乗る。
import json, math
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
A,B=best[1],best[2]
mx,my=(A[0]+B[0])/2,(A[1]+B[1])/2
th=-math.atan2(B[1]-A[1],B[0]-A[0]); ct,st=math.cos(th),math.sin(th)
def tr(p):
    x,y=p[0]-mx,p[1]-my
    return (x*ct-y*st, x*st+y*ct)
POLY=[[tr(p) for p in q] for q in POLY]
SG=[(tr(s),[tr(p) for p in g]) for s,g in SG]
S=[s for s,_ in SG]
A,B=tr(A),tr(B)
V=24.0
RA=sorted({round(math.hypot(p[0]-A[0],p[1]-A[1]),5) for p in S})
RB=sorted({round(math.hypot(p[0]-B[0],p[1]-B[1]),5) for p in S})
RA=[r for r in RA if 0<r<=V*1.5]; RB=[r for r in RB if 0<r<=V*1.5]
o=['<svg xmlns="http://www.w3.org/2000/svg" viewBox="%.1f %.1f %.1f %.1f" width="950" height="950">'%(-V,-V,2*V,2*V),
   '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#0e1014"/>'%(-V,-V,2*V,2*V),
   '<clipPath id="c"><rect x="%.1f" y="%.1f" width="%.1f" height="%.1f"/></clipPath>'%(-V,-V,2*V,2*V)]
o.append('<g stroke="#39424e" stroke-width="0.055" fill="none">')
for p in POLY:
    if all(abs(x)<V+2 and abs(y)<V+2 for x,y in p):
        o.append('<polygon points="%s"/>'%(" ".join("%.4f,%.4f"%v for v in p)))
o.append('</g>')
o.append('<g fill="#c8912f" fill-opacity="0.26">')
for s,g in SG:
    if abs(s[0])<V+2 and abs(s[1])<V+2:
        o.append('<polygon points="%s"/>'%(" ".join("%.4f,%.4f"%v for v in g)))
o.append('</g>')
o.append('<g clip-path="url(#c)" fill="none" stroke-width="0.09">')
for r in RA: o.append('<circle cx="%.5f" cy="%.5f" r="%.5f" stroke="#e0b050" stroke-opacity="0.60"/>'%(A[0],A[1],r))
for r in RB: o.append('<circle cx="%.5f" cy="%.5f" r="%.5f" stroke="#6fa8dc" stroke-opacity="0.60"/>'%(B[0],B[1],r))
o.append('</g>')
o.append('<g fill="#ffffff">')
n=0
for s in S:
    if abs(s[0])<V and abs(s[1])<V: o.append('<circle cx="%.4f" cy="%.4f" r="0.26"/>'%s); n+=1
o.append('</g>')
o.append('<g fill="#ff6b6b"><circle cx="%.4f" cy="%.4f" r="0.45"/><circle cx="%.4f" cy="%.4f" r="0.45"/></g>'%(A[0],A[1],B[0],B[1]))
o.append('<g font-family="Menlo,Consolas,monospace" fill="#8fa0b4" font-size="1.15">')
o.append('<text x="%.1f" y="%.1f">中心の異なる二つの同心円族。焦点は隣り合う五芒星（赤）、距離 φ⁴</text>'%(-V+1.0,-V+2.4))
o.append('<text x="%.1f" y="%.1f" fill="#e0b050">A のまわり %d 本</text>'%(-V+1.0,-V+4.2,len(RA)))
o.append('<text x="%.1f" y="%.1f" fill="#6fa8dc">B のまわり %d 本　半径はどちらも √(黄金整数)</text>'%(-V+1.0,-V+5.9,len(RB)))
o.append('<text x="%.1f" y="%.1f">この窓の五芒星 %d 個は、すべて二族の交点に乗る</text>'%(-V+1.0,V-1.4,n))
o.append('</g></svg>')
open("/mnt/user-data/outputs/fig_A4_two_centers.svg","w").write("\n".join(o))
# 交点のうち五芒星になるものの割合
hit=tot=0
for ra in RA:
    for rb in RB:
        d=math.hypot(A[0]-B[0],A[1]-B[1])
        if ra+rb<d or abs(ra-rb)>d: continue
        x=(d*d+ra*ra-rb*rb)/(2*d); y2=ra*ra-x*x
        if y2<=0: continue
        y=math.sqrt(y2)
        for sg in (1,-1):
            px,py=A[0]+x*(B[0]-A[0])/d - sg*y*(B[1]-A[1])/d, A[1]+x*(B[1]-A[1])/d + sg*y*(B[0]-A[0])/d
            if abs(px)<V and abs(py)<V:
                tot+=1
                if any(abs(px-s[0])<1e-4 and abs(py-s[1])<1e-4 for s in S): hit+=1
import os
print("size",os.path.getsize("/mnt/user-data/outputs/fig_A4_two_centers.svg"))
print("A のまわり %d 本 / B のまわり %d 本 / 窓の五芒星 %d 個" % (len(RA),len(RB),n))
print("窓の中の交点 %d 個、うち五芒星 %d 個（%.1f%%）" % (tot,hit,100*hit/tot))
