# 図A3：二つの五芒星を焦点にした二極の格子をペンローズに重ねる
#   j = r1²+r2²（中点まわりの円）  k = r1²−r2²（軸に垂直な直線）  どちらも黄金整数
#   すべての五芒星が「円 × 直線」の交点に乗る
import json, math
import b13_chain_units as U
phi=(1+5**0.5)/2
cells={tuple(int(x) for x in k.split(",")):a for k,a in json.load(open("carrier_300.json"))["cells"].items()}
polys=[[U.xy(U.zadd(q,U.zt(a+2*i))) for i in range(5)] for q,a in cells.items()]
gaps=[cyc for area,cyc in U.gaps(cells) if abs(area-2.9389)<0.01]
SC=[]
for cyc in gaps:
    c=(0,0,0,0)
    for p in cyc: c=U.zadd(c,p)
    SC.append(((U.xy(c)[0]/10,U.xy(c)[1]/10),[U.xy(p) for p in cyc]))
N=phi**4
B=min((s for s,_ in SC),key=lambda s:abs(math.hypot(*s)-N))
th=-math.atan2(B[1],B[0]); ct,st=math.cos(th),math.sin(th)
rot=lambda p:(p[0]*ct-p[1]*st, p[0]*st+p[1]*ct)
polys=[[rot(p) for p in q] for q in polys]
SC=[(rot(s),[rot(p) for p in g]) for s,g in SC]
R=52.0; X0,Y0,W=-R,-R,2*R
o=['<svg xmlns="http://www.w3.org/2000/svg" viewBox="%.1f %.1f %.1f %.1f" width="950" height="950">'%(X0,Y0,W,W),
   '<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#0e1014"/>'%(X0,Y0,W,W)]
o.append('<g stroke="#39424e" stroke-width="0.07" fill="none">')
for p in polys: o.append('<polygon points="%s"/>'%(" ".join("%.4f,%.4f"%v for v in p)))
o.append('</g>')
o.append('<g fill="#c8912f" fill-opacity="0.28">')
for _,g in SC: o.append('<polygon points="%s"/>'%(" ".join("%.4f,%.4f"%v for v in g)))
o.append('</g>')
J=set(); K=set()
for (x,y),_ in SC:
    r1=x*x+y*y; r2=(x-N)**2+y*y
    J.add(round(r1+r2,6)); K.add(round(r1-r2,6))
cx=N/2
o.append('<clipPath id="c"><rect x="%.1f" y="%.1f" width="%.1f" height="%.1f"/></clipPath>'%(X0,Y0,W,W))
o.append('<g clip-path="url(#c)" fill="none" stroke="#6fa8dc" stroke-width="0.10" stroke-opacity="0.55">')
for j in sorted(J):
    rr=(j-N*N/2)/2
    if rr>1e-9: o.append('<circle cx="%.6f" cy="0" r="%.6f"/>'%(cx,math.sqrt(rr)))
o.append('</g>')
o.append('<g clip-path="url(#c)" stroke="#7ee0a0" stroke-width="0.10" stroke-opacity="0.55">')
for k in sorted(K):
    x=(k+N*N)/(2*N)
    o.append('<line x1="%.6f" y1="%.1f" x2="%.6f" y2="%.1f"/>'%(x,Y0,x,Y0+W))
o.append('</g>')
o.append('<g fill="#ffffff">')
for (x,y),_ in SC: o.append('<circle cx="%.4f" cy="%.4f" r="0.34"/>'%(x,y))
o.append('</g>')
o.append('<g fill="#ff6b6b"><circle cx="0" cy="0" r="0.62"/><circle cx="%.4f" cy="0" r="0.62"/></g>'%N)
o.append('<g font-family="Menlo,Consolas,monospace" fill="#8fa0b4" font-size="1.9">')
o.append('<text x="%.1f" y="%.1f">焦点＝中心の五芒星と第1殻の五芒星（赤）　距離 N=φ⁴</text>'%(X0+1.5,Y0+4.0))
o.append('<text x="%.1f" y="%.1f" fill="#6fa8dc">円　j = r₁²+r₂²（中点まわり）　%d 本</text>'%(X0+1.5,Y0+7.0,len(J)))
o.append('<text x="%.1f" y="%.1f" fill="#7ee0a0">直線 k = r₁²−r₂²（軸に垂直）　%d 本</text>'%(X0+1.5,Y0+10.0,len(K)))
o.append('<text x="%.1f" y="%.1f">j も k も黄金整数。五芒星 %d 個すべてが交点に乗る</text>'%(X0+1.5,Y0+13.0,len(SC)))
o.append('</g></svg>')
open("/mnt/user-data/outputs/fig_A3_bipolar.svg","w").write("\n".join(o))
import os; print("size",os.path.getsize("/mnt/user-data/outputs/fig_A3_bipolar.svg"),"円",len(J),"直線",len(K),"星",len(SC))
