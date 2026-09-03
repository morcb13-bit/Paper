# 図A2：二つの五芒星を焦点にした黄金円。交点が楕円と双曲線を作り、離心率が φ の冪になる
import json, math
import b13_chain_units as U
phi=(1+5**0.5)/2
cells={tuple(int(x) for x in k.split(",")):a for k,a in json.load(open("carrier_300.json"))["cells"].items()}
polys=[[U.xy(U.zadd(q,U.zt(a+2*i))) for i in range(5)] for q,a in cells.items()]
stars=[]
for area,cyc in U.gaps(cells):
    if abs(area-2.9389)<0.01:
        c=(0,0,0,0)
        for p in cyc: c=U.zadd(c,p)
        stars.append((U.xy(c)[0]/10,U.xy(c)[1]/10))
N=phi**4
tgt=min(stars,key=lambda s:abs(math.hypot(*s)-N)+ (0 if s[1]>=0 else 0))
th=-math.atan2(tgt[1],tgt[0]); ct,st=math.cos(th),math.sin(th)
rot=lambda p:(p[0]*ct-p[1]*st, p[0]*st+p[1]*ct)
polys=[[rot(p) for p in poly] for poly in polys]
stars=[rot(s) for s in stars]
A=(0.0,0.0); B=(N,0.0)
X0,X1,Y0,Y1=-15.5,22.5,-16.0,16.0
o=[]
o.append('<svg xmlns="http://www.w3.org/2000/svg" viewBox="%.2f %.2f %.2f %.2f" width="1000" height="842">'%(X0,Y0,X1-X0,Y1-Y0))
o.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="#0e1014"/>'%(X0,Y0,X1-X0,Y1-Y0))
o.append('<g stroke="#242a33" stroke-width="0.045" fill="none">')
for p in polys:
    if all(X0-2<x<X1+2 and Y0-2<y<Y1+2 for x,y in p):
        o.append('<polygon points="%s"/>'%(" ".join("%.4f,%.4f"%v for v in p)))
o.append('</g>')
o.append('<g fill="#c8912f" fill-opacity="0.16"><circle cx="0" cy="0" r="0.9"/><circle cx="%.4f" cy="0" r="0.9"/></g>'%N)
# 黄金円
RD=[phi**k for k in range(0,6)]
o.append('<g fill="none" stroke-width="0.075">')
for r in RD:
    o.append('<circle cx="0" cy="0" r="%.6f" stroke="#e0b050" stroke-opacity="0.55"/>'%r)
    o.append('<circle cx="%.6f" cy="0" r="%.6f" stroke="#6fa8dc" stroke-opacity="0.55"/>'%(N,r))
o.append('</g>')
# 円錐曲線
def ell(j,col,lab):
    a=j/2; c=N/2; b=math.sqrt(max(a*a-c*c,0))
    o.append('<ellipse cx="%.6f" cy="0" rx="%.6f" ry="%.6f" fill="none" stroke="%s" stroke-width="0.16"/>'%(c,a,b,col))
    o.append('<text x="%.3f" y="%.3f" fill="%s" font-size="1.25" text-anchor="middle">%s</text>'%(c,-b-0.55,col,lab))
def hyp(k,col,lab):
    a=k/2; c=N/2; b=math.sqrt(max(c*c-a*a,0)); cx=N/2
    for s in (1,-1):
        pts=[]
        t=-2.05
        while t<=2.05:
            pts.append((cx+s*a*math.cosh(t), b*math.sinh(t))); t+=0.05
        pts=[(x,y) for x,y in pts if X0<x<X1 and Y0<y<Y1]
        if len(pts)>1:
            o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="0.16"/>'%(" ".join("%.4f,%.4f"%p for p in pts),col))
    o.append('<text x="%.3f" y="%.3f" fill="%s" font-size="1.25" text-anchor="middle">%s</text>'%(cx+a,1.55,col,lab))
o.append('<g font-family="Menlo,Consolas,monospace">')
ell(phi**5,"#7ee0a0","e=φ⁻¹")
ell(phi**6,"#7ee0a0","e=φ⁻²")
ell(phi**7,"#7ee0a0","e=φ⁻³")
hyp(phi**1,"#e08a8a","e=φ³")
hyp(phi**2,"#e08a8a","e=φ²")
hyp(phi**3,"#e08a8a","e=φ")
o.append('<line x1="0" y1="0" x2="%.6f" y2="0" stroke="#f0f0f0" stroke-width="0.20" stroke-dasharray="0.5 0.35"/>'%N)
o.append('<text x="%.3f" y="-0.6" fill="#f0f0f0" font-size="1.25" text-anchor="middle">境界 e=1（和 φ²+φ³=φ⁴）</text>'%(N/2))
o.append('</g>')
# 交点
o.append('<g fill="#ffffff">')
seen=set()
for m in RD:
    for n in RD:
        d=N
        if m+n<d or abs(m-n)>d: continue
        x=(d*d+m*m-n*n)/(2*d); y2=m*m-x*x
        if y2<=0: continue
        y=math.sqrt(y2)
        for s in (1,-1):
            if X0<x<X1 and Y0<s*y<Y1: o.append('<circle cx="%.4f" cy="%.4f" r="0.16"/>'%(x,s*y))
o.append('</g>')
o.append('<g font-family="Menlo,Consolas,monospace" fill="#8fa0b4" font-size="1.5">')
o.append('<text x="%.1f" y="%.1f">焦点＝二つの五芒星の中心　距離 N=φ⁴=2+3φ</text>'%(X0+0.8,Y0+2.4))
o.append('<text x="%.1f" y="%.1f" fill="#e0b050">半径 1, φ, φ², φ³, φ⁴ の黄金円を両方の焦点に</text>'%(X0+0.8,Y0+4.4))
o.append('<text x="%.1f" y="%.1f" fill="#7ee0a0">隣り合う円の和 φᵏ+φᵏ⁺¹=φᵏ⁺² → 楕円　e=φ²⁻ᵏ</text>'%(X0+0.8,Y0+6.4))
o.append('<text x="%.1f" y="%.1f" fill="#e08a8a">一つ飛ばした円の差 φᵏ⁺²−φᵏ=φᵏ⁺¹ → 双曲線　e=φ³⁻ᵏ</text>'%(X0+0.8,Y0+8.4))
o.append('</g></svg>')
open("/mnt/user-data/outputs/fig_A2_two_foci.svg","w").write("\n".join(o))
import os; print("size",os.path.getsize("/mnt/user-data/outputs/fig_A2_two_foci.svg"))
print("焦点 B = (%.6f, 0)  ＝ 第1殻の五芒星" % N)
