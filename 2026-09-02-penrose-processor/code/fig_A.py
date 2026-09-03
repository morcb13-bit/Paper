# 図A：五芒星の殻は r ではなく r² で刻むと乗る
import json, math
from collections import defaultdict
import b13_chain_units as U
phi=(1+5**0.5)/2
cells={tuple(int(x) for x in k.split(",")):a for k,a in json.load(open("carrier_300.json"))["cells"].items()}
TC={0:(2,0),1:(0,1),2:(-1,1),3:(1,-1)}
def norm2(v):
    A=B=0
    for j in range(4):
        for k in range(4):
            a,b=TC[abs(j-k)]; A+=v[j]*v[k]*a; B+=v[j]*v[k]*b
    return A//2,B//2
polys=[]
for q,a in cells.items():
    polys.append([U.xy(U.zadd(q,U.zt(a+2*i))) for i in range(5)])
gaps=U.gaps(cells)
stars=[c for area,c in gaps if abs(area-2.9389)<0.01]
shells={}; SC=[]
for cyc in stars:
    c=(0,0,0,0)
    for p in cyc: c=U.zadd(c,p)
    A,B=norm2(c); r=math.hypot(*U.xy(c))/10
    SC.append((U.xy(c)[0]/10,U.xy(c)[1]/10,[U.xy(p) for p in cyc]))
    shells[round(r,6)]=(A//100,B//100)
S=sorted(shells)
R=58.0
out=[]
w=2*R
out.append('<svg xmlns="http://www.w3.org/2000/svg" viewBox="%.2f %.2f %.2f %.2f" width="900" height="900">'%(-R,-R,w,w))
out.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="#0e1014"/>'%(-R,-R,w,w))
out.append('<g stroke="#2c333d" stroke-width="0.06" fill="none">')
for p in polys:
    out.append('<polygon points="%s"/>'%(" ".join("%.4f,%.4f"%(x,y) for x,y in p)))
out.append('</g>')
# 五芒星の隙間
out.append('<g fill="#c8912f" fill-opacity="0.30" stroke="none">')
for _,_,cyc in SC:
    out.append('<polygon points="%s"/>'%(" ".join("%.4f,%.4f"%(x,y) for x,y in cyc)))
out.append('</g>')
# 殻の円
out.append('<g fill="none" stroke="#e0b050" stroke-width="0.16" stroke-opacity="0.85">')
for r in S:
    if r>0: out.append('<circle cx="0" cy="0" r="%.6f"/>'%r)
out.append('</g>')
# 中心の点
out.append('<g fill="#f2e2b8">')
for x,y,_ in SC: out.append('<circle cx="%.4f" cy="%.4f" r="0.30"/>'%(x,y))
out.append('</g>')
# ラベル（左上へ伸ばす）
out.append('<g font-family="Menlo,Consolas,monospace" fill="#f2e2b8" font-size="1.55">')
for i,r in enumerate(S):
    if r==0: continue
    A,B=shells[r]
    th=math.radians(112+ (0 if i%2 else 6))
    x,y=r*math.cos(th), r*math.sin(th)
    out.append('<text x="%.3f" y="%.3f" text-anchor="middle">r²=%d+%dφ</text>'%(x,-y-0.7,A,B))
out.append('</g>')
out.append('<g font-family="Menlo,Consolas,monospace" fill="#8fa0b4" font-size="2.0">')
out.append('<text x="%.1f" y="%.1f">半径の二乗が黄金整数のところに円を置くと、五芒星の殻に乗る</text>'%(-R+2.5,-R+4.5))
out.append('<text x="%.1f" y="%.1f" fill="#c8912f">13+21φ=φ⁸   47+76φ=φ⁸+φ¹⁰   89+144φ=φ¹²</text>'%(-R+2.5,-R+7.6))
out.append('<text x="%.1f" y="%.1f">%d殻・五芒星%d個・五角形%d枚</text>'%(-R+2.5,R-2.5,len(S)-1,len(SC),len(cells)))
out.append('</g></svg>')
open("/mnt/user-data/outputs/fig_A_shells.svg","w").write("\n".join(out))
print("殻 %d 個" % (len(S)-1))
for r in S:
    if r: print("  r=%9.6f   r²=%d+%dφ   検算 %.6f" % (r,*shells[r],(shells[r][0]+shells[r][1]*phi)))
import os; print("size", os.path.getsize("/mnt/user-data/outputs/fig_A_shells.svg"))
