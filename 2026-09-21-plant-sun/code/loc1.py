# 観測の深さが住所の精度になるか：半径 r の近傍が同じ点（候補）の数を、ペンローズ担体と四角の格子（周期）で数える
#  近傍＝中心からの距離² ≤ (rφ)² の点の差の集合（Z[ζ] の整数、判定は phi_lt の整数比較）
#  事前の基準  OK：ペンローズは r とともに候補が減る／四角（周期）は全点のまま  NG：ペンローズでも減らない
#  向き：羅針盤あり（そのまま）と、なし（ζ^k で回して重なるものを同一視、四角は 90° 回し）の二通り
#  内側：r_max の近傍が担体の中に収まる点だけ数える（場所の選別にだけ浮動小数）
import math, collections, b13_chain_units as U
exec(open('gen10.py').read().split("pid={}")[0])
cells=U.fits([place(c,k) for (w,r,c,k) in rings]); Q=list(cells); z0=(2,-2,0,-3)
def pol(z): x,y=U.xy(U.zsub(z,z0)); return math.hypot(x,y), math.atan2(y,x)
bins=collections.defaultdict(float)
for q in Q:
    d,a=pol(q); b=int((a+math.pi)/(2*math.pi)*72)%72; bins[b]=max(bins[b],d)
Rfill=min(bins.values()); RMAX=6; PHI=(1+5**.5)/2
inner=[q for q in Q if pol(q)[0]+RMAX*PHI+1e-9<Rfill]
print(f"担体の五角形 {len(Q)}  すき間なく覆う半径 {Rfill:.1f}  内側の点 {len(inner)}")
grid=collections.defaultdict(list)
for q in Q:
    x,y=U.xy(q); grid[(int(x//4),int(y//4))].append(q)
def ball(q,r):
    th=(r*r,r*r)                                   # (rφ)² = r²(1+φ)
    x,y=U.xy(q); R=int(r*PHI//4)+2; out=[]
    for dx in range(-R,R+1):
        for dy in range(-R,R+1):
            for p in grid.get((int(x//4)+dx,int(y//4)+dy),()):
                d=U.zsub(p,q); n=U.norm2(d)
                if not U.phi_lt(th,n): out.append(d)
    return out
def sig(offs,rot):
    if not rot: return tuple(sorted(offs))
    return min(tuple(sorted(U.zrot(d,k) for d in offs)) for k in range(10))
def count(nodes,nb,rot):
    g=collections.Counter(); s={}
    for q in nodes: s[q]=sig(nb(q),rot); g[s[q]]+=1
    c=[g[s[q]] for q in nodes]
    return len(nodes), sum(c)/len(c), max(c), sum(1 for x in c if x==1), min(c)
# 四角の格子（周期）：点数を内側の点数にそろえる
N=int(len(inner)**.5)
def sq_ball(r):
    return [(dx,dy) for dx in range(-r*2,r*2+1) for dy in range(-r*2,r*2+1) if dx*dx+dy*dy<=r*r*3]  # 半径 r·φ 相当（φ²≈2.6→3）
def sq_count(r,rot):
    offs=sq_ball(r); g=collections.Counter(); s={}
    for i in range(N):
        for j in range(N):
            pts=[(dx,dy) for dx,dy in offs]      # 周期の格子：近傍は折り返した先の点で、差は (dx,dy) のまま
            t=tuple(sorted(pts))
            if rot: t=min(tuple(sorted(((a,b) if k==0 else (-b,a) if k==1 else (-a,-b) if k==2 else (b,-a)) for a,b in pts)) for k in range(4))
            s[(i,j)]=t; g[t]+=1
    c=[g[v] for v in s.values()]
    return N*N, sum(c)/len(c), max(c), sum(1 for x in c if x==1), min(c)
print("r   | ペンローズ 羅針盤あり（点数・平均候補・最大候補・一意）| 羅針盤なし | 四角 羅針盤あり | 四角 なし")
for r in range(1,RMAX+1):
    a=count(inner,lambda q:ball(q,r),False); b=count(inner,lambda q:ball(q,r),True)
    c=sq_count(r,False); d=sq_count(r,True)
    f=lambda t:f"{t[1]:8.1f} {t[2]:5d} {t[3]:5d} min{t[4]}"
    print(f"r={r} |{f(a)} |{f(b)} |{f(c)} |{f(d)}", flush=True)
