# 一辺3の網（N・O殻）の20点を球へ膨らませると、正十二面体の20頂点に一対一で着くか
# 網の点と内積は整数。正十二面体の照合（最も近い頂点を選ぶ）だけ浮動小数
import math, itertools
V=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
add=lambda *vs:tuple(sum(x) for x in zip(*vs))
mul=lambda k,v:tuple(k*x for x in v)
dot=lambda u,v:sum(a*b for a,b in zip(u,v))
def net(n):
    pts={}
    for w in itertools.product(range(n+1),repeat=4):
        if sum(w)==n and 0 in w:
            p=add(*[mul(w[i],V[i]) for i in range(4)])
            kind=("頂点" if max(w)==n else "面の中" if w.count(0)==1 else "辺")
            pts[p]=kind
    return pts
for n in (1,2,3):
    P=net(n); kinds={}
    for p,k in P.items(): kinds.setdefault(k,set()).add(dot(p,p))
    print(f"一辺{n}: 点{len(P)}  種類と長さ²:",{k:(sum(1 for q in P.values() if q==k),sorted(v)) for k,v in kinds.items()})
# 一辺2の辺の点 = 芯の正八面体の頂点か
print("一辺2の辺の点:",sorted(p for p,k in net(2).items() if k=="辺"))
# 正十二面体（二つの向き）
f=(1+5**.5)/2
def dodeca(o):
    D=[(a,b,c) for a in(1,-1) for b in(1,-1) for c in(1,-1)]
    for s1 in(1,-1):
        for s2 in(1,-1):
            base=[(0,s1/f,s2*f),(s2*f,0,s1/f),(s1/f,s2*f,0)] if o==0 else [(0,s1*f,s2/f),(s2/f,0,s1*f),(s1*f,s2/f,0)]
            D+=base
    return D
unit=lambda v:tuple(x/math.sqrt(dot(v,v)) for x in v)
P=net(3)
for o in (0,1):
    D=[unit(d) for d in dodeca(o)]
    m={}; worst=0
    for p,k in P.items():
        u=unit(p); j=max(range(20),key=lambda i:dot(u,D[i])); m[p]=j
        worst=max(worst,math.degrees(math.acos(min(1,dot(u,D[j])))))
    ok=len(set(m.values()))==20
    by={}
    for p,k in P.items():
        a=math.degrees(math.acos(min(1,dot(unit(p),D[m[p]])))); by.setdefault(k,set()).add(round(a,2))
    print(f"向き{o}: 一対一={ok}  ずれ（度）:",by)
# 網（膨らむ前）の辺の点どうしの組：和が対角線にそろう三本
E=[p for p,k in P.items() if k=="辺"]
tri=set()
for a,b,c in itertools.combinations(E,3):
    s=add(a,b,c)
    if abs(s[0])==abs(s[1])==abs(s[2])!=0 and dot(a,b)==dot(b,c)==dot(c,a):
        tri.add((dot(a,b),dot(a,a),s))
print("辺の点で作る対称な三本（内積,長さ²,和）:",sorted(tri))
