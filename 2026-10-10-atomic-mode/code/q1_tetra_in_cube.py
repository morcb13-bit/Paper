# 立方体＋各面の対角線の図の中にできる四面体を数える（整数のみ）
from itertools import combinations
from collections import Counter
V=[(x,y,z) for x in(-2,2) for y in(-2,2) for z in(-2,2)]
F=[(2,0,0),(-2,0,0),(0,2,0),(0,-2,0),(0,0,2),(0,0,-2)]
P=V+F
sub=lambda a,b:tuple(i-j for i,j in zip(a,b))
d2=lambda a,b:sum(t*t for t in sub(a,b))
cross=lambda a,b:(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
dot=lambda a,b:sum(i*j for i,j in zip(a,b))
# 引いた線分：立方体の辺12本（長さ²16）と面の対角線12本（長さ²32）
segs=[(a,b) for a,b in combinations(V,2) if d2(a,b) in(16,32)]
def on(p,a,b):  # p が線分 ab 上（端含む）
    ab,ap=sub(b,a),sub(p,a)
    return cross(ab,ap)==(0,0,0) and 0<=dot(ap,ab)<=dot(ab,ab)
def joined(p,q): return any(on(p,a,b) and on(q,a,b) for a,b in segs)
tets=[]
for t in combinations(P,4):
    if not all(joined(p,q) for p,q in combinations(t,2)):continue
    a,b,c,d=t
    if dot(sub(b,a),cross(sub(c,a),sub(d,a)))==0:continue
    tets.append(t)
kind=Counter()
for t in tets:
    L=tuple(sorted(d2(p,q) for p,q in combinations(t,2)))
    nv=sum(p in V for p in t)
    kind[(L,nv)]+=1
print('四面体の総数',len(tets))
for (L,nv),n in sorted(kind.items()):print(n,'個  辺²',L,' 立方体の頂点',nv,'・面の中心',4-nv)
