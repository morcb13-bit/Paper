# 二酸化ケイ素（教科書の図の形）：ダイヤモンドの Si–Si 結合の真ん中に O。座標は 2 倍して整数にする
from itertools import combinations
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
DA=[(6,6,6),(6,-6,-6),(-6,6,-6),(-6,-6,6)]; DB=[tuple(-x for x in d) for d in DA]
half=lambda d:tuple(x//2 for x in d)
Si={(0,0,0):("A",0)}; O={}; front=[(0,0,0)]; shells=[1]; oshell=[0]
for k in range(1,7):
    nxt=[];no=0
    for p in front:
        s=Si[p][0]
        for d in (DA if s=="A" else DB):
            o=tuple(p[i]+d[i]//2 for i in range(3)); q=tuple(p[i]+d[i] for i in range(3))
            if o not in O: O[o]=k; no+=1
            if q not in Si: Si[q]=("B" if s=="A" else "A",k); nxt.append(q)
    shells.append(len(nxt)); oshell.append(no); front=nxt
print("Si の世代",shells)
cum=[]
for k in range(5):
    Os=set(tuple(p[i]+d[i]//2 for i in range(3)) for p in Si if Si[p][1]<=k for d in (DA if Si[p][0]=="A" else DB)); cum.append(len(Os))
print("世代 k までの Si の正四面体にある O（累計）",cum)
T=[half(d) for d in DA]
print("SiO4 の正四面体：頂点",T,"辺の長さ²",sorted(set(dot(tuple(a[i]-b[i] for i in range(3)),tuple(a[i]-b[i] for i in range(3))) for a,b in combinations(T,2))),"Si–O 長さ²",dot(T[0],T[0]),"O–Si–O 内積",dot(T[0],T[1]))
# 各 O は何個の正四面体に属するか（内側）
inner=[p for p in Si if Si[p][1]<=4]
cnt={}
for p in Si:
    for d in (DA if Si[p][0]=="A" else DB):
        o=tuple(p[i]+d[i]//2 for i in range(3)); cnt[o]=cnt.get(o,0)+1
inO=[o for o in O if O[o]<=4]
print("内側の O が属する正四面体の数",set(cnt[o] for o in inO))
# 同じ向きどうしが接するか
Apos=[p for p in Si if Si[p][0]=="A"]; Bpos=[p for p in Si if Si[p][0]=="B"]
verts=lambda p:set(tuple(p[i]+d[i]//2 for i in range(3)) for d in (DA if Si[p][0]=="A" else DB))
aa=sum(1 for a,b in combinations(Apos[:200],2) if verts(a)&verts(b)); print("A どうしで頂点を共有する組",aa)
# 単位格子（辺 24）あたり
inside=lambda p:all(0<=x<24 for x in p)
# 周期的に作る：格子 24 の面心立方の点＋(6,6,6) 型
from itertools import product
S=set()
for c in product(range(-1,2),repeat=3):
    for f in [(0,0,0),(0,12,12),(12,0,12),(12,12,0)]:
        p=tuple(24*c[i]+f[i] for i in range(3)); S.add(p); S.add(tuple(p[i]+6 for i in range(3)))
Os=set(tuple((a[i]+b[i])//2 for i in range(3)) for a in S for b in S if dot(tuple(a[i]-b[i] for i in range(3)),tuple(a[i]-b[i] for i in range(3)))==108)
print("単位格子あたり Si",sum(inside(p) for p in S),"O",sum(inside(o) for o in Os))
