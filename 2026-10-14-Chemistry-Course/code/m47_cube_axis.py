# への字 3 つの六角形と、トリプチセンの 3 分岐：星型正八面体の 8 個の中心（立方体の角）で数える。整数だけ
from itertools import product, combinations
dot=lambda a,b:sum(x*y for x,y in zip(a,b)); sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
CORN=[tuple(3*x for x in s) for s in product((1,-1),repeat=3)]
edge=lambda p,q:dot(sub(p,q),sub(p,q))==36          # 正四面体どうしが辺を分け合う＝立方体の辺
RING=[(3,3,-3),(3,-3,-3),(3,-3,3),(-3,-3,3),(-3,3,3),(-3,3,-3)]
print("ベンゼンの高さ（(1,1,1) の向き）",[sum(c) for c in RING],"→ への字 3 つ")
P,Q=RING[0],RING[3]; d=sub(Q,P)
print("向かい合う炭素",P,Q,"高さ",sum(P),sum(Q),"（への字の頂点と底）／ P–Q は立方体の対角線か",dot(d,d)==108)
NP=[c for c in CORN if edge(c,P)]; NQ=[c for c in CORN if edge(c,Q)]
print("P の隣 3 つ",NP,"Q の隣 3 つ",NQ)
w=[tuple(dot(d,d)*x-dot(sub(n,P),d)*y for x,y in zip(sub(n,P),d)) for n in NP]   # 軸に垂直な成分（整数倍）
print("P から出る 3 本を軸に沿って見た内積×2 ＋ 長さ²",[2*dot(a,b)+dot(a,a) for a,b in combinations(w,2)],"（0 なら 120 度）")
paths=[(n,m) for n in NP for m in NQ if edge(n,m)]
mt=[S for S in combinations(paths,3) if len({p[0] for p in S})==3 and len({p[1] for p in S})==3]
print("P→隣→隣→Q の道",len(paths),"／ 重ならない 3 本の選び方",len(mt))
for S in mt: print("  3 本",S)
caps=[c for c in CORN if c not in RING]
print("ベンゼンのふた",caps,"は、P の隣と Q の隣に 1 つずつ入っている",[c in NP for c in caps],[c in NQ for c in caps])
