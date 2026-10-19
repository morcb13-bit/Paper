# 炭素を正四面体の中心に置く（ダイヤモンドの置き方）。A の 3 個・B の 3 個が大きな正四面体をつくるか、二つで星型正八面体になるか
from itertools import combinations
add=lambda a,b:tuple(x+y for x,y in zip(a,b)); sub=lambda a,b:tuple(x-y for x,y in zip(a,b)); neg=lambda a:tuple(-x for x in a)
d2=lambda a,b:sum((x-y)**2 for x,y in zip(a,b))
DIA=[(3,-3,-3),(-3,3,-3),(-3,-3,3),(3,3,3)]
C=[(4,-2,-2),(2,2,-4),(-2,4,-2),(-4,2,2),(-2,-2,4),(2,-4,2)]
CC=[add(p,(1,1,1)) if i%2==0 else p for i,p in enumerate(C)]          # m38 の椅子形
CC=[tuple(2*x for x in p) for p in CC]; CC=[sub(p,(1,1,1)) for p in CC]  # 2 倍して中心を原点へ（中心 (1/2,1/2,1/2)）
D2=[tuple(2*x for x in d) for d in DIA]
TET=[[add(p,d) if i%2==0 else sub(p,d) for d in D2] for i,p in enumerate(CC)]
print("炭素（2 倍）",CC,"高さ",[sum(p) for p in CC])
print("正四面体の辺²",sorted({d2(a,b) for T in TET for a,b in combinations(T,2)}),"中心＝炭素か",all(tuple(sum(x) for x in zip(*T))==tuple(4*c for c in p) for T,p in zip(TET,CC)))
print("隣の炭素は自分の正四面体の頂点か",all(CC[(i+1)%6] in TET[i] and CC[i-1] in TET[i] for i in range(6)))
for s,idx in [("A",(0,2,4)),("B",(1,3,5))]:
    V=[v for i in idx for v in TET[i]]; from collections import Counter; cnt=Counter(V)
    print(s,"頂点の数",len(cnt),"2 個で分け合う頂点",[v for v,k in cnt.items() if k==2])
    # 大きな正四面体の角：どの小さな正四面体にも 1 個だけ属し、ほかの 2 個の頂点と辺² 1152 (=288×4)
    corners=[v for v,k in cnt.items() if k==1]
    big=[a for a in corners if sum(1 for b in corners if d2(a,b)==1152)==3]
    print("  辺² 1152 で結ばれる角",big)
print("----")
from collections import Counter
BIG={}
for s,idx in [("A",(0,2,4)),("B",(1,3,5))]:
    cnt=Counter(v for i in idx for v in TET[i]); U=[v for v,k in cnt.items() if k==1]
    print(s,"1 個だけに属する頂点",U)
    pairs=[(a,b) for a,b in combinations(U,2) if d2(a,b)==1152]
    tri=sorted({v for p in pairs for v in p}); print("  辺² 1152 で結ばれる頂点",tri,"組",len(pairs))
    # 4 つ目の角：3 角から等距離 1152
    if len(tri)==3:
        g=tuple(sum(x) for x in zip(*tri))
        for k in range(-60,61):
            t=add(tuple(x//3 for x in g),(k,k,k)) if all(x%3==0 for x in g) else None
            if t and all(d2(t,c)==1152 for c in tri): print("  4 つ目の角",t); BIG[s]=tri+[t]
for s,V in BIG.items(): print(s,"大きな正四面体の重心×4",tuple(sum(x) for x in zip(*V)))
if len(BIG)==2:
    P=BIG["A"]+BIG["B"]; print("8 頂点が立方体か：互いの距離²",sorted(Counter(d2(a,b) for a,b in combinations(P,2)).items()))
    print("B は A の反転か",sorted(neg(v) for v in BIG["A"])==sorted(BIG["B"]))
    # 4 つ目の小さな正四面体（空いている角）
    for s,idx in [("A",(0,2,4)),("B",(1,3,5))]:
        top=BIG[s][3]; mids=[tuple((a+b)//2 for a,b in zip(top,c)) for c in BIG[s][:3]]
        print(s,"空いた角の小さな正四面体：中点",mids,"が既存の頂点か",[m in {v for i in idx for v in TET[i]} for m in mids])
