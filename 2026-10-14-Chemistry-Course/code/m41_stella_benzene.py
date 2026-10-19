# ベンゼン＝星型正八面体？ 炭素は正四面体の中心、共有結合は辺を分け合う（v288 の決め）。整数だけ
from itertools import combinations, product
from collections import Counter
add=lambda a,b:tuple(x+y for x,y in zip(a,b)); neg=lambda a:tuple(-x for x in a)
d2=lambda a,b:sum((x-y)**2 for x,y in zip(a,b))
E=lambda T:{frozenset(p) for p in combinations(T,2)}
# 星型正八面体：辺 12 の立方体（v288 の枠）。正八面体の頂点は面の中心 (±6,0,0) など、各面に正四面体を立てる
OCT=[tuple(6*s if i==k else 0 for i in range(3)) for k in range(3) for s in (1,-1)]
SG=list(product((1,-1),repeat=3))
TET={s:[tuple(6*s[k] if i==k else 0 for i in range(3)) for k in range(3)]+[tuple(6*x for x in s)] for s in SG}
cen={s:tuple(sum(v[i] for v in T)//4 for i in range(3)) for s,T in TET.items()}
print("小さな正四面体の辺²",sorted({d2(a,b) for T in TET.values() for a,b in combinations(T,2)}),"正八面体の辺²",sorted({d2(a,b) for a,b in combinations(OCT,2) if d2(a,b)<144}))
par=lambda s:"A" if s[0]*s[1]*s[2]>0 else "B"
for p in "AB":
    big=[tuple(6*x for x in s) for s in SG if par(s)==p]; print(p,"大きな正四面体",big,"辺²",{d2(a,b) for a,b in combinations(big,2)},"重心×4",tuple(sum(x) for x in zip(*big)))
# 軸 (1,1,1)：ふた＝ (+,+,+) と (−,−,−)、残る 6 個が帯
CAP=[(1,1,1),(-1,-1,-1)]; BELT=[s for s in SG if s not in CAP]
# 帯を軸のまわりに並べる（辺を分け合う順）
ring=[BELT[0]]
while len(ring)<6:
    for s in BELT:
        if s not in ring and len(E(TET[ring[-1]])&E(TET[s]))==1: ring.append(s); break
print("帯の 6 個（辺を分け合う順に一周）",ring,"閉じるか",len(E(TET[ring[-1]])&E(TET[ring[0]]))==1)
C=[cen[s] for s in ring]
print("炭素＝正四面体の中心",C,"種類",[par(s) for s in ring],"高さ",[sum(c) for c in C])
print("炭素どうしの距離²",[d2(C[i],C[(i+1)%6]) for i in range(6)],"軸に沿って見た距離²（3 倍）",[3*d2(C[i],C[(i+1)%6])-(sum(C[i])-sum(C[(i+1)%6]))**2 for i in range(6)])
for i,s in enumerate(ring):
    sh=[(t,len(E(TET[s])&E(TET[t]))) for t in SG if t!=s and E(TET[s])&E(TET[t])]
    apex=TET[s][3]; free=[v for v in TET[s] if sum(1 for t in SG if t!=s and v in TET[t])==0]
    print(" ",par(s),s,"辺を分け合う相手",[("ふた" if t in CAP else "帯")+par(t) for t,_ in sh],"誰とも分けない頂点",free,"中心からの距離²",[d2(v,(0,0,0)) for v in free],"C–頂点²",d2(apex,C[i]))
for s in CAP:
    nb=[t for t in ring if E(TET[s])&E(TET[t])]; print("ふた",par(s),"中心",cen[s],"辺を分け合う帯",[par(t) for t in nb],"その炭素の高さ",[sum(cen[t]) for t in nb])
# 鏡映・60 度で全体が自分に重なるか
allT=set(frozenset(T) for T in TET.values())
for n,f in {"鏡映 x↔y":lambda p:(p[1],p[0],p[2]),"60 度（巡回して反転）":lambda p:(-p[2],-p[0],-p[1])}.items():
    print(n,"自分に重なる",set(frozenset(f(v) for v in T) for T in TET.values())==allT,"ふたはふたへ",{frozenset(f(v) for v in TET[c]) for c in CAP}=={frozenset(TET[c]) for c in CAP})
