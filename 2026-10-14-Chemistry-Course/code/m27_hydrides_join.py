# 第2周期の水素化物（枠の B の席に H を入れる）と、枠どうしのつなぎ方
from itertools import product, combinations
BITS=list(product([0,1],repeat=3))
odd=lambda n:sum(n)%2==1
def tet(o,e,want_odd): return [tuple(o[i]+b[i]*e for i in range(3)) for b in BITS if odd(b)==want_odd]
def edges(T): return {frozenset(p) for p in combinations(T,2)}
def frame(off,mirror=False):
    # 枠：大きさ12。A＝番地の和が奇数の小立方体に奇の正四面体。mirror なら偶の小立方体に偶の正四面体
    parts=[];seats=[]
    for n in BITS:
        o=tuple(off[i]+6*n[i] for i in range(3))
        if odd(n)!=mirror: parts.append(tet(o,6,not mirror))
        else: seats.append(tet(o,6,mirror))
    return parts,seats
P,S=frame((0,0,0))
print("== 水素化物：B の席 4 つに H を k 個 ==")
for name,k in [("Ne",0),("HF",1),("H2O",2),("NH3",3),("CH4",4)]:
    sh=[sum(1 for A in P if len(edges(A)&edges(H))==1) for H in S[:k]]
    allsh=set().union(*[edges(A)&edges(H) for H in S[:k] for A in P]) if k else set()
    print(name,"H",k,"空き",4-k,"Hごとの共有辺",sh,"共有辺(重複なし)",len(allsh))
pts=set(v for T in P+S for v in T); print("CH4 の頂点数",len(pts),"辺(重複なし)",len(set().union(*[edges(T) for T in P+S])))
print("O の A の辺",len(set().union(*[edges(T) for T in P])))
print("\n== 枠どうし（2個目は鏡像）のつなぎ方 ==")
for s in [(2,0,0),(1,1,0),(1,1,1),(2,1,1)]:
    off=tuple(6*x for x in s)
    for mir in [False,True]:
        P2,S2=frame(off,mir)
        e1=set().union(*[edges(T) for T in P]); e2=set().union(*[edges(T) for T in P2])
        v1=set(v for T in P for v in T); v2=set(v for T in P2 for v in T)
        same=sum(1 for A in P for B in P2 if set(A)==set(B))
        fill12=sum(1 for B in P2 for T in S if set(B)==set(T)); fill21=sum(1 for A in P for T in S2 if set(A)==set(T))
        print(s,"鏡像" if mir else "同じ向き","共有頂点",len(v1&v2),"共有辺",len(e1&e2),"同じ正四面体",same,"席を埋める 2→1",fill12,"1→2",fill21)
