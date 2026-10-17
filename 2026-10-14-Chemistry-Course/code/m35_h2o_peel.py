# H2O の H の席：O の正四面体と分け合う辺・点、H どうしの共有、残る辺の候補
from itertools import product, combinations
BITS=list(product([0,1],repeat=3)); odd=lambda n:sum(n)%2==1
def tet(o,e,want_odd): return [tuple(o[i]+b[i]*e for i in range(3)) for b in BITS if odd(b)==want_odd]
E=lambda T:{frozenset(p) for p in combinations(T,2)}
A={n:tet(tuple(6*x for x in n),6,True) for n in BITS if odd(n)}
SEAT={n:tet(tuple(6*x for x in n),6,False) for n in BITS if not odd(n)}
H=[(0,0,0),(1,1,0)]
for h in H:
    T=SEAT[h]; print("H",h,"頂点",T)
    for a,At in A.items():
        s=E(T)&E(At)
        if s: print("  A",a,"共有辺",[sorted(x) for x in s])
hv=set(SEAT[H[0]])&set(SEAT[H[1]]); print("H どうしの共有点",hv)
for h in H:
    for a,At in A.items():
        for e in E(SEAT[h])&E(At):
            both=sum(1 for g in H if E(SEAT[g])&E(At))
            print("H",h,"辺",sorted(e),"A",a,"この A と辺を分ける H の数",both,"H-H 共有点に触れる",bool(e&hv))
