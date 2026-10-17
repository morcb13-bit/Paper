# 水 2 個：水1の H と 水2 の空いた席（孤立電子対）がどう接するか
from itertools import product, combinations
B=list(product([0,1],repeat=3)); odd=lambda n:sum(n)%2==1
def tet(o,e,w): return [tuple(o[i]+b[i]*e for i in range(3)) for b in B if odd(b)==w]
E=lambda T:{frozenset(p) for p in combinations(T,2)}
def water(off,mirror,hseats):
    cubes={}
    parts=[];seats={}
    for n in B:
        o=tuple(off[i]+6*n[i] for i in range(3))
        if odd(n)!=mirror: parts.append(tet(o,6,not mirror))
        else: seats[n]=tet(o,6,mirror)
    H=[seats[n] for n in hseats]; LP=[seats[n] for n in seats if n not in hseats]
    cubes={tuple(off[i]//6+n[i] for i in range(3)) for n in B}
    return parts,H,LP,cubes
def contact(X,Y):
    v=len(set(a for T in X for a in T)&set(a for T in Y for a in T)); e=len(set().union(*map(E,X))&set().union(*map(E,Y))) if X and Y else 0
    return (v,e)
P1,H1,L1,C1=water((0,0,0),False,[(0,0,0),(1,1,0)])
print("水1 内部：O-H の共有（頂点,辺）",[contact(P1,[h]) for h in H1])
# 水2：鏡像、(-12,-12,-12) に置き、席 (1,1,1) を空ける
P2,H2,L2,C2=water((-12,-12,-12),True,[(1,0,0),(0,1,0)])
print("水2 が鏡像のとき  H1-LP2",contact(H1,L2),"H1-O2",contact(H1,P2),"O1-O2",contact(P1,P2),"O1-LP2",contact(P1,L2),"H1-H2",contact(H1,H2),"LP1-?",contact(L1,P2+H2+L2))
print("共有点",set(a for T in H1 for a in T)&set(a for T in L2 for a in T))
# 一覧：重ならない隣の置き方すべて（ずらし -2..2、鏡像あり/なし、H の席は全 6 通り）
from collections import Counter
cnt=Counter()
for s in product(range(-2,3),repeat=3):
    if all(abs(x)<2 for x in s): continue
    if max(abs(x) for x in s)>2: continue
    for mir in (False,True):
        seats=[n for n in B if odd(n)==mir]
        for hs in combinations(seats,2):
            P2,H2,L2,C2=water(tuple(6*x for x in s),mir,list(hs))
            if C1&C2: continue
            r=(contact(H1,L2),contact(H1,P2),contact(P1,P2),contact(H1,H2),contact(P1,L2))
            if r[0]!=(0,0) and r[1:]==((0,0),)*4: cnt[("鏡像" if mir else "同じ向き",r[0])]+=1
print("H1 と LP2 だけが接する置き方",dict(cnt))
