# 転がる電子（正四面体・一辺 6）の向きは、入っている立方体の番地で決まるか ── 整数だけ
# 6 辺すべてで倒れる手を 6 歩先まで全部たどり、(立方体, 向き) の組を集める。
from itertools import product, combinations
BITS=list(product([0,1],repeat=3)); odd=lambda n:sum(n)%2==1
def tet(o,e,w): return [tuple(o[i]+b[i]*e for i in range(3)) for b in BITS if odd(b)==w]
def cube(T): return tuple(min(p[i] for p in T)//6 for i in range(3))
def orient(T):  # 立方体の角 (0,0,0) を含むなら B（偶の角）、含まなければ A
    c=tuple(6*x for x in cube(T)); return 'B' if c in T else 'A'
start=tet((0,0,0),6,False)
seen={}; front=[start]; moves=set()
for step in range(6):
    nxt=[]
    for T in front:
        c0=cube(T)
        for u,v in combinations(T,2):
            S=[tuple(u[i]+v[i]-p[i] for i in range(3)) if p not in (u,v) else p for p in T]
            c=cube(S); d=tuple(c[i]-c0[i] for i in range(3)); moves.add(d)
            key=c
            seen.setdefault(key,set()).add(orient(S)); nxt.append(S)
    front={tuple(sorted(S)):S for S in nxt}.values()
bad=[c for c,o in seen.items() if len(o)>1 or o!={'B' if sum(c)%2==0 else 'A'}]
print('一歩で動く向き', sorted(moves))
print('たどった立方体', len(seen), ' 向きが番地の偶奇で決まらない立方体', len(bad))
