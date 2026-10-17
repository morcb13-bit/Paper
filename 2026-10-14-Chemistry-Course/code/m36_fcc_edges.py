# 面心立方：小立方体すべてに正四面体（奇は A 向き、偶は B 向き）。内側の正四面体が分け合う辺・点を数える
from itertools import product, combinations
BITS=list(product([0,1],repeat=3)); odd=lambda n:sum(n)%2==1
def tet(n): return [tuple(6*n[i]+6*b[i] for i in range(3)) for b in BITS if odd(b)==odd(n)]
E=lambda T:{frozenset(p) for p in combinations(T,2)}
N=5; T={n:tet(n) for n in product(range(N),repeat=3)}
c=(2,2,2)
sh=[(m,len(E(T[c])&E(T[m])),len(set(T[c])&set(T[m]))) for m in T if m!=c]
print("内側の正四面体：辺を分け合う相手",sum(1 for m,e,v in sh if e),"分け合う辺",sum(e for m,e,v in sh),"点だけで接する相手",sum(1 for m,e,v in sh if v and not e))
from collections import Counter
cnt=Counter(e for t in T.values() for e in E(t))
inner=[e for e in E(T[c])]; print("その 6 辺それぞれを持つ正四面体の数",[cnt[e] for e in inner])
