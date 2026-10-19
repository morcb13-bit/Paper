# 速い道具：辺の長さ² 288（2 倍の単位）の向きの表から正四面体を組む
from itertools import product, combinations
dot=lambda a,b:sum(x*y for x,y in zip(a,b)); add=lambda a,b:tuple(x+y for x,y in zip(a,b)); sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
neg=lambda a:tuple(-x for x in a)
D=[d for d in product(range(-16,17),repeat=4) if dot(d,d)==288]
def tets_with_edge(p,q):
    """辺 p–q を持つ正四面体：(中心, 頂点の frozenset)"""
    R=[add(p,d) for d in D if dot(sub(add(p,d),q),sub(add(p,d),q))==288]
    out=set()
    for r,t in combinations(R,2):
        if dot(sub(r,t),sub(r,t))!=288: continue
        s=add(add(p,q),add(r,t))
        if all(x%4==0 for x in s): out.add((tuple(x//4 for x in s),frozenset((p,q,r,t))))
    return list(out)
