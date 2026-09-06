"""小さい担体を作る。rings_integer.json の先頭 k 環。座標は素性の確認にしか使わない。"""
import json, collections
from b13_two_tilings import zadd, zt, to_xy
from b13_layers import LayeredFigure

R=[tuple(int(x) for x in c) for c in json.load(open('rings_integer.json'))['rings']]
def build(k):
    F=LayeredFigure()
    for c in R[:k]: F.add_ring(c)
    pent={q:[zadd(q,zt(2*j+(a%2))) for j in range(5)] for q,a in F.cells.items()}
    V={}; adj={}
    def vid(v):
        if v not in V: V[v]=len(V); adj[V[v]]=set()
        return V[v]
    for q,cor in pent.items():
        ids=[vid(c) for c in cor]
        for i in range(5):
            a,b=ids[i],ids[(i+1)%5]; adj[a].add(b); adj[b].add(a)
    A=[sorted(adj[i]) for i in range(len(V))]
    deg=collections.Counter(len(a) for a in A)
    return A, len(F.cells), deg
for k in (1,2,3,5,10,21,60):
    A,nc,deg=build(k)
    E=sum(len(a) for a in A)
    print(f"環{k:>3}: 五角形{nc:>4} 頂点{len(A):>5} 有向辺{E:>6} 次数{dict(sorted(deg.items()))} "
          f"24を割らない次数: {[d for d in deg if 24%d]}")
    json.dump({'adj':A}, open(f'cage_{k}.json','w'))
