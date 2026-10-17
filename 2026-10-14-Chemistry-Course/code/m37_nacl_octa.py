# 塩化ナトリウム：Cl を中心とする正八面体（頂点は Na 6 個）が隣と分け合う辺
from itertools import product, combinations
def octa(c): return [tuple(c[i]+d*(i==k)*6 for i in range(3)) for k in range(3) for d in (1,-1)]
def E(O):
    return {frozenset((p,q)) for p,q in combinations(O,2) if sum((p[i]-q[i])**2 for i in range(3))==72}
cl=[c for c in product(range(-12,31,6),repeat=3) if (sum(c)//6)%2==1]
c0=(6,6,6); O0=E(octa(c0))
nb=[(c,len(O0&E(octa(c)))) for c in cl if c!=c0]
print("正八面体の辺",len(O0),"辺を分け合う隣",sum(1 for c,k in nb if k),"相手ごとの本数",sorted({k for c,k in nb if k}))
print("Na-Cl（立方体の辺）の長さ²",36,"正八面体の辺（Na-Na）の長さ²",72)
