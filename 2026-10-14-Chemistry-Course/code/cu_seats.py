# シェルピンスキーの銅原子を導線に並べたとき、1 原子に新しく座れる点の数（整数だけ）
from itertools import product
mid=lambda a,b:tuple((a[i]+b[i])//2 for i in range(3))
split=lambda t:[[t[i]]+[mid(t[i],t[j]) for j in range(4) if j!=i] for i in range(4)]
def even(o,s): return [tuple(o[i]+s*c[i] for i in range(3)) for c in product((0,1),repeat=3) if (sum(o)//s+sum(c))%2==0]
def points(o,s,m):
    T=[even(o,s)]
    for _ in range(m-1): T=[u for v in T for u in split(v)]
    return {p for t in T for p in t}
S=96
for m in (1,2,3,4):
    A=points((0,0,0),S,m); B=points((S,0,0),S,m)
    shared=len(A&B)
    new=len(B-A)                       # 2 個目の原子が新しく持ち込む点
    print(f"{m} 段：点 {len(A):3d}　隣と分け合う点 {shared}　並べたとき 1 原子に新しく座れる点 {new}"+("　← 銅の電子 29" if new==29 else ""))
