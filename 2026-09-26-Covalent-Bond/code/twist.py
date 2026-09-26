# 外接球が互いの中心を通る二つの正十二面体を、頂点で接するようにひねって組む（Z[φ] の整数のみ）
from itertools import product, combinations
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def sub(x,y): return (x[0]-y[0],x[1]-y[1])
def mul(x,y): a,b=x;c,d=y; return (a*c+b*d, a*d+b*c+b*d)
def neg(x): return (-x[0],-x[1])
def sc(k,x): return (k*x[0],k*x[1])
def dot(u,v):
    s=(0,0)
    for i in range(3): s=add(s,mul(u[i],v[i]))
    return s
def vsub(u,v): return tuple(sub(u[i],v[i]) for i in range(3))
def cross(a,b): return (sub(mul(a[1],b[2]),mul(a[2],b[1])),sub(mul(a[2],b[0]),mul(a[0],b[2])),sub(mul(a[0],b[1]),mul(a[1],b[0])))
Z=(0,0);O1=(1,0);P=(0,1);IP=(-1,1)
V=set()
for s in product([O1,neg(O1)],repeat=3): V.add(s)
for a in [IP,neg(IP)]:
    for b in [P,neg(P)]:
        base=[Z,a,b]
        for k in range(3): V.add(tuple(base[(i-k)%3] for i in range(3)))
V=sorted(V); R2=(3,0)
# 共有できる頂点は二つの外接球の両方に乗る＝中心を結ぶ線の垂直二等分面の上。
# その面は A の中心から R/2 の距離：(n·v)² ×4 ＝ 3|n|²
best={}
for a,b,c in combinations(range(20),3):
    n=cross(vsub(V[b],V[a]),vsub(V[c],V[a]))
    if n==(Z,Z,Z): continue
    cc=dot(n,V[a]); nn=dot(n,n)
    if sc(4,mul(cc,cc))!=sc(3,nn): continue
    on=frozenset(i for i in range(20) if dot(n,V[i])==cc)
    best[on]=(n,cc)
sizes=sorted({len(s) for s in best})
print("R/2 の面に乗る頂点の数の種類:",sizes, " 面の数:",len(best))
for k in sizes: print(f"  {k} 頂点の面: {sum(len(s)==k for s in best)} 枚")

# 2点で接する場合：g = v1·v2 ごとに、陽子Bの中心 u が存在するか（|u|²=3, u·v1=u·v2=3/2）。
# 存在 ⇔ 9+6g > 0（Z[φ] の符号判定）。作り方は、中心を結ぶ線の垂直二等分面での鏡映（正十二面体は鏡映で自分と同じ形）。
def sgn(x):
    u=2*x[0]+x[1]; v=x[1]
    if u>=0 and v>=0: return 0 if u==0 and v==0 else 1
    if u<=0 and v<=0: return -1
    return 1 if (u*u>5*v*v)==(u>0) else -1
from collections import Counter
cnt=Counter()
for i,j in combinations(range(20),2):
    g=dot(V[i],V[j]); d2=sub((6,0),sc(2,g)); cnt[(d2,g)]+=1
print("\n頂点の組（距離の二乗 d², 内積 g）ごと：")
for (d2,g),n in sorted(cnt.items(),key=lambda t:t[0][0][0]+t[0][0][1]*1.618):
    ok=sgn(add((9,0),sc(6,g)))>0
    print(f"  d²={d2}  g={g}  組の数 {n:2d}  2点で接する配置: {'あり' if ok else 'なし'}")

# 表示用の確認（浮動小数）：各種類から一組取って作り、重なる頂点がちょうど2個か
import numpy as np
phi=(1+5**0.5)/2; F=lambda x:x[0]+x[1]*phi
A=np.array([[F(c) for c in v] for v in V])
done=set()
for i,j in combinations(range(20),2):
    g=dot(V[i],V[j])
    if g in done or sgn(add((9,0),sc(6,g)))<=0: continue
    done.add(g)
    v1,v2=A[i],A[j]; gg=v1@v2
    al=3/(2*(3+gg)); c=np.cross(v1,v2); be=np.sqrt((3-al*al*(6+2*gg))/(c@c))
    u=al*(v1+v2)+be*c
    B=np.array([x-2*((x@u-1.5)/3)*u for x in A])
    sh=sum(np.min(np.linalg.norm(B-a,axis=1))<1e-9 for a in A)
    print(f"  g={F(g):+.3f}: |u|²={u@u:.12f}  B の頂点と中心の距離²={np.unique(np.round(np.sum((B-u)**2,1),9))}  重なる頂点 {sh}")
np.save('twist_u.npy',u)
