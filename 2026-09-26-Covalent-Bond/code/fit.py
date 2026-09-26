# どんな立体なら「外接球が互いの中心を通る（中心間＝R）」ときに頂点がカチッと揃うか。
# 共有できる頂点は、中心から R/2 の面の上。立体を鏡に写した相手 B は、その面の上の頂点を全部共有する。
# 数は Z[φ]（a+bφ を整数の組）で持ち、判定はすべて整数。
from itertools import product, combinations, permutations
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
def vadd(u,v): return tuple(add(u[i],v[i]) for i in range(3))
def cross(a,b): return (sub(mul(a[1],b[2]),mul(a[2],b[1])),sub(mul(a[2],b[0]),mul(a[0],b[2])),sub(mul(a[0],b[1]),mul(a[1],b[0])))
N=lambda k:(k,0); F=lambda a,b:(a,b)   # a+bφ
Z=N(0); PHI=F(0,1); IPHI=F(-1,1)
def signs(t):
    out=set()
    for s in product([1,-1],repeat=3): out.add(tuple(sc(s[i],t[i]) for i in range(3)))
    return out
def cyc(t): return {t,(t[1],t[2],t[0]),(t[2],t[0],t[1])}
def allperm(t): return set(permutations(t))
def build(gens,perm):
    S=set()
    for g in gens:
        for t in signs(g):
            for c in perm(t): S.add(c)
    return sorted(S)
shapes={
 "正四面体":[t for t in build([(N(1),N(1),N(1))],cyc) if (t[0][0]*t[1][0]*t[2][0])>0],
 "立方体":build([(N(1),N(1),N(1))],cyc),
 "正八面体":build([(N(1),Z,Z)],cyc),
 "正十二面体":build([(N(1),N(1),N(1)),(Z,IPHI,PHI)],cyc),
 "正二十面体":build([(Z,N(1),PHI)],cyc),
 "立方八面体":build([(N(1),N(1),Z)],allperm),
 "切頂八面体":build([(Z,N(1),N(2))],allperm),
 "二十・十二面体":build([(Z,Z,sc(2,PHI)),(N(1),PHI,F(1,1))],cyc),      # 2倍した座標
 "切頂二十面体":build([(Z,N(1),F(0,3)),(N(1),F(2,1),F(0,2)),(PHI,N(2),F(1,2))],cyc),
 "斜方二十・十二面体":build([(N(1),N(1),F(1,2)),(F(1,1),PHI,F(0,2)),(F(2,1),Z,F(1,1))],cyc),
 "切頂十二面体":build([(Z,IPHI,F(2,1)),(IPHI,PHI,F(0,2)),(PHI,N(2),F(1,1))],cyc),
}
phi=(1+5**0.5)/2; val=lambda x:x[0]+x[1]*phi
def study(name,V):
    R2=dot(V[0],V[0]); assert all(dot(v,v)==R2 for v in V), name
    best=None; planes={}
    for a,b,c in combinations(range(len(V)),3):
        n=cross(vsub(V[b],V[a]),vsub(V[c],V[a]))
        if n==(Z,Z,Z): continue
        cc=dot(n,V[a]); nn=dot(n,n)
        if sc(4,mul(cc,cc))!=mul(R2,nn): continue          # 面の距離² = R²/4
        on=frozenset(i for i in range(len(V)) if dot(n,V[i])==cc)
        planes[on]=(n,cc,nn)
    if not planes: return (name,len(V),0,None)
    on,(n,cc,nn)=max(planes.items(),key=lambda t:len(t[0]))
    # 相手の中心 u = 2(cc/nn) n。平行移動だけで重なるか（鏡映＝平行移動）も見る：A+u の頂点が鏡像の頂点と一致するか
    # 整数のまま：座標を nn 倍して比べる
    k=len(on)
    return (name,len(V),k,(n,cc,nn,on,R2))
res=[study(n,V) for n,V in shapes.items()]
for name,nv,k,info in sorted(res,key=lambda r:-r[2]):
    print(f"{name:10s} 頂点{nv:3d}  R/2 の面に乗る頂点（＝共有できる頂点）の最大 {k}")
import pickle; pickle.dump((shapes,res),open('fit.pkl','wb'))
