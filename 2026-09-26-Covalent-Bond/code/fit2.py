import pickle
from fractions import Fraction as Q
exec(open('fit.py').read().split("res=[study")[0])
shapes,res=pickle.load(open('fit.pkl','rb'))
def qmul(x,y): a,b=x;c,d=y; return (a*c+b*d, a*d+b*c+b*d)
def qinv(x):
    a,b=x; n=a*a+a*b-b*b; return (Q(a+b,n),Q(-b,n))
def qadd(x,y): return (x[0]+y[0],x[1]+y[1])
def qsub(x,y): return (x[0]-y[0],x[1]-y[1])
def qdot(u,v):
    s=(Q(0),Q(0))
    for i in range(3): s=qadd(s,qmul(u[i],v[i]))
    return s
def Qv(v): return tuple((Q(c[0]),Q(c[1])) for c in v)
for name,nv,k,info in res:
    if k<3: continue
    n,cc,nn,on,R2=info
    V=[Qv(v) for v in shapes[name]]; n=Qv(n)
    t=qmul((Q(2),Q(0)),qmul((Q(cc[0]),Q(cc[1])),qinv(nn)))
    u=tuple(qmul(t,c) for c in n)
    R2q=(Q(R2[0]),Q(R2[1]))
    # 鏡像 B：v - 2((v·n - cc)/nn) n
    B=[]
    for v in V:
        s=qmul((Q(2),Q(0)),qmul(qsub(qdot(v,n),(Q(cc[0]),Q(cc[1]))),qinv(nn)))
        B.append(tuple(qsub(v[i],qmul(s,n[i])) for i in range(3)))
    T=[tuple(qadd(v[i],u[i]) for i in range(3)) for v in V]          # 平行移動した相手
    Vs=set(V)
    print(f"\n== {name}")
    print("  中心間² = R²:", qdot(u,u)==R2q, "  相手の中心は自分の頂点:", u in Vs)
    print("  鏡像の頂点すべてが相手の中心から R:", all(qdot(tuple(qsub(b[i],u[i]) for i in range(3)),tuple(qsub(b[i],u[i]) for i in range(3)))==R2q for b in B))
    print("  共有する頂点（鏡像）:", len(Vs&set(B)), "  鏡像＝ただの平行移動（ひねり不要）:", set(B)==set(T))
    # 辺の長さと R の比
    E2=min((qdot(tuple(qsub(a[i],b[i]) for i in range(3)),tuple(qsub(a[i],b[i]) for i in range(3))) for a in V for b in V if a!=b),key=lambda x:float(x[0])+float(x[1])*1.618033988749895)
    print("  R² と 辺² :", tuple(map(str,R2q)), tuple(map(str,E2)), "  R＝辺:", R2q==E2)
    shared=[v for v in V if v in set(B)]
    print("  共有する頂点どうしの距離²:", sorted({tuple(map(str,qdot(tuple(qsub(a[i],b[i]) for i in range(3)),tuple(qsub(a[i],b[i]) for i in range(3))))) for a in shared for b in shared if a!=b}))
