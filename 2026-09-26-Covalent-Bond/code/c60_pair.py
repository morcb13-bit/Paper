# 切頂二十面体（C60）を六角形の中心について反対へ写して組む。判定は Z[φ] の整数のみ。
from itertools import product
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def sub(x,y): return (x[0]-y[0],x[1]-y[1])
def mul(x,y): a,b=x;c,d=y; return (a*c+b*d, a*d+b*c+b*d)
def sc(k,x): return (k*x[0],k*x[1])
def dot(u,v):
    s=(0,0)
    for i in range(3): s=add(s,mul(u[i],v[i]))
    return s
def sign(x):   # a+bφ の符号を整数だけで：2(a+bφ) = (2a+b) + b√5
    u=2*x[0]+x[1]; v=x[1]
    if u>=0 and v>=0: return 0 if (u==0 and v==0) else 1
    if u<=0 and v<=0: return -1
    return 1 if (u*u>5*v*v)==(u>0) else -1
Z=(0,0);O=(1,0);P=(0,1)
def cyc(t): return [t,(t[1],t[2],t[0]),(t[2],t[0],t[1])]
def signs(t):
    out=set()
    for s in product([1,-1],repeat=3):
        out.add(tuple(sc(s[i],t[i]) for i in range(3)))
    return out
def c60():
    V=set()
    for base in [(Z,O,sc(3,P)),(O,(2,1),sc(2,P)),(P,(2,0),(1,2))]:   # 2+φ, φ³=1+2φ
        for t in signs(base):
            for c in cyc(t): V.add(c)
    return sorted(V), (4,0)            # 辺の長さ2 → 二乗4
def dodeca():
    V=set()
    for s in product([O,sc(-1,O)],repeat=3): V.add(tuple(s))
    IP=(-1,1)
    for t in signs((Z,IP,P)):
        for c in cyc(t): V.add(c)
    return sorted(V), (8,-4)
def graph(V,E2):
    E=set()
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            w=tuple(sub(V[i][k],V[j][k]) for k in range(3))
            if dot(w,w)==E2: E.add((i,j))
    adj={i:set() for i in range(len(V))}
    for i,j in E: adj[i].add(j); adj[j].add(i)
    return E,adj
def cycles(adj,L):
    res=set()
    def dfs(s,cur,path):
        for w in adj[cur]:
            if len(path)==L:
                if w==s: res.add(frozenset(path))
            elif w>s and w not in path: dfs(s,w,path+[w])
    for s in adj: dfs(s,s,[s])
    return list(res)

def test(name,V,E2,L):
    print(f"\n=== {name}")
    E,adj=graph(V,E2)
    faces=cycles(adj,L)
    print(f"照合 頂点{len(V)} 辺{len(E)} 長さ{L}の輪{len(faces)}")
    F=sorted(faces,key=lambda f:sorted(f))[0]
    S=sorted(F)
    # 6倍の座標で中心の2倍 = 2·6·(Σx/L) を整数に保つ：L=6 なら 2Σx、L=5 なら 12Σx/5（割り切れを検査）
    s=(Z,Z,Z)
    for i in S: s=tuple(add(s[k],V[i][k]) for k in range(3))
    num=tuple(sc(60,c) for c in s)
    if any(c[0]%L or c[1]%L for c in num):
        print("中心が6倍座標の格子に乗らない（写しても頂点は格子外）→ T1 NG"); return
    twoP=tuple((c[0]//L,c[1]//L) for c in num)
    X=[tuple(sc(30,c) for c in v) for v in V]
    Y=[tuple(sub(twoP[k],x[k]) for k in range(3)) for x in X]
    idx={x:i for i,x in enumerate(X)}
    Fset=set(S)
    img=[idx.get(y) for y in Y]
    T1=all(img[i] is not None and img[i] in Fset for i in S)
    common=[i for i,y in enumerate(Y) if y in idx]
    T2=set(idx[Y[i]] for i in common)==Fset and len(common)==len(S)
    print("T1 面の頂点が面の頂点へ写る:", "OK" if T1 else "NG")
    print("T2 A と B で重なる頂点:", len(common), "OK" if T2 else "NG")
    if not (T1 and T2): return
    # 合成した頂点集合：A の 60 個＋B の新しい 54 個
    allv=list(X); bidx={}
    for i,y in enumerate(Y):
        if y in idx: bidx[i]=idx[y]
        else: bidx[i]=len(allv); allv.append(y)
    EA={(i,j) for i,j in E}
    EAll=set(EA)|{tuple(sorted((bidx[i],bidx[j]))) for i,j in E}
    N=len(allv)
    # 電子の写像 σ：A の頂点 i → B の bidx[i]、B の頂点 bidx[i] → i
    sig=[None]*N
    for i in range(len(V)): sig[i]=bidx[i]; sig[bidx[i]]=i
    T3=all(sig[sig[v]]==v for v in range(N)) and all(tuple(sorted((sig[a],sig[b]))) in EAll for a,b in EAll)
    fixed=[v for v in range(N) if sig[v]==v]
    print(f"合成 頂点{N} 辺{len(EAll)}")
    print("T3 辺を辺へ・二回で元へ:", "OK" if T3 else "NG")
    print("T4 動かない頂点:", len(fixed), "OK" if not fixed else "NG")
    # T5 法線（中心→面の中心 = s 方向）で高さを比べる。面の高さ c = s·x（面の頂点で共通）
    c=dot(s,X[S[0]])
    assert all(dot(s,X[i])==c for i in S)
    a_side=[sign(sub(dot(s,x),c)) for x in X]
    b_side=[sign(sub(dot(s,y),c)) for y in Y]
    T5=all(v<=0 for v in a_side) and all(v>=0 for v in b_side) \
       and sum(v==0 for v in a_side)==len(S) and sum(v==0 for v in b_side)==len(S)
    print("T5 共有面以外で重ならない:", "OK" if T5 else "NG",
          f"(A 面上{sum(v==0 for v in a_side)} 片側{sum(v<0 for v in a_side)} / B 面上{sum(v==0 for v in b_side)} 反対側{sum(v>0 for v in b_side)})")
    # 共有面の上での電子の相手
    print("共有六角形の上の対:", sorted({tuple(sorted((i,sig[i]))) for i in S}))

V,E2=c60()
test("C60 × 六角形の中心",V,E2,6)
test("C60 × 五角形の中心（対照）",V,E2,5)
D,E2d=dodeca()
test("正12面体 × 五角形の中心（対照）",D,E2d,5)
