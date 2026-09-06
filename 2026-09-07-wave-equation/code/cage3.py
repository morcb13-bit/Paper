import json, numpy as np, sys, math
def mat(k,m,coef):
    A=[list(a) for a in json.load(open(f'cage_{k}.json'))['adj']]
    N=len(A); DEG=[len(a) for a in A]
    E={}; el=[]
    for u in range(N):
        for v in A[u]: E[(u,v)]=len(el); el.append((u,v))
    IN=[[] for _ in range(N)]
    for i,(u,v) in enumerate(el): IN[v].append(i)
    REV=[E[(v,u)] for (u,v) in el]
    L=len(el); M=np.zeros((L,L),dtype=np.int64)
    for v in range(N):
        c=coef[DEG[v]]
        for i in IN[v]:
            r=REV[i]
            for j in IN[v]: M[r,j]=(M[r,j]+c)%m
            M[r,i]=(M[r,i]-12)%m
    return M%m, L
def inv_mod(M,m):
    L=len(M); A=np.concatenate([M%m, np.eye(L,dtype=np.int64)],axis=1)%m
    for c in range(L):
        p=None
        for r in range(c,L):
            if A[r,c]%m: p=r; break
        assert p is not None,"特異"
        if p!=c: A[[c,p]]=A[[p,c]]
        A[c]=(A[c]*pow(int(A[c,c]),m-2,m))%m
        col=A[:,c].copy(); col[c]=0
        A=(A-np.outer(col,A[c]))%m
    return A[:,L:]%m
def mpow(M,e,m):
    R=np.eye(len(M),dtype=np.int64); B=M%m
    while e:
        if e&1: R=(R@B)%m
        B=(B@B)%m; e>>=1
    return R
def bsgs(M,v,m,B=40000):
    Mi=inv_mod(M,m)
    tab={}; x=(v%m).copy()
    for j in range(B):
        tab.setdefault(x.tobytes(), j)
        x=(Mi@x)%m
    G=mpow(M,B,m); y=(v%m).copy()
    for i in range(0,B+1):
        j=tab.get(y.tobytes())
        if j is not None:
            t=i*B+j
            if t>0: return t
        y=(G@y)%m
    return None
C={2:12,3:8,4:6}; FLAT={2:12,3:12,4:12}
def fac(n):
    f={};d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return " · ".join(f"{p}^{e}" if e>1 else str(p) for p,e in f.items())
which=sys.argv[1]; ks=[int(x) for x in sys.argv[2].split(',')]
coef = C if which=='C' else FLAT
print(f"係数 {'24/d' if which=='C' else '全次数12（対照）'}")
print(f"{'環':>4} {'有向辺':>6} {'m':>4} {'周期':>14}  素因数分解")
for k in ks:
    for m in (5,7,11,13):
        M,L=mat(k,m,coef)
        v=np.zeros(L,dtype=np.int64); v[0]=1
        try:
            p=bsgs(M,v,m)
        except AssertionError:
            print(f"{k:>4} {L:>6} {m:>4} {'特異（逆が無い）':>14}"); continue
        print(f"{k:>4} {L:>6} {m:>4} {str(p):>14}  {fac(p) if p else ''}")
        sys.stdout.flush()
