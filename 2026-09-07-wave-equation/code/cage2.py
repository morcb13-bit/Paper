import json, numpy as np, sys
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
def krylov_dim(M,v,m):
    L=len(v); rows=[]; piv={}
    x=v.copy()%m; d=0
    for _ in range(L+1):
        r=x.copy()%m
        for p,pr in piv.items():
            if r[p]: r=(r-r[p]*pr)%m
        nz=np.nonzero(r)[0]
        if len(nz)==0: break
        p=int(nz[0]); inv=pow(int(r[p]),m-2,m); r=(r*inv)%m
        piv[p]=r; d+=1
        x=(M@x)%m
    return d
def period(M,v,m,cap):
    init=(v%m).copy(); x=init.copy(); 
    for t in range(1,cap+1):
        x=(M@x)%m
        if np.array_equal(x,init): return t
    return None
C={2:12,3:8,4:6}; FLAT={2:12,3:12,4:12}
print(f"{'環':>4} {'有向辺':>6} {'m':>4} {'クリロフ次元':>12} {'周期':>14}")
for k in (1,2,3):
    for m in (5,7,11,13):
        M,L=mat(k,m,C)
        v=np.zeros(L,dtype=np.int64); v[0]=1
        d=krylov_dim(M,v,m)
        p=period(M,v,m,20000)
        print(f"{k:>4} {L:>6} {m:>4} {d:>12} {str(p) if p else '>2e5':>14}")
    sys.stdout.flush()
