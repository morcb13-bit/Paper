import json, numpy as np, sys
exec(open('cage3.py').read().split("C={2:12")[0])
def reduce_min(M,v,m,t):
    f={};n=t;d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    for p in list(f):
        while t%p==0 and np.array_equal((mpow(M,t//p,m)@(v%m))%m, v%m):
            t//=p
    return t
def facs(n):
    f={};d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return " · ".join(f"{p}^{e}" if e>1 else str(p) for p,e in f.items())
def kdim(M,v,m):
    L=len(v); piv={}; x=(v%m).copy(); d=0
    for _ in range(L+1):
        r=x.copy()%m
        for p,pr in piv.items():
            if r[p]: r=(r-int(r[p])*pr)%m
        nz=np.nonzero(r)[0]
        if len(nz)==0: break
        p=int(nz[0]); r=(r*pow(int(r[p]),m-2,m))%m; piv[p]=r; d+=1
        x=(M@x)%m
    return d
C={2:12,3:8,4:6}; FLAT={2:12,3:12,4:12}
print("係数 24/d ── 1環（有向辺80、頂点30、次数 2:10 / 3:20）")
print(f"{'m':>4} {'クリロフ次元':>10} {'最小周期':>12}  素因数分解")
for m in (5,7,11,13,17,19,23):
    M,L=mat(1,m,C); v=np.zeros(L,dtype=np.int64); v[0]=1
    d=kdim(M,v,m)
    t=bsgs(M,v,m,B=40000)
    p=reduce_min(M,v,m,t) if t else None
    print(f"{m:>4} {d:>10} {str(p):>12}  {facs(p) if p else '>1.6e9'}")
    sys.stdout.flush()
print("\n対照 全次数 c=12 ── 1環")
for m in (5,7,11,13):
    M,L=mat(1,m,FLAT); v=np.zeros(L,dtype=np.int64); v[0]=1
    d=kdim(M,v,m)
    try: t=bsgs(M,v,m,B=40000)
    except AssertionError: print(f"{m:>4} {d:>10} 逆写像なし"); continue
    p=reduce_min(M,v,m,t) if t else None
    print(f"{m:>4} {d:>10} {str(p):>12}  {facs(p) if p else '>1.6e9'}")
    sys.stdout.flush()
print("\n初期の辺を変える（24/d, 1環）")
for m in (5,7):
    M,L=mat(1,m,C)
    for e in (0,1,7,23,50):
        v=np.zeros(L,dtype=np.int64); v[e]=1
        t=bsgs(M,v,m,B=40000); p=reduce_min(M,v,m,t) if t else None
        print(f"  m={m} 辺{e:>3}: クリロフ次元 {kdim(M,v,m):>3}  最小周期 {p}")
    sys.stdout.flush()
