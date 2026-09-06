"""檻の大きさ：剰余 m での一歩の周期。
   gcd(m,12)=1 で一歩は可逆 → 全状態が巡回。測るのは周期（整数）。
   必ず落ちる設定：c を全次数 12 にすると周期が変わること。"""
import json, collections, math, sys

def load(k):
    A=[list(a) for a in json.load(open(f'cage_{k}.json'))['adj']]
    N=len(A); DEG=[len(a) for a in A]
    E={}; el=[]
    for u in range(N):
        for v in A[u]: E[(u,v)]=len(el); el.append((u,v))
    IN=[[] for _ in range(N)]
    for i,(u,v) in enumerate(el): IN[v].append(i)
    REV=[E[(v,u)] for (u,v) in el]
    return A,N,DEG,el,IN,REV

def period(k, m, e0, coef, cap=10**7):
    A,N,DEG,el,IN,REV=load(k)
    L=len(el)
    psi=[0]*L; psi[e0]=1%m
    init=tuple(psi)
    t=0
    while t<cap:
        new=[0]*L
        for v in range(N):
            S=0
            for i in IN[v]: S+=psi[i]
            S%=m
            c=coef[DEG[v]]
            for i in IN[v]: new[REV[i]]=(c*S-12*psi[i])%m
        psi=new; t+=1
        if tuple(psi)==init: return t
    return None

C={2:12,3:8,4:6}
FLAT={2:12,3:12,4:12}
print("周期（初期＝辺0に1、他0）")
print(f"{'環':>4} {'有向辺':>6} " + " ".join(f"m={m:<10}" for m in (5,7,11,13)))
for k in (1,2,3,5,10):
    A,N,DEG,el,IN,REV=load(k)
    row=[]
    for m in (5,7,11,13):
        p=period(k,m,0,C,cap=2_000_000)
        row.append(str(p) if p else ">2e6")
    print(f"{k:>4} {len(el):>6} " + " ".join(f"{x:<12}" for x in row))
