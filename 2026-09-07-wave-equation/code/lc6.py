"""検定LC6 τ層上の整数波伝播。座標は使わない。
   辺は頭（行き先）の τ で層に割り振る。"""
import json, collections
from fractions import Fraction as Fr
G=json.load(open('carrier_1245_graph.json'))
ADJ=[list(a) for a in G['adj']]; N=len(ADJ)
DEG=[len(a) for a in ADJ]; C={2:12,3:8,4:6}
src=int(open('center.txt').read())
d=[-1]*N; d[src]=0; q=collections.deque([src])
while q:
    u=q.popleft()
    for v in ADJ[u]:
        if d[v]<0: d[v]=d[u]+1; q.append(v)
E={}; el=[]
for u in range(N):
    for v in ADJ[u]: E[(u,v)]=len(el); el.append((u,v))
IN=[[] for _ in range(N)]; OUT=[[] for _ in range(N)]
for k,(u,v) in enumerate(el): IN[v].append(k); OUT[u].append(k)
REV=[E[(v,u)] for (u,v) in el]
HEAD=[d[v] for (u,v) in el]

def fwd(psi):
    new=[0]*len(el)
    for v in range(N):
        S=0; nz=False
        for k in IN[v]:
            if psi[k]: nz=True
            S+=psi[k]
        if not nz: continue
        c=C[DEG[v]]
        for k in IN[v]: new[REV[k]]=c*S-12*psi[k]
    return new

def bwd(psi):
    old=[0]*len(el)
    for v in range(N):
        S12=sum(psi[k] for k in OUT[v])       # Σout = 12·S
        if S12==0 and all(psi[k]==0 for k in OUT[v]): continue
        assert S12 % 12 == 0
        S=S12//12; c=C[DEG[v]]
        for k in OUT[v]:
            num=c*S-psi[k]
            assert num % 12 == 0
            old[REV[k]]=num//12
    return old

e0=E[(src,ADJ[src][0])]
psi=[0]*len(el); psi[e0]=1
init=psi[:]
T=24
print(f"出発 辺 {el[e0]}  τ(頭)={HEAD[e0]}")
print(f"{'t':>3} {'非零辺':>7} {'全体の二乗和':>16} {'/144^t':>8} {'全偶数':>6}  層ごとの (τ:本数)")
for t in range(0,T+1):
    if t: psi=fwd(psi)
    nz=sum(1 for x in psi if x)
    ss=sum(x*x for x in psi)
    lay=collections.Counter(HEAD[k] for k in range(len(el)) if psi[k])
    ok = (ss == 144**t)
    ev = all(x%2==0 for x in psi) if t else "-"
    s=" ".join(f"{k}:{lay[k]}" for k in sorted(lay))
    print(f"{t:>3} {nz:>7} {ss:>16} {str(ok):>8} {str(ev):>6}  {s}")

# 層ごとの二乗和（最終歩）
lay2=collections.defaultdict(int)
for k in range(len(el)):
    if psi[k]: lay2[HEAD[k]]+=psi[k]*psi[k]
tot=sum(lay2.values())
print(f"\nt={T} 層ごとの二乗和（合計 {tot}、144^{T} と一致 {tot==144**T}）")
for k in sorted(lay2):
    print(f"  τ={k:>3}  本数 {sum(1 for j in range(len(el)) if psi[j] and HEAD[j]==k):>5}  二乗和 {lay2[k]}")

# 逆向き再構成
back=psi[:]
for _ in range(T): back=bwd(back)
print(f"\n逆向きに {T} 歩戻して初期状態と厳密に一致: {back==init}")
