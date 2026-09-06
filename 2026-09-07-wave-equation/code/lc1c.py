"""LC1c：前線のさざなみの絶対幅が半径とともに増えるかを見る。
   OK なら：幅が一定 → 前線は円（ずれは有限）。
   NG なら：幅が半径に比例 → 前線は多角形。"""
import json, math, collections

G = json.load(open('carrier_1245_graph.json'))
XY = [complex(x, y) for x, y in G['xy']]
ADJ = [list(a) for a in G['adj']]
N = len(XY)
center = min(range(N), key=lambda i: abs(XY[i]))

def bfs(src, adj):
    d = [-1]*len(adj); d[src]=0
    q = collections.deque([src])
    while q:
        u=q.popleft()
        for v in adj[u]:
            if d[v]<0: d[v]=d[u]+1; q.append(v)
    return d

def profile(tau, rad, xy, nb, tlist):
    """歩数 t の前線半径 R(θ,t) を方位区間ごとに取り、絶対幅と平均を返す"""
    out=[]
    for t in tlist:
        R=[0.0]*nb
        for i in range(len(tau)):
            if 0<=tau[i]<=t and rad[i]>R[i2b[id(xy)][i]] if False else False: pass
        out.append(None)
    return out

def binidx(xy, nb):
    return [int((math.atan2(p.imag,p.real)%(2*math.pi))/(2*math.pi)*nb)%nb for p in xy]

def ripple(tau, rad, bins, nb, tlist):
    rows=[]
    for t in tlist:
        R=[0.0]*nb
        for i in range(len(tau)):
            if 0<=tau[i]<=t and rad[i]>R[bins[i]]: R[bins[i]]=rad[i]
        me=sum(R)/nb
        rows.append((t, me, max(R)-min(R), (max(R)-min(R))/me))
    return rows

NB=20
TAU=bfs(center,ADJ); RAD=[abs(p) for p in XY]; B=binidx(XY,NB)
print("ペンローズ担体（区間20）")
print(f"{'t':>4} {'平均R':>8} {'絶対幅':>8} {'相対幅':>8}")
for t,me,ab,rl in ripple(TAU,RAD,B,NB,[20,30,40,50,60,70,80,90,100]):
    print(f"{t:>4} {me:8.2f} {ab:8.3f} {rl:8.4f}")

def lattice_square(S):
    def sid(x,y): return (x+S)*(2*S+1)+(y+S)
    adj=[[] for _ in range((2*S+1)**2)]; xy=[0j]*len(adj)
    for x in range(-S,S+1):
        for y in range(-S,S+1):
            i=sid(x,y); xy[i]=complex(x,y)
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                if -S<=x+dx<=S and -S<=y+dy<=S: adj[i].append(sid(x+dx,y+dy))
    return adj, xy, sid(0,0)

sadj,sxy,sc = lattice_square(110)
stau=bfs(sc,sadj); srad=[abs(p) for p in sxy]; sb=binidx(sxy,NB)
print("\n負の対照 正方格子（区間20）")
print(f"{'t':>4} {'平均R':>8} {'絶対幅':>8} {'相対幅':>8}")
for t,me,ab,rl in ripple(stau,srad,sb,NB,[20,30,40,50,60,70,80,90,100]):
    print(f"{t:>4} {me:8.2f} {ab:8.3f} {rl:8.4f}")

# c(θ) の形：10回対称か
def cprof(tau,rad,xy,nb,rmin,rmax):
    c=[0.0]*nb
    bb=binidx(xy,nb)
    for i in range(len(tau)):
        if tau[i]>0 and rmin<=rad[i]<=rmax: c[bb[i]]=max(c[bb[i]],rad[i]/tau[i])
    return c
c=cprof(TAU,RAD,XY,40,60,107)
print("\nc(θ) の形（40区間・9度刻み・r=60〜107）")
for k in range(0,40,4):
    print("  " + " ".join(f"{c[j]:.4f}" for j in range(k,k+4)) + f"   θ={k*9}〜{(k+3)*9+8}度")
mn=min(c); print(f"  最小 {mn:.4f} が出る方位: {[j*9 for j in range(40) if c[j]<mn+0.004]}")
print(f"  全体の最大 r/τ = {max(max(RAD[i]/TAU[i] for i in range(N) if TAU[i]>0),0):.6f}")

# 遅れの中身
DEG=[len(a) for a in ADJ]; C={2:12,3:8,4:6}
E={}; el=[]
for u in range(N):
    for v in ADJ[u]: E[(u,v)]=len(el); el.append((u,v))
IN=[[] for _ in range(N)]
for k,(u,v) in enumerate(el): IN[v].append(k)
REV=[E[(v,u)] for (u,v) in el]
e0=E[(center,ADJ[center][0])]
psi=[0]*len(el); psi[e0]=1
firstv=[-1]*N; firstv[el[e0][1]]=0
STEPS=40
for t in range(1,STEPS+1):
    new=[0]*len(el)
    for v in range(N):
        S_=0; nz=False
        for k in IN[v]:
            if psi[k]: nz=True
            S_+=psi[k]
        if not nz: continue
        cc=C[DEG[v]]
        for k in IN[v]: new[REV[k]]=cc*S_-12*psi[k]
    psi=new
    for k,(u,v) in enumerate(el):
        if psi[k] and firstv[v]<0: firstv[v]=t
tau0=bfs(el[e0][1],ADJ)
dl=[(v,tau0[v],firstv[v],DEG[v]) for v in range(N) if tau0[v]<=STEPS-1 and firstv[v]!=tau0[v]]
print(f"\n遅れた頂点 {len(dl)} 個: " + ", ".join(f"τ={a} 到達={b} 次数={d}" for _,a,b,d in dl))
