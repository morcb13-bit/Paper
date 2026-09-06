"""LC1b：前線の速さ c(θ) = max{ r(v)/τ(v) } を方位区間ごとに出す。
   LC2b：頂点の側で、振幅が最初に非零になる歩数と歩数 τ を比べる。"""
import json, math, collections

G = json.load(open('carrier_1245_graph.json'))
XY = [complex(x, y) for x, y in G['xy']]
ADJ = [list(a) for a in G['adj']]
N = len(XY)
center = min(range(N), key=lambda i: abs(XY[i]))

def bfs(src, adj):
    d = [-1] * len(adj); d[src] = 0
    q = collections.deque([src])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if d[v] < 0:
                d[v] = d[u] + 1; q.append(v)
    return d

TAU = bfs(center, ADJ)
RAD = [abs(p) for p in XY]

def speed_by_bin(nb, rmin, rmax, tau, rad, xy):
    c = [0.0] * nb
    cnt = [0] * nb
    for i in range(len(tau)):
        if tau[i] <= 0 or not (rmin <= rad[i] <= rmax):
            continue
        b = int((math.atan2(xy[i].imag, xy[i].real) % (2*math.pi)) / (2*math.pi) * nb) % nb
        cnt[b] += 1
        c[b] = max(c[b], rad[i] / tau[i])
    return c, cnt

for nb in (10, 20, 40):
    c, cnt = speed_by_bin(nb, 40, 100, TAU, RAD, XY)
    mx, mn, me = max(c), min(c), sum(c)/nb
    print(f"区間{nb:>3}: 前線の速さ c(θ) = {mn:.5f}〜{mx:.5f} 平均{me:.5f} "
          f"ばらつき {(mx-mn)/me:.5f}  各区間の標本数 {min(cnt)}〜{max(cnt)}")

# 半径帯ごとに（縁の影響を見る）
print("\n半径帯ごと（区間20）")
for lo, hi in ((20,40),(40,60),(60,80),(80,100),(100,107)):
    c, cnt = speed_by_bin(20, lo, hi, TAU, RAD, XY)
    mx, mn, me = max(c), min(c), sum(c)/20
    print(f"  r={lo:>3}〜{hi:>3}: c = {mn:.5f}〜{mx:.5f} ばらつき {(mx-mn)/me:.5f}")

# 負の対照：正方格子
S = 120
def sid(x,y): return (x+S)*(2*S+1)+(y+S)
sadj=[[] for _ in range((2*S+1)**2)]; sxy=[0j]*len(sadj)
for x in range(-S,S+1):
    for y in range(-S,S+1):
        i=sid(x,y); sxy[i]=complex(x,y)
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
            if -S<=x+dx<=S and -S<=y+dy<=S: sadj[i].append(sid(x+dx,y+dy))
stau = bfs(sid(0,0), sadj); srad=[abs(p) for p in sxy]
c, cnt = speed_by_bin(20, 40, 100, stau, srad, sxy)
mx,mn,me = max(c),min(c),sum(c)/20
print(f"\n負の対照 正方格子: c = {mn:.5f}〜{mx:.5f} ばらつき {(mx-mn)/me:.5f}（理論 1−1/√2 の形）")

# 三角格子（6近傍）も
T=120; tri={}; tadj=[]; txy=[]
def tid(a,b):
    if (a,b) not in tri:
        tri[(a,b)]=len(txy); txy.append(complex(a+b*0.5, b*math.sqrt(3)/2)); tadj.append([])
    return tri[(a,b)]
for a in range(-T,T+1):
    for b in range(-T,T+1):
        if abs(complex(a+b*0.5,b*math.sqrt(3)/2))<=T: tid(a,b)
for (a,b),i in list(tri.items()):
    for da,db in ((1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)):
        if (a+da,b+db) in tri: tadj[i].append(tri[(a+da,b+db)])
ttau = bfs(tri[(0,0)], tadj); trad=[abs(p) for p in txy]
c,cnt = speed_by_bin(20, 40, 100, ttau, trad, txy)
mx,mn,me=max(c),min(c),sum(c)/20
print(f"負の対照 三角格子: c = {mn:.5f}〜{mx:.5f} ばらつき {(mx-mn)/me:.5f}")

# ── LC2b 頂点の側
DEG=[len(a) for a in ADJ]; C={2:12,3:8,4:6}
E={}; elist=[]
for u in range(N):
    for v in ADJ[u]:
        E[(u,v)]=len(elist); elist.append((u,v))
IN=[[] for _ in range(N)]
for k,(u,v) in enumerate(elist): IN[v].append(k)
REV=[E[(v,u)] for (u,v) in elist]

def run(steps, e0):
    psi=[0]*len(elist); psi[e0]=1
    firstv=[-1]*N; firstv[elist[e0][1]]=0
    ss=[1]
    for t in range(1,steps+1):
        new=[0]*len(elist)
        for v in range(N):
            S_=0; nz=False
            for k in IN[v]:
                if psi[k]: nz=True
                S_+=psi[k]
            if not nz: continue
            c=C[DEG[v]]
            for k in IN[v]: new[REV[k]]=c*S_-12*psi[k]
        psi=new
        for k,(u,v) in enumerate(elist):
            if psi[k] and firstv[v]<0: firstv[v]=t
        ss.append(sum(x*x for x in psi))
    return firstv, ss

e0=E[(center,ADJ[center][0])]
STEPS=40
firstv, ss = run(STEPS,e0)
tau0 = bfs(elist[e0][1], ADJ)
ag=late=never=0
for v in range(N):
    if tau0[v]>STEPS-1: continue
    if firstv[v]==tau0[v]: ag+=1
    elif firstv[v]<0: never+=1
    else: late+=1
print(f"\nLC2b 頂点の側（{STEPS}歩）: 歩数と一致 {ag} / 遅れ {late} / 未到達 {never}")
print(f"  二乗和 毎歩144倍: {all(ss[i+1]==144*ss[i] for i in range(len(ss)-1))}")
