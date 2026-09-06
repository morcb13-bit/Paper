"""検定LC3b 遅れは担体か、初期条件の鏡映対称か
     OK（初期条件を非対称な位置に移すと S=0 と遅れが消える）→ 鏡映の産物。担体の性質ではない
     NG（移しても残る）→ 担体の局所配置が規則に干渉している
   必ず落ちる設定：元の中心では 8 件・5点が再現すること
   検定LC4 三角格子の半径依存"""
import json, math, collections
G=json.load(open('carrier_1245_graph.json'))
XY=[complex(x,y) for x,y in G['xy']]; ADJ=[list(a) for a in G['adj']]; N=len(XY)
DEG=[len(a) for a in ADJ]; C={2:12,3:8,4:6}
def bfs(s,adj):
    d=[-1]*len(adj); d[s]=0; q=collections.deque([s])
    while q:
        u=q.popleft()
        for v in adj[u]:
            if d[v]<0: d[v]=d[u]+1; q.append(v)
    return d
E={}; el=[]
for u in range(N):
    for v in ADJ[u]: E[(u,v)]=len(el); el.append((u,v))
IN=[[] for _ in range(N)]
for k,(u,v) in enumerate(el): IN[v].append(k)
REV=[E[(v,u)] for (u,v) in el]

def trial(src, nb_idx, steps, label):
    e0=E[(src,ADJ[src][nb_idx])]; start=el[e0][1]
    TAU=bfs(start,ADJ)
    psi=[0]*len(el); psi[e0]=1
    firstv=[-1]*N; firstv[start]=0
    sz=[]
    for t in range(1,steps+1):
        new=[0]*len(el)
        for v in range(N):
            S_=0; nz=False
            for k in IN[v]:
                if psi[k]: nz=True
                S_+=psi[k]
            if not nz: continue
            if S_==0: sz.append((t-1,v))
            cc=C[DEG[v]]
            for k in IN[v]: new[REV[k]]=cc*S_-12*psi[k]
        psi=new
        for k,(u,v) in enumerate(el):
            if psi[k] and firstv[v]<0: firstv[v]=t
    late=[v for v in range(N) if 0<TAU[v]<=steps-1 and firstv[v]!=TAU[v]]
    reach=sum(1 for v in range(N) if 0<TAU[v]<=steps-1)
    ang=sorted({round(math.degrees(math.atan2(XY[v].imag,XY[v].real))%360,2) for _,v in sz})
    print(f"{label}: 到達{reach} 遅れ{len(late)} S=0が{len(sz)}件 "
          f"S=0の方位={ang}")
    return late, sz

c0=min(range(N),key=lambda i:abs(XY[i]))
trial(c0,0,32,"中心 出発（元の設定）")
trial(c0,1,32,"中心 出発・辺を変える")
for tgt in (complex(30,0), complex(0,30), complex(21,21), complex(-17,9)):
    s=min(range(N),key=lambda i:abs(XY[i]-tgt))
    trial(s,0,32,f"出発を ({tgt.real:.0f},{tgt.imag:.0f}) 付近へ")

# ── LC4 三角格子の半径依存
def ripple(tau,rad,xy,nb,tl):
    b=[int((math.atan2(p.imag,p.real)%(2*math.pi))/(2*math.pi)*nb)%nb for p in xy]
    out=[]
    for t in tl:
        R=[0.0]*nb
        for i in range(len(tau)):
            if 0<=tau[i]<=t and rad[i]>R[b[i]]: R[b[i]]=rad[i]
        me=sum(R)/nb; out.append((t,me,max(R)-min(R),(max(R)-min(R))/me))
    return out
T=130; tri={}; txy=[]; tadj=[]
def tid(a,b):
    if (a,b) not in tri:
        tri[(a,b)]=len(txy); txy.append(complex(a+b*0.5,b*math.sqrt(3)/2)); tadj.append([])
    return tri[(a,b)]
for a in range(-T,T+1):
    for b in range(-T,T+1):
        if abs(complex(a+b*0.5,b*math.sqrt(3)/2))<=T: tid(a,b)
for (a,b),i in list(tri.items()):
    for da,db in ((1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)):
        if (a+da,b+db) in tri: tadj[i].append(tri[(a+da,b+db)])
tt=bfs(tri[(0,0)],tadj); tr=[abs(p) for p in txy]
print("\n検定LC4 三角格子（6近傍・区間20）")
print(f"{'t':>4} {'平均R':>8} {'絶対幅':>8} {'相対幅':>8}")
for t,me,ab,rl in ripple(tt,tr,txy,20,[20,30,40,50,60,70,80,90,100]):
    print(f"{t:>4} {me:8.2f} {ab:8.3f} {rl:8.4f}")
