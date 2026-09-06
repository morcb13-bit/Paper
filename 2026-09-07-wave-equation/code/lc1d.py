import json, math, collections
G=json.load(open('carrier_1245_graph.json'))
XY=[complex(x,y) for x,y in G['xy']]; ADJ=[list(a) for a in G['adj']]; N=len(XY)
def bfs(s,adj):
    d=[-1]*len(adj); d[s]=0; q=collections.deque([s])
    while q:
        u=q.popleft()
        for v in adj[u]:
            if d[v]<0: d[v]=d[u]+1; q.append(v)
    return d
RAD=[abs(p) for p in XY]
c0=min(range(N),key=lambda i:RAD[i])
TAU=bfs(c0,ADJ)
for lo in (20,40,60,80):
    m=max(RAD[i]/TAU[i] for i in range(N) if TAU[i]>=lo)
    print(f"τ≥{lo:>3} の最大 r/τ = {m:.6f}")
# 中心を変えて c(θ) の遅い偏りが動くか
def cprof(src,nb=20,rmin=60,rmax=100):
    t=bfs(src,ADJ); o=XY[src]
    c=[0.0]*nb
    for i in range(N):
        d=XY[i]-o; r=abs(d)
        if t[i]>0 and rmin<=r<=rmax:
            b=int((math.atan2(d.imag,d.real)%(2*math.pi))/(2*math.pi)*nb)%nb
            c[b]=max(c[b],r/t[i])
    return c
for src,lab in ((c0,"中心"),(min(range(N),key=lambda i:abs(XY[i]-complex(30,0))),"x=30の頂点"),
                (min(range(N),key=lambda i:abs(XY[i]-complex(0,30))),"y=30の頂点")):
    c=cprof(src)
    k=min(range(20),key=lambda j:c[j])
    print(f"{lab:>10}: c={min(c):.4f}〜{max(c):.4f} ばらつき{(max(c)-min(c))/(sum(c)/20):.4f} 最小方位 {k*18}度")
