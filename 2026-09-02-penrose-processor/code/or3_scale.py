# 検定OR3-S  担体を大きくすると |W| と W の種類数はどう伸びるか
#   盤面は 1245環に固定し、原点に使う五芒星の半径だけを広げる（境界の効き方を一定に保つ）
#   A) 受理率 |W|/出発数     B) 相異なる W の種類数 N_W    C) |W| の分布（二山が残るか）
import json, math, random, sys
from collections import defaultdict, Counter
import b13_chain_units as U
N=12; M=int(sys.argv[1]) if len(sys.argv)>1 else 2000
cells={tuple(int(x) for x in k.split(",")):a for k,a in json.load(open("carrier_1245.json"))["cells"].items()}
adj=defaultdict(set)
for q,a in cells.items():
    vs=[U.zadd(q,U.zt(a+2*i)) for i in range(5)]
    for i in range(5):
        u,w=vs[i],vs[(i+1)%5]; adj[u].add(w); adj[w].add(u)
XY={p:U.xy(p) for p in adj}
print("盤面：頂点 %d／次数 %s" % (len(adj), dict(sorted(Counter(len(v) for v in adj.values()).items()))))
O={}
for cur in adj:
    for prev in adj[cur]:
        outs=[w for w in adj[cur] if w!=prev]
        if not outs: O[(prev,cur)]={}; continue
        back=math.atan2(XY[prev][1]-XY[cur][1],XY[prev][0]-XY[cur][0])
        ws=[w for _,w in sorted(((math.atan2(XY[w][1]-XY[cur][1],XY[w][0]-XY[cur][0])-back)%(2*math.pi),w) for w in outs)]
        if len(ws)==1: O[(prev,cur)]={0:ws[0]}
        elif len(ws)==2: O[(prev,cur)]={1:ws[0],-1:ws[1]}
        else: O[(prev,cur)]={1:ws[0],0:ws[1],-1:ws[2]}
stars=[cyc for area,cyc in U.gaps(cells) if abs(area-2.9389)<0.01]
SS=[]
for cyc in stars:
    c10=(0,0,0,0)
    for p in cyc: c10=U.zadd(c10,p)
    cx,cy=U.xy(c10)[0]/10,U.xy(c10)[1]/10
    P=[U.xy(p) for p in cyc]
    ar=sum(P[i][0]*P[(i+1)%10][1]-P[(i+1)%10][0]*P[i][1] for i in range(10))
    if ar<0: cyc=cyc[::-1]
    if any(len(adj[p])<2 for p in cyc): continue
    SS.append((round(math.hypot(cx,cy),4),[(p,cyc[(i+1)%10]) for i,p in enumerate(cyc)]))
SS.sort()
print("五芒星 %d 個（最大半径 %.2f）" % (len(SS), SS[-1][0]))
rng=random.Random(20260903)
words=[tuple(rng.choice((-1,0,1)) for _ in range(N)) for _ in range(M)]
CUTS=[(21,"60環相当"),(61,""),(106,""),(226,""),(441,"")]
ALL=[st for R,sts in SS for st in sts]
def accepts(starts,w):
    W=[]
    for i,(a,b) in enumerate(starts):
        prev,cur=a,b; ok=True
        for x in w:
            o=O[(prev,cur)]
            if x not in o: ok=False; break
            prev,cur=cur,o[x]
        if ok: W.append(i)
    return W
print("\n 星数  出発数   受理語  受理率  <|W|>  |W|/出発  相異なるW  W中央  W最大")
res={}
for k,_ in CUTS:
    if k>len(SS): continue
    starts=[st for R,sts in SS[:k] for st in sts]
    Ws=[frozenset(accepts(starts,w)) for w in words]
    nz=[W for W in Ws if W]
    sz=sorted(len(W) for W in nz)
    res[k]=sz
    print("%5d %7d %7d %6.2f%% %6.1f %8.4f %9d %6d %6d"
          % (k,len(starts),len(nz),100*len(nz)/M,sum(sz)/len(sz),
             sum(sz)/len(sz)/len(starts),len(set(nz)),sz[len(sz)//2],sz[-1]))
print("\n■ |W| の分布（出発数で割った受理率の帯ごと）")
for k in sorted(res):
    n=k*10; h=Counter(min(int(20*s/n),19) for s in res[k])
    print("  星%3d  %s" % (k, " ".join("%d%%:%d"%(5*b,h[b]) for b in sorted(h))))
