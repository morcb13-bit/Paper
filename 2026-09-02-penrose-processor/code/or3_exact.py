# 長さ12の全語 3^12 について W(w) を厳密に数える（生きている接頭辞だけ伸ばす）
import json, math, sys
from collections import defaultdict, Counter
import b13_chain_units as U
N=12
K=int(sys.argv[1]) if len(sys.argv)>1 else 21
cells={tuple(int(x) for x in k.split(",")):a for k,a in json.load(open("carrier_1245.json"))["cells"].items()}
adj=defaultdict(set)
for q,a in cells.items():
    vs=[U.zadd(q,U.zt(a+2*i)) for i in range(5)]
    for i in range(5):
        u,w=vs[i],vs[(i+1)%5]; adj[u].add(w); adj[w].add(u)
XY={p:U.xy(p) for p in adj}
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
starts=[st for R,sts in SS[:K] for st in sts]
idx={st:i for i,st in enumerate(starts)}
print("五芒星 %d 個／出発 %d 本" % (K,len(starts)))
# 生きている接頭辞だけ伸ばす。状態＝(出発の番号, 現在の有向辺)
live=[(frozenset((i,st) for i,st in enumerate(starts)),)]
cur=[[(i,st) for i,st in enumerate(starts)]]
for step in range(N):
    nx=[]
    for al in cur:
        for x in (-1,0,1):
            n2=[(i,O[e][x]and (e[1],O[e][x])) for i,e in al if x in O[e]]
            n2=[(i,(e[1],O[e][x])) for i,e in al if x in O[e]]
            if n2: nx.append(n2)
    cur=nx
    print("  桁%2d：生きている語 %d 本" % (step+1,len(cur)))
Ws=[frozenset(i for i,_ in al) for al in cur]
sz=sorted(len(W) for W in Ws)
print("\n長さ%d の全語 %d 本のうち 受理された語 %d 本（%.3f%%）" % (N,3**N,len(Ws),100*len(Ws)/3**N))
print("|W| 最小 %d／中央 %d／最大 %d／平均 %.1f（出発%d本）" % (sz[0],sz[len(sz)//2],sz[-1],sum(sz)/len(sz),len(starts)))
print("相異なる W … %d 種" % len(set(Ws)))
h=Counter(int(20*s/len(starts)) for s in sz)
print("|W|/出発 の分布 %s" % {"%d%%"%(5*b):h[b] for b in sorted(h)})
json.dump({"K":K,"nstart":len(starts),"acc":len(Ws),"NW":len(set(Ws)),
           "hist":{str(b):h[b] for b in h}},open("or3_exact_%d.json"%K,"w"))
