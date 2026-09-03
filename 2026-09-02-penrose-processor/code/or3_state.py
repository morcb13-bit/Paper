# 検定OR3-S2  w ↦ W(w) は「状態」か（深さ優先・全接頭辞を保持しない）
#   測る量（n ごと）
#     A(n)   受理語
#     N_W(n) 相異なる受理集合 W
#     N_M(n) 相異なる「生存写像」M（どの出発が いま どの有向辺にいるか）
#   N_M = N_W なら W が M を決める＝W が状態。W 空間上に有限機械が書ける
#   N_M > N_W なら W は M の粗い像にすぎず、W から次の W は決まらない
#   包含 W(wd) ⊆ W(w) は定義から自動（接頭辞が流れないと語は流れない）。測らない。
#   代わりに測るもの：一桁足しても W が縮まらない割合（W(wd)=W(w) の頻度）
import json, math, sys
from collections import defaultdict, Counter
import b13_chain_units as U
K=int(sys.argv[1]); NMAX=int(sys.argv[2])
cells={tuple(int(x) for x in k.split(",")):a for k,a in json.load(open("carrier_1245.json"))["cells"].items()}
adj=defaultdict(set)
for q,a in cells.items():
    vs=[U.zadd(q,U.zt(a+2*i)) for i in range(5)]
    for i in range(5):
        u,w=vs[i],vs[(i+1)%5]; adj[u].add(w); adj[w].add(u)
XY={p:U.xy(p) for p in adj}
E=[]; EI={}
for cur in adj:
    for prev in adj[cur]:
        EI[(prev,cur)]=len(E); E.append((prev,cur))
T=[[-1,-1,-1] for _ in E]
for (prev,cur),ei in EI.items():
    outs=[w for w in adj[cur] if w!=prev]
    if not outs: continue
    back=math.atan2(XY[prev][1]-XY[cur][1],XY[prev][0]-XY[cur][0])
    ws=[w for _,w in sorted(((math.atan2(XY[w][1]-XY[cur][1],XY[w][0]-XY[cur][0])-back)%(2*math.pi),w) for w in outs)]
    if len(ws)==1: mp={1:ws[0]} if False else {0:ws[0]}
    elif len(ws)==2: mp={1:ws[0],-1:ws[1]}
    else: mp={1:ws[0],0:ws[1],-1:ws[2]}
    for x,w in mp.items(): T[ei][x+1]=EI[(cur,w)]
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
starts=[EI[st] for R,sts in SS[:K] for st in sts]; NS=len(starts)
print("五芒星 %d 個／出発 %d 本／有向辺 %d" % (K,NS,len(E)))
root=tuple((i,e) for i,e in enumerate(starts))
A=[0]*(NMAX+1); WS=[set() for _ in range(NMAX+1)]; MS=[set() for _ in range(NMAX+1)]
SAME=[0]*(NMAX+1); FULL=[0]*(NMAX+1); TOT=[0]*(NMAX+1)
def dfs(al,d):
    if d:
        A[d]+=1; TOT[d]+=len(al)
        WS[d].add(hash(tuple(i for i,_ in al)))
        MS[d].add(hash(al))
        if len(al)==NS: FULL[d]+=1
    if d==NMAX: return
    for x in (0,1,2):
        n2=tuple((i,t) for i,e in al if (t:=T[e][x])>=0)
        if n2:
            if len(n2)==len(al): SAME[d+1]+=1
            dfs(n2,d+1)
sys.setrecursionlimit(NMAX+50)
dfs(root,0)
print("\n  n      A(n)      N_W(n)     N_M(n)  N_M/N_W   A/N_W   縮まらない子語  100%語  <|W|>/出発")
for n in range(1,NMAX+1):
    print("%3d %10d %11d %10d %8.4f %7.1f %12.4f %7d %9.4f"
          % (n,A[n],len(WS[n]),len(MS[n]),len(MS[n])/len(WS[n]),A[n]/len(WS[n]),
             SAME[n]/A[n],FULL[n],TOT[n]/A[n]/NS))
