# or3_state.py と同じ測定を、C5 の代表だけで走らせる（5倍速）
#   前提：担体は 72° 回転で自分に写り、桁の割り当ても折れ角から決まるので同変。
#         よって生存集合は必ず軌道の和になり、代表42本で W も M も決まる。
#   この前提自体を n=12 の既知の値（A=54045, N_W=1116, N_M=15962）で検算する。
import json, math, sys
from collections import defaultdict
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
    for prev in adj[cur]: EI[(prev,cur)]=len(E); E.append((prev,cur))
T=[[-1,-1,-1] for _ in E]
for (prev,cur),ei in EI.items():
    outs=[w for w in adj[cur] if w!=prev]
    if not outs: continue
    back=math.atan2(XY[prev][1]-XY[cur][1],XY[prev][0]-XY[cur][0])
    ws=[w for _,w in sorted(((math.atan2(XY[w][1]-XY[cur][1],XY[w][0]-XY[cur][0])-back)%(2*math.pi),w) for w in outs)]
    mp={0:ws[0]} if len(ws)==1 else ({1:ws[0],-1:ws[1]} if len(ws)==2 else {1:ws[0],0:ws[1],-1:ws[2]})
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
allst=[st for R,sts in SS[:K] for st in sts]
# 72°回転で軌道に分ける
def rot(p): return U.zrot(p,2)
pos={st:i for i,st in enumerate(allst)}
seen=set(); reps=[]; osz=[]
for i,st in enumerate(allst):
    if i in seen: continue
    orb=[]; s=st
    for _ in range(5):
        j=pos.get(s)
        if j is None: orb=None; break
        orb.append(j); s=(rot(s[0]),rot(s[1]))
    if orb is None: print("回転で閉じない出発がある"); sys.exit(1)
    seen|=set(orb); reps.append(i); osz.append(len(set(orb)))
print("出発 %d 本 → C5 の軌道 %d 個（軌道の大きさ %s）" % (len(allst),len(reps),sorted(set(osz))))
root=tuple((k,EI[allst[i]]) for k,i in enumerate(reps)); NR=len(reps)
A=[0]*(NMAX+1); WS=[set() for _ in range(NMAX+1)]; MS=[set() for _ in range(NMAX+1)]
SAME=[0]*(NMAX+1); FULL=[0]*(NMAX+1); TOT=[0]*(NMAX+1)
def dfs(al,d):
    if d:
        A[d]+=1; TOT[d]+=len(al)
        WS[d].add(hash(tuple(i for i,_ in al))); MS[d].add(hash(al))
        if len(al)==NR: FULL[d]+=1
    if d==NMAX: return
    for x in (0,1,2):
        n2=tuple((i,t) for i,e in al if (t:=T[e][x])>=0)
        if n2:
            if len(n2)==len(al): SAME[d+1]+=1
            dfs(n2,d+1)
sys.setrecursionlimit(NMAX+50)
dfs(root,0)
print("\n  n      A(n)      N_W(n)     N_M(n)  N_M/N_W   A/N_W   縮まらない  100%語  <|W|>/出発")
for n in range(1,NMAX+1):
    print("%3d %10d %11d %10d %8.4f %8.1f %10.4f %7d %9.4f"
          % (n,A[n],len(WS[n]),len(MS[n]),len(MS[n])/len(WS[n]),A[n]/len(WS[n]),
             SAME[n]/A[n],FULL[n],TOT[n]/A[n]/NR))
