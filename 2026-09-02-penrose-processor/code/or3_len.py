# 検定OR3-L  語長 n に対して N_W(n) を測る
#   出すもの：n=1..NMAX について 受理語数 A(n)／相異なる W の種類数 N_W(n)／
#             W が全出発になる語（100%語）の本数とその共通接頭辞
#   ケースA：N_W が急増 → 語長が照合パターンを増やしている
#   ケースB：N_W が頭打ち → 盤面の規則が許す照合型が有限個
#   ケースC：A だけ増えて N_W が増えない → 既存の型を持つ語が増えているだけ
#   必ず落ちる設定：全部0の語がどの n でも受理されないこと
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
starts=[st for R,sts in SS[:K] for st in sts]; NS=len(starts)
print("五芒星 %d 個／出発 %d 本" % (K,NS))
print("\n  n    全語数     受理語 A(n)    A/3^n     N_W(n)   A/N_W   100%%語  <|W|>/出発")
front=[([(i,st) for i,st in enumerate(starts)],())]
for n in range(1,NMAX+1):
    nx=[]
    for al,w in front:
        for x in (-1,0,1):
            n2=[(i,(e[1],O[e][x])) for i,e in al if x in O[e]]
            if n2: nx.append((n2,w+(x,)))
    front=nx
    sigs=set(); full=[]; tot=0
    for al,w in front:
        sigs.add(tuple(i for i,_ in al)); tot+=len(al)
        if len(al)==NS: full.append(w)
    A=len(front)
    print("%3d %10d %12d %9.5f %10d %7.1f %7d %9.4f"
          % (n,3**n,A,A/3**n,len(sigs),A/len(sigs),len(full),tot/A/NS))
    if full:
        p=full[0]
        for w in full:
            p=p[:next((i for i,(a,b) in enumerate(zip(p,w)) if a!=b), len(p))]
        print("      100%%語の共通接頭辞 %s（%d桁）" % (",".join("%+d"%x for x in p) if p else "なし", len(p)))
