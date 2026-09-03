# 検定OR3-Q  (a) 増加率 1.695 は何で決まるか  (b) 型の中に何個の状態が入るか
#   (a) 代表の本数 r を 1,2,3,5,8,13,21,42 と増やし、N_M(n) の増加率が
#       r とともに動くか止まるかを見る。止まるなら 1.695 は少数の出発で決まる数、
#       動くなら出発の集まり方が決めている数
#   (b) n を固定し、W ごとに何個の M が入るかの分布を出す
import json, math, sys
from collections import defaultdict, Counter
import b13_chain_units as U
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 16
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
allst=[st for R,sts in SS[:21] for st in sts]
pos={st:i for i,st in enumerate(allst)}; seen=set(); reps=[]
for i,st in enumerate(allst):
    if i in seen: continue
    orb=[]; s=st
    for _ in range(5): orb.append(pos[s]); s=(U.zrot(s[0],2),U.zrot(s[1],2))
    seen|=set(orb); reps.append(i)
print("C5 の代表 %d 本" % len(reps))

def run(rr,nmax,keepWM=None):
    root=tuple((k,EI[allst[i]]) for k,i in enumerate(rr)); NR=len(rr)
    A=[0]*(nmax+1); WS=[set() for _ in range(nmax+1)]; MS=[set() for _ in range(nmax+1)]
    pairs=[] if keepWM else None
    def dfs(al,d):
        if d:
            A[d]+=1; WS[d].add(hash(tuple(i for i,_ in al))); MS[d].add(hash(al))
            if keepWM==d: pairs.append((hash(tuple(i for i,_ in al)),hash(al)))
        if d==nmax: return
        for x in (0,1,2):
            n2=tuple((i,t) for i,e in al if (t:=T[e][x])>=0)
            if n2: dfs(n2,d+1)
    sys.setrecursionlimit(nmax+50); dfs(root,0)
    return A,[len(s) for s in WS],[len(s) for s in MS],pairs

print("\n■ (a) 代表の本数 r を変えたとき N_M(n) の増加率")
print("   r   N_M(%d)   最後の3つの増加率" % NMAX)
for r in (1,2,3,5,8,13,21,42):
    rr=reps[:r]
    A,W,M,_=run(rr,NMAX)
    rat=[M[n+1]/M[n] for n in range(NMAX-3,NMAX)]
    print("%4d %9d   %s" % (r,M[NMAX]," ".join("%.3f"%x for x in rat)))

print("\n■ (b) 一つの W の中に入る M の個数（代表42本・長さ%d）" % NMAX)
A,W,M,pairs=run(reps,NMAX,keepWM=NMAX)
g=defaultdict(set)
for w,m in pairs: g[w].add(m)
c=sorted(len(v) for v in g.values())
print("   W %d 種／M %d 種／1つの W に入る M：最小 %d 中央 %d 平均 %.1f 最大 %d"
      % (len(g),len(M and set()) if False else M[NMAX],c[0],c[len(c)//2],sum(c)/len(c),c[-1]))
h=Counter(1 if x==1 else (2 if x<=4 else (3 if x<=16 else (4 if x<=64 else 5))) for x in c)
lab={1:"1個",2:"2〜4",3:"5〜16",4:"17〜64",5:"65以上"}
print("   分布 %s" % {lab[k]:h[k] for k in sorted(h)})
print("   M が1個しかない W（型が状態を隠していない）… %d / %d（%.1f%%）"
      % (h[1],len(g),100*h[1]/len(g)))
