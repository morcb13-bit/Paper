# 複眼の効き：一つの型紙が何箇所に嵌るか、個眼を増やすと1に絞れるか
#   一眼   W(w)      = その語が流れる五芒星の集合
#   二眼   W ∩ (W−b) = 基線 b だけ離した二つの個眼の両方で嵌る位置
#   三眼   さらに c を足す
#   必ず落ちる設定：b=0 なら二眼でも一眼と同じ
#   負の対照：無作為な語で同じことをして絞りが効かないこと
import json, math, random, sys
from collections import defaultdict, Counter
import b13_chain_units as U
phi=(1+5**0.5)/2
N=12
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
        O[(prev,cur)]={0:ws[0]} if len(ws)==1 else ({1:ws[0],-1:ws[1]} if len(ws)==2 else {1:ws[0],0:ws[1],-1:ws[2]})
ST=[]; CEN=[]
for area,cyc in U.gaps(cells):
    if abs(area-2.9389)<0.01:
        c=(0,0,0,0)
        for p in cyc: c=U.zadd(c,p)
        if any(len(adj[p])<2 for p in cyc): continue
        CEN.append(c); ST.append([(p,cyc[(i+1)%10]) for i,p in enumerate(cyc)])
n=len(CEN); print("五芒星 %d 個" % n)
key={c:i for i,c in enumerate(CEN)}
def flows(w,starts):
    for a,b in starts:
        prev,cur=a,b; ok=True
        for x in w:
            o=O[(prev,cur)]
            if x not in o: ok=False; break
            prev,cur=cur,o[x]
        if ok: return True
    return False
rng=random.Random(20260903)
words=[]
while len(words)<600:
    w=tuple(rng.choice((-1,0,1)) for _ in range(N))
    words.append(w)
Ws=[]
for w in words:
    S={i for i in range(n) if flows(w,ST[i])}
    if S: Ws.append((w,S))
print("受理された語 %d / %d" % (len(Ws),len(words)))
# 基線の候補：五芒星どうしの平行移動ベクトル（担体の中で多くの星に効くもの）
cand=[]
for c in CEN[:40]:
    if all(x%10==0 for x in c):
        t=tuple(x//10 for x in c); d=math.hypot(*U.xy(t))
        if 0<d<26: cand.append((round(d,4),t))
cand=sorted(set(cand))
def shift(S,t):
    out=set()
    for i in S:
        j=key.get(U.zadd(CEN[i],tuple(10*x for x in t)))
        if j is not None: out.add(j)
    return out
print()
print("基線 b      一眼|W|中央/最小   二眼|W|中央/最小   1点に絞れた語")
sz1=sorted(len(S) for _,S in Ws)
print("  （b=0）      %4d / %2d          %4d / %2d          %d"
      % (sz1[len(sz1)//2],sz1[0],sz1[len(sz1)//2],sz1[0],sum(1 for _,S in Ws if len(S)==1)))
for d,t in cand[:6]:
    two=[S & shift(S,t) for _,S in Ws]
    two=[x for x in two if x]
    if not two: continue
    s2=sorted(len(x) for x in two)
    print("  %8.4f    %4d / %2d          %4d / %2d          %d"
          % (d,sz1[len(sz1)//2],sz1[0],s2[len(s2)//2],s2[0],sum(1 for x in two if len(x)==1)))
