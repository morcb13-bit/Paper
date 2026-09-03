# 検定TR1：ルートを座標の列ではなく「語」で持てるか
#   一つの語を繰り返して歩かせ、どこまで進むか・閉じるか・何枚を覆うかを測る
#   OK なら：短い語の繰り返しで広い範囲を覆う経路が書ける（＝ルートが語で持てる）
#   NG なら：どの語も数回で止まるか、すぐ閉じて同じ所を回る
#   必ず落ちる設定：全部0の語は一歩も進まない
#   負の対照：桁と出口の対応を無作為に付け替えたとき、覆いが同じなら幾何は効いていない
import json, math, random, sys
from collections import defaultdict, Counter
import b13_chain_units as U
cells={tuple(int(x) for x in k.split(",")):a for k,a in json.load(open("carrier_1245.json"))["cells"].items()}
adj=defaultdict(set)
for q,a in cells.items():
    vs=[U.zadd(q,U.zt(a+2*i)) for i in range(5)]
    for i in range(5):
        u,w=vs[i],vs[(i+1)%5]; adj[u].add(w); adj[w].add(u)
XY={p:U.xy(p) for p in adj}
def build(shuffle,seed=3):
    rng=random.Random(seed); O={}
    for cur in adj:
        for prev in adj[cur]:
            outs=[w for w in adj[cur] if w!=prev]
            if not outs: O[(prev,cur)]={}; continue
            back=math.atan2(XY[prev][1]-XY[cur][1],XY[prev][0]-XY[cur][0])
            ws=[w for _,w in sorted(((math.atan2(XY[w][1]-XY[cur][1],XY[w][0]-XY[cur][0])-back)%(2*math.pi),w) for w in outs)]
            if shuffle: rng.shuffle(ws)
            O[(prev,cur)]={0:ws[0]} if len(ws)==1 else ({1:ws[0],-1:ws[1]} if len(ws)==2 else {1:ws[0],0:ws[1],-1:ws[2]})
    return O
CEN=[];RIM=[]
for area,cyc in U.gaps(cells):
    if abs(area-2.9389)<0.01:
        c=(0,0,0,0)
        for p in cyc: c=U.zadd(c,p)
        if any(len(adj[p])<2 for p in cyc): continue
        CEN.append((U.xy(c)[0]/10,U.xy(c)[1]/10)); RIM.append(cyc)
start=[(RIM[0][i],RIM[0][(i+1)%10]) for i in range(10)]
L=12; REP=40
def run(O,words):
    res=[]
    for w in words:
        for st in start:
            prev,cur=st; seen={cur}; edges={st}; steps=0; closed=None
            ok=True
            for r in range(REP):
                for x in w:
                    o=O[(prev,cur)]
                    if x not in o: ok=False; break
                    prev,cur=cur,o[x]; steps+=1; seen.add(cur)
                    if (prev,cur) in edges: closed=steps; break
                    edges.add((prev,cur))
                if not ok or closed: break
            if steps>0:
                res.append((steps,len(seen),closed is not None,
                            max(math.hypot(*XY[p]) for p in seen)))
    return res
rng=random.Random(7)
words=[tuple(rng.choice((-1,0,1)) for _ in range(L)) for _ in range(3000)]
print("必ず落ちる設定：全部0の語  →  %s" % ("進んだ" if run(build(False),[(0,)*L]) else "一歩も進まない"))
for tag,sh in (("本番",False),("負の対照（桁の付け替え）",True)):
    r=run(build(sh),words)
    if not r: print(tag,"進む語なし"); continue
    st=sorted(x[0] for x in r); cov=sorted(x[1] for x in r)
    cl=sum(1 for x in r if x[2]); far=sorted(x[3] for x in r)
    full=sum(1 for x in r if x[0]==L*REP)
    print("\n== %s ==" % tag)
    print("  歩けた出発 %d 通り（語3000本×10角のうち）" % len(r))
    print("  歩数 中央 %d / 最大 %d（上限 %d）  最後まで歩いた %d 通り" % (st[len(st)//2],st[-1],L*REP,full))
    print("  覆った頂点 中央 %d / 最大 %d" % (cov[len(cov)//2],cov[-1]))
    print("  途中で閉じた（同じ有向辺に戻った） %d / %d" % (cl,len(r)))
    print("  到達半径 中央 %.1f / 最大 %.1f" % (far[len(far)//2],far[-1]))
