# 検定OR3  w ↦ W(w)（語がどの出発から流れるか＝受理集合）
#
#   出すもの：長さNの語 w に対し W(w) ⊂ 210本の出発（21星×10角）
#   OK なら：W(w) が殻や向きの中で閉じる／語ごとに W が違う（受理集合が語を分ける）
#            ＝盤面は語の受理を選別する装置
#   NG なら：W(w) が 42 の C₅ 軌道の和にしかならず、しかも配置換えグラフと
#            同じ分布になる ＝ 受理は次数の話で、幾何は効いていない
#   必ず落ちる設定：全部 0 の語は W=∅
#   負の対照：同じ次数列の配置換えグラフ（幾何なし）で |W| の分布と W の種類数を比べる
import json, math, random, sys
from collections import defaultdict, Counter
import b13_chain_units as U
N = int(sys.argv[1]) if len(sys.argv)>1 else 12
M = int(sys.argv[2]) if len(sys.argv)>2 else 4000
d=json.load(open("rings_integer.json"))
cells={tuple(eval(k)):v for k,v in d["cells"].items()}
adj=defaultdict(set)
for q,a in cells.items():
    vs=[U.zadd(q,U.zt(a+2*i)) for i in range(5)]
    for i in range(5):
        u,w=vs[i],vs[(i+1)%5]; adj[u].add(w); adj[w].add(u)
XY={p:U.xy(p) for p in adj}
def opts(A,xy):
    O={}
    for cur in A:
        for prev in A[cur]:
            outs=[w for w in A[cur] if w!=prev]
            if not outs: O[(prev,cur)]={}; continue
            if xy is None: ws=sorted(outs)
            else:
                back=math.atan2(xy[prev][1]-xy[cur][1],xy[prev][0]-xy[cur][0])
                ws=[w for _,w in sorted(((math.atan2(xy[w][1]-xy[cur][1],xy[w][0]-xy[cur][0])-back)%(2*math.pi),w) for w in outs)]
            if len(ws)==1: O[(prev,cur)]={0:ws[0]}
            elif len(ws)==2: O[(prev,cur)]={1:ws[0],-1:ws[1]}
            else: O[(prev,cur)]={1:ws[0],0:ws[1],-1:ws[2]}
    return O
O=opts(adj,XY)
stars=[cyc for area,cyc in U.gaps(cells) if abs(area-2.9389)<0.01]
STARS=[]
for cyc in stars:
    c10=(0,0,0,0)
    for p in cyc: c10=U.zadd(c10,p)
    cx,cy=U.xy(c10)[0]/10,U.xy(c10)[1]/10
    P=[U.xy(p) for p in cyc]
    ar=sum(P[i][0]*P[(i+1)%10][1]-P[(i+1)%10][0]*P[i][1] for i in range(10))
    if ar<0: cyc=cyc[::-1]
    # 星の向き：中心から角へ向かう方位が 18+36k のどちらの組か
    a0=(math.degrees(math.atan2(XY[cyc[0]][1]-cy,XY[cyc[0]][0]-cx))-18)%72
    STARS.append({"c10":c10,"R":round(math.hypot(cx,cy),4),
                  "or":0 if min(a0,72-a0)<1e-6 else 1,
                  "starts":[(p,cyc[(i+1)%10]) for i,p in enumerate(cyc)]})
STARS.sort(key=lambda s:(s["R"],math.atan2(U.xy(s["c10"])[1],U.xy(s["c10"])[0])))
ST=[(si,ci,st) for si,s in enumerate(STARS) for ci,st in enumerate(s["starts"])]
print("五芒星の向き（画面方位の2種）ごとの個数 %s" % dict(Counter(s["or"] for s in STARS)))

def accepts(Om, starts, w):
    W=[]
    for i,(si,ci,st) in enumerate(starts):
        prev,cur=st
        ok=True
        for x in w:
            o=Om[(prev,cur)]
            if x not in o: ok=False; break
            prev,cur=cur,o[x]
        if ok: W.append(i)
    return frozenset(W)

rng=random.Random(20260903)
words=[tuple(rng.choice((-1,0,1)) for _ in range(N)) for _ in range(M)]
print("\n■ 必ず落ちる設定：全部0の語 → |W| = %d" % len(accepts(O,ST,(0,)*N)))

Ws=[accepts(O,ST,w) for w in words]
nz=[W for W in Ws if W]
print("\n■ 本番（長さ%d・無作為な語 %d 本）" % (N,M))
print("  受理された語 %d 本（%.2f%%）" % (len(nz),100*len(nz)/M))
sz=sorted(len(W) for W in nz)
print("  |W| … 最小 %d／中央 %d／最大 %d／平均 %.1f" % (sz[0],sz[len(sz)//2],sz[-1],sum(sz)/len(sz)))
print("  |W| が 5 の倍数でないもの %d 本（C₅対称なら 0 のはず）" % sum(1 for s in sz if s%5))
print("  相異なる W %d 種 / 受理語 %d 本" % (len({W for W in nz}), len(nz)))
# W は殻の中で閉じるか／向きの中で閉じるか
sh=[STARS[si]["R"] for si,ci,st in ST]; orr=[STARS[si]["or"] for si,ci,st in ST]
c1=sum(1 for W in nz if len({sh[i] for i in W})==1)
c2=sum(1 for W in nz if len({orr[i] for i in W})==1)
print("  W が一つの殻に収まる語 %d / %d（%.1f%%）" % (c1,len(nz),100*c1/len(nz)))
print("  W が一つの向きに収まる語 %d / %d（%.1f%%）" % (c2,len(nz),100*c2/len(nz)))
print("  殻ごとの受理のされやすさ（その殻の出発が W に入る割合）")
for R in sorted(set(sh)):
    idx=[i for i in range(210) if sh[i]==R]
    tot=sum(len([i for i in W if sh[i]==R]) for W in nz)
    print("    殻 %-8s 出発%3d本  平均 %.2f 本が受理（%.1f%%）" % (R,len(idx),tot/len(nz),100*tot/len(nz)/len(idx)))

# 負の対照
def config_graph(seed):
    r=random.Random(seed); stubs=[]
    for p in adj: stubs+=[p]*len(adj[p])
    for _ in range(300):
        r.shuffle(stubs); B={p:set() for p in adj}; ok=True
        for i in range(0,len(stubs),2):
            a,b=stubs[i],stubs[i+1]
            if a==b or b in B[a]: ok=False; break
            B[a].add(b); B[b].add(a)
        if ok and all(len(B[p])==len(adj[p]) for p in adj): return B
G=config_graph(11); OG=opts(G,None)
gs=[]
r2=random.Random(5); vs=list(G)
while len(gs)<210:
    p=r2.choice(vs); q=r2.choice(sorted(G[p]))
    gs.append((0,0,(p,q)))
Wg=[accepts(OG,gs,w) for w in words]; nzg=[W for W in Wg if W]
szg=sorted(len(W) for W in nzg)
print("\n■ 負の対照（同じ次数列の配置換えグラフ・幾何なし）")
print("  受理された語 %d 本（%.2f%%）" % (len(nzg),100*len(nzg)/M))
print("  |W| … 最小 %d／中央 %d／最大 %d／平均 %.1f" % (szg[0],szg[len(szg)//2],szg[-1],sum(szg)/len(szg)))
print("  相異なる W %d 種 / 受理語 %d 本" % (len({W for W in nzg}), len(nzg)))
