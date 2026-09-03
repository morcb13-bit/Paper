# 検定PI1（カーナビ方式）：離れた二つの五芒星から両側に広げて出会わせる
#   出すもの：最短の歩数 d／最短経路の本数／その桁和 S の分布／片側と両側の探索費用
#   OK なら：両側から広げると費用が √ になり、出会いの集合が経路の全体を与える
#   NG なら：出会いが起きない（非後退の規則で届かない）か、費用が下がらない
#   必ず落ちる設定：A=B なら d=0
import json, math, sys
from collections import defaultdict, Counter, deque
import b13_chain_units as U
phi=(1+5**0.5)/2
cells={tuple(int(x) for x in k.split(",")):a for k,a in json.load(open("carrier_1245.json"))["cells"].items()}
adj=defaultdict(set)
for q,a in cells.items():
    vs=[U.zadd(q,U.zt(a+2*i)) for i in range(5)]
    for i in range(5):
        u,w=vs[i],vs[(i+1)%5]; adj[u].add(w); adj[w].add(u)
XY={p:U.xy(p) for p in adj}
DIG={p:{2:-1,3:0,4:1}[len(adj[p])] for p in adj}
O={}
for cur in adj:
    for prev in adj[cur]:
        outs=[w for w in adj[cur] if w!=prev]
        if not outs: O[(prev,cur)]={}; continue
        back=math.atan2(XY[prev][1]-XY[cur][1],XY[prev][0]-XY[cur][0])
        ws=[w for _,w in sorted(((math.atan2(XY[w][1]-XY[cur][1],XY[w][0]-XY[cur][0])-back)%(2*math.pi),w) for w in outs)]
        O[(prev,cur)]={0:ws[0]} if len(ws)==1 else ({1:ws[0],-1:ws[1]} if len(ws)==2 else {1:ws[0],0:ws[1],-1:ws[2]})
CEN=[]; RIM=[]
for area,cyc in U.gaps(cells):
    if abs(area-2.9389)<0.01:
        c=(0,0,0,0)
        for p in cyc: c=U.zadd(c,p)
        if any(len(adj[p])<2 for p in cyc): continue
        CEN.append((U.xy(c)[0]/10,U.xy(c)[1]/10)); RIM.append(cyc)
print("五芒星 %d 個" % len(CEN))
def pick(target):
    return min(range(len(CEN)),key=lambda i:abs(math.hypot(*CEN[i])-target))
ia=pick(12.0); ib=None
best=None
for j in range(len(CEN)):
    d=math.hypot(CEN[j][0]-CEN[ia][0],CEN[j][1]-CEN[ia][1])
    if best is None or abs(d-40)<abs(best[0]-40): best=(d,j)
ib=best[1]
DA=math.hypot(CEN[ib][0]-CEN[ia][0],CEN[ib][1]-CEN[ia][1])
print("A=%s  B=%s  実距離 %.4f" % (tuple(round(x,3) for x in CEN[ia]),tuple(round(x,3) for x in CEN[ib]),DA))
SA=[(RIM[ia][i],RIM[ia][(i+1)%10]) for i in range(10)]
SB=set(RIM[ib])
# 片側の広がり（有向辺の集合）と、最短で B の縁に着く歩数
front={e:1 for e in SA}; ssum={e:0 for e in SA}
sizes=[]; step=0; hit=None
while step<60:
    step+=1
    nf=Counter(); ns=defaultdict(Counter)
    for e,c in front.items():
        for x,nxt in O[e].items():
            ne=(e[1],nxt); nf[ne]+=c
    front=dict(nf); sizes.append(len(front))
    if any(e[1] in SB for e in front):
        hit=step; break
print("最短の歩数 d = %s" % hit)
print("片側の前線の大きさ:", sizes)
# 最短経路の本数と桁和
front={e:{0:1} for e in SA}
for k in range(hit):
    nf=defaultdict(Counter)
    for e,ss in front.items():
        for x,nxt in O[e].items():
            ne=(e[1],nxt)
            for s,c in ss.items(): nf[ne][s+DIG[nxt]]+=c
    front={e:dict(v) for e,v in nf.items()}
tot=Counter()
for e,ss in front.items():
    if e[1] in SB:
        for s,c in ss.items(): tot[s]+=c
n=sum(tot.values())
print("最短経路 %d 本   桁和 S の分布 %s" % (n,dict(sorted(tot.items()))))
import cmath
w=cmath.exp(2j*math.pi/3)
amp=sum(c*w**s for s,c in tot.items())
print("符号つきの和 Σ ω^S = %.4f%+.4fi   |Σ| = %.4f   （打ち消しなしなら %d）" % (amp.real,amp.imag,abs(amp),n))
print("Σ(−1)^S = %d" % sum(c*(-1)**s for s,c in tot.items()))
# 両側から広げたときの費用
h=hit//2
fa={e:1 for e in SA}
for _ in range(h):
    nf=Counter()
    for e in fa:
        for x,nxt in O[e].items(): nf[(e[1],nxt)]+=1
    fa=dict(nf)
fb={(RIM[ib][(i+1)%10],RIM[ib][i]) for i in range(10)}
for _ in range(hit-h):
    nf=set()
    for e in fb:
        for x,nxt in O[e].items(): nf.add((e[1],nxt))
    fb=nf
print()
print("片側で %d 歩  → 前線 %d 個" % (hit,sizes[-1]))
print("両側で %d + %d 歩 → 前線 %d + %d = %d 個" % (h,hit-h,len(fa),len(fb),len(fa)+len(fb)))
