#  検定PL14  東から西へ渡る太陽で枝分かれがくり返されるか（14行）
#   太陽は十分遠いとして平行光：まわり5枚の明るさ＝太陽方向への射影（丸めて比べる）
#   1日＝13コマ：方位 0°(東の地平)→18°刻みで上を通り→180°(西の地平) の11コマ＋夜2コマ
#   （画面座標 y下向き。上＝270°）
#   事前の予想（走らせる前に書き直した）
#     平行光では、二つの最大が同点になるのは太陽が星の尖りの方向にあるとき
#     a の尖り 18,90,162,234,306 → 空にあるのは 306(朝)・234(午後)
#     b の尖り 54,126,198,270,342 → 空にあるのは 342(日の出直後)・270(正午)・198(日没前)
#     → 1日の中で b,a,b,a,b の順に分岐の機会が来る
#   OK   二本に分かれる分岐が、a では 306・234、b では 342・270・198 のときだけ起きる
#   NG   それ以外の方位で分かれる
#   対照 もやし：太陽を真上(270)に止めたまま同じ日数
#   空試験 夜だけ：伸びない
import math, itertools, os
from collections import Counter, defaultdict
src=open('pl11.py').read(); exec(src[:src.index('A=Counter()')])
DAY=[0,342,324,306,288,270,252,234,216,198,180,None,None]
def tips_par(i,deg):
    ux,uy=math.cos(math.radians(deg)),math.sin(math.radians(deg))
    v=[round(p[0]*ux+p[1]*uy,6) for q,p in AROUND[i]]
    m=max(v); return [AROUND[i][t][1] for t,x in enumerate(v) if x==m]
def step_to(i,target):
    cx,cy=SC[i]; th=math.atan2(target[1]-cy,target[0]-cx)
    best=bd=None
    for j in LK[i]:
        d=abs((math.atan2(SC[j][1]-cy,SC[j][0]-cx)-th+math.pi)%(2*math.pi)-math.pi)
        if bd is None or d<bd: bd,best=d,j
    return best if best is not None and bd<math.pi/2 else None
def run(start, sched, days):
    live={start}; ev=[]; perday=[]
    for d in range(days):
        for f,deg in enumerate(sched):
            if deg is None: continue
            nxt=set()
            for i in live:
                tg=tips_par(i,deg)
                if len(tg)>=3: nxt.add(i); continue      # 伸びない（その場に留まる）
                dest={step_to(i,t) for t in tg}; dest.discard(None)
                if len(tg)==2 and len(dest)==2: ev.append((d,f,deg,ori[i]))
                if dest: nxt|=dest
                else: nxt.add(i)                          # 行き先なし＝その場に留まる
            live=nxt
        perday.append(len(live))
    return ev, perday
ymax=max(y for _,y in SC)
ground=[i for i in range(len(SC)) if abs(SC[i][1]-ymax)<1e-6]
for name,sched in (("東→西",DAY),("もやし(真上固定)",[270]*11+[None,None]),("夜だけ",[None]*13)):
    E=Counter(); tips=[]
    for s in ground:
        ev,pd=run(s,sched,3); tips.append(pd)
        for d,f,deg,o_ in ev: E[(o_,deg)]+=1
    print(f"{name}  地面の星{len(ground)}個・3日  分岐（向き,方位）: {dict(sorted(E.items()))}")
    print(f"   日ごとの先端の数（出発ごと）: {tips}")
    if name=="東→西":
        okset={('a',306),('a',234),('b',342),('b',270),('b',198)}
        print("   判定", "OK" if E and set(E)<=okset else "NG")
