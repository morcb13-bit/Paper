# 14行の三角形担体で、二つのスリット（五芒星の対）から出た波の着地を整数だけで数える
import os, sys
os.environ['B13_ROWS']='14'
sys.path.insert(0,'/home/claude/Paper/2026-09-21-plant-sun/code')
import b13_chain_units as U
from collections import deque, defaultdict
rows,place,offs=U.build_stack()
cells=U.fits(sum(place,[]))
Q=list(cells)
# 隣（辺を共有する五角形）：中心差のノルムが NCELL と整数で一致
by=defaultdict(list)
for q in Q:
    p=U.xy(q); by[(int(p[0]//2),int(p[1]//2))].append(q)
adj=defaultdict(list)
for q in Q:
    p=U.xy(q); gx,gy=int(p[0]//2),int(p[1]//2)
    for dx in (-1,0,1):
        for dy in (-1,0,1):
            for r in by.get((gx+dx,gy+dy),()):
                if r!=q and U.norm2(U.zsub(q,r))==U.NCELL: adj[q].append(r)
# 五芒星：5枚の五角形に距離φで囲まれる点（整数で判定）
S=set(Q); stars={}
for q in Q:
    for m in range(10):
        g=U.zsub(q,U.zmul(U.PHI,U.zt(m)))
        for k in (0,1):
            ring=[U.zadd(g,U.zmul(U.PHI,U.zt(k+2*i))) for i in range(5)]
            if all(r in S for r in ring): stars[g]=ring
print('五角形',len(Q),'五芒星',len(stars))
# 鏡：三角形の軸（頂点の円環の中心を通る縦線）で折る。対応づけにだけ描画座標を使う
x0=U.xy(place[0][0])[0]
key={ (round(U.xy(q)[0],4), round(U.xy(q)[1],4)) : q for q in Q}
def mir(q):
    x,y=U.xy(q); return key.get((round(2*x0-x,4), round(y,4)))
mirror_ok=sum(1 for q in Q if mir(q) is not None)
print('鏡の相手がある五角形',mirror_ok,'/',len(Q))
def bfs(src):
    d={s:0 for s in src}; dq=deque(src)
    while dq:
        u=dq.popleft()
        for w in adj[u]:
            if w not in d: d[w]=d[u]+1; dq.append(w)
    return d
bottom_y=max(U.xy(q)[1] for q in Q)
screen=[q for q in Q if U.xy(q)[1]>bottom_y-1.0]
BASE=3120; STEP=780                                # 一歩で位相 780（4歩で一周）
def status(dA,dB,delta):
    ph=(STEP*(dA-dB)+delta)%BASE
    return 'B' if ph==0 else ('D' if ph==BASE//2 else '-')
def trial(sa,sb,label):
    dA=bfs(stars[sa]); dB=bfs(stars[sb])
    T1=all(dA[q]==dB[mir(q)] for q in screen if mir(q) is not None)
    st0={q:status(dA[q],dB[q],0) for q in screen}
    st1={q:status(dA[q],dB[q],1560) for q in screen}
    T2=all(st0[q]==st0[mir(q)] for q in screen if mir(q) is not None)
    eq=[q for q in screen if dA[q]==dB[q]]
    T3=all(st0[q]=='B' and st1[q]=='D' for q in eq)
    flip=sum(1 for q in screen if (st0[q],st1[q]) in (('B','D'),('D','B')))
    print(f'{label}: T1 鏡で距離が入れ替わる={T1}  T2 縞が左右対称={T2}  T3 等距離{len(eq)}個が明→暗={T3}  '
          f'明{sum(v=="B" for v in st0.values())} 暗{sum(v=="D" for v in st0.values())} → ずらすと 明{sum(v=="B" for v in st1.values())} 暗{sum(v=="D" for v in st1.values())}  反転{flip}/{len(screen)}')
    return st0,st1,dA,dB
# 鏡の対になる五芒星を選ぶ（中ほどの高さ）
cands=[]
for g in stars:
    x,y=U.xy(g)
    if x<x0-1e-6:
        gm=[h for h in stars if abs(U.xy(h)[0]-(2*x0-x))<1e-4 and abs(U.xy(h)[1]-y)<1e-4]
        if gm: cands.append((y,g,gm[0]))
cands.sort()
print('鏡の対になる五芒星の組',len(cands))
y,sa,sb=cands[len(cands)//3]
st0,st1,dA,dB=trial(sa,sb,'鏡の対')
# 負の対照：鏡の対でない二つの五芒星
other=[h for h in stars if h not in (sa,sb) and abs(U.xy(h)[1]-y)>0.5]
trial(sa,other[0],'対照（鏡の対でない）')
import pickle; pickle.dump(dict(Q=Q,cells=cells,adj=dict(adj),stars=stars,sa=sa,sb=sb,screen=screen,st0=st0,st1=st1),open('/home/claude/prop/slit14.pkl','wb'))
