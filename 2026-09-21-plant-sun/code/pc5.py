# 設計 PC5：10枚の担体で平衡5進5桁の加算器
#  主の扇 w=0,2,4,6,8 が桁 0..4。各桁は自分の扇の五角形だけで動く（翼0の間隔4の鏡の対を ζ^w で回した対）
#  桁の計算（PS1 の構造オートマトン）：一歩で指しを ζ¹、表裏を入れ替え（裏→表で符号反転）
#    入力 x,y と繰り上がり c_in：A を -2(x+c_in) 歩、B を +2y 歩 遅らせる
#    s ＝ 指しが一致する着地の段 d/2、繰り上がりの大きさ＝表裏の不一致
#    繰り上がりの符号：s≠0 なら -sign(s)。s=0 で不一致なら c_in の符号を通す（通過の規則）
#  繰り上がりの線：主 w → φ² → ずらし w+1 の中を歩く → φ² → 主 w+2。符号は振幅 ±1（半周）で運ぶ
#  最上位（主8）の繰り上がりはずらし9を渡った先で読み、5⁵ の桁とする（輪はここで切る）
# 事前の基準
#   D1 5つの扇すべてで、一桁の表（x,y,c_in の 5×5×3=75通り）が同じ s,c になり、x+y+c_in=5c+s
#   D2 繰り上がりの線：5本とも主 w から主 w+2 へ届き、歩数が5本で等しく、符号 ±1 がそのまま届く
#   D3 5桁：端の組（全部+2 同士・全部−2 同士・交互）と無作為 20000 組で X+Y = 3125·c_top + Σ 5^i s_i
#   対照 ねじれなし／通過の規則なし（s=0 なら c=0）→ D1・D3 が崩れる
#   NG 一つでも食い違う
import sys,pickle,random; sys.path.insert(0,'.')
import b13_chain_units as U
from collections import defaultdict,deque
d=pickle.load(open('/home/claude/prop/g10.pkl','rb')); Q=d['Q']; cw=d['cell_w']; ADJ=d['adj']; z0=d['z0']; S=set(Q)
s14=pickle.load(open('/home/claude/prop/slit14.pkl','rb'))
def rotw(p,k): return U.zadd(U.zrot(U.zsub(p,z0),k),z0)
wing={w:{q for q in Q if w in cw[q]} for w in range(10)}
def ring5(g): 
    for kk in (0,1):
        r=[U.zadd(g,U.zmul(U.PHI,U.zt(kk+2*i))) for i in range(5)]
        if all(x in S for x in r): return r
sa0,sb0=s14['sa'],s14['sb']; scr0=s14['screen']
assert all(q in wing[0] for q in ring5(sa0)+ring5(sb0)+scr0)
TW=lambda s:(-s[1],s[0]); NT=lambda s:(s[1],s[0]); ZR=lambda z:U.zmul(z,U.zt(1))
def auto(src,cells,init,op):
    st={q:init for q in src}; fr=list(src)
    while fr:
        nx=[]
        for u in fr:
            for v in ADJ.get(u,()):
                if v in cells and v not in st: st[v]=op(st[u]); nx.append(v)
        fr=nx
    return st
def rep(f,z,n):
    for _ in range(n): z=f(z)
    return z
# 各桁の扇の準備
D={}
for i,w in enumerate((0,2,4,6,8)):
    sa=rotw(sa0,w); sb=rotw(sb0,w); scr=[rotw(q,w) for q in scr0]
    base=U.zsub(ring5(sa)[0],sa)
    rA,rB=ring5(sa),ring5(sb)
    D[i]=dict(scr=scr,
      PA=auto(rA,wing[w],base,ZR), PB=auto(rB,wing[w],base,ZR),
      SA=auto(rA,wing[w],(1,0),TW), SB=auto(rB,wing[w],(1,0),TW),
      SAn=auto(rA,wing[w],(1,0),NT), SBn=auto(rB,wing[w],(1,0),NT),
      dA=auto(rA,wing[w],0,lambda n:n+1), dB=auto(rB,wing[w],0,lambda n:n+1))
def digit(i,x,y,cin,twist=True,passrule=True):
    g=D[i]; a=-2*(x+cin); b=2*y
    SA,SB,op=(g['SA'],g['SB'],TW) if twist else (g['SAn'],g['SBn'],NT)
    L=[q for q in g['scr'] if rep(ZR,g['PA'][q],a%10)==rep(ZR,g['PB'][q],b%10)]
    lv={(g['dA'][q]-g['dB'][q])//2 for q in L}
    if len(lv)!=1: return None
    s=lv.pop(); same={rep(op,SA[q],a%4)==rep(op,SB[q],b%4) for q in L}
    if len(same)!=1: return None
    if same.pop(): return s,0
    if s!=0: return s,(-1 if s>0 else 1)
    return (s,(1 if cin>0 else -1)) if (passrule and cin!=0) else (s,0)
# D1
R=range(-2,3)
tabs=[]
for i in range(5):
    t={(x,y,c):digit(i,x,y,c) for x in R for y in R for c in (-1,0,1)}
    ok=sum(1 for (x,y,c),v in t.items() if v and x+y+c==5*v[1]+v[0])
    tabs.append(t); print(f'D1 桁{i}（主の扇{2*i}） {ok}/75')
print('D1 5つの扇の表が同じ',all(tabs[i]==tabs[0] for i in range(5)))
# D2 繰り上がりの線：φ² の繋がり
by=defaultdict(list)
for q in Q:
    p=U.xy(q); by[(int(p[0]//3),int(p[1]//3))].append(q)
L2=defaultdict(set)
for q in Q:
    p=U.xy(q); gx,gy=int(p[0]//3),int(p[1]//3)
    for dx in (-1,0,1):
        for dy in (-1,0,1):
            for r in by.get((gx+dx,gy+dy),()):
                if (min(cw[q])%2)!=(min(cw[r])%2) and U.norm2(U.zsub(q,r))==(2,3): L2[q].add(r)
lines=[]
for i,w in enumerate((0,2,4,6,8)):
    m=(w+1)%10; t=(w+2)%10
    starts=[(q,r) for q in wing[w] for r in L2[q] if r in wing[m]]
    # 桁の扇の端（φ² の口）から、ずらしの扇の中を歩き、φ² で次の主の扇に出るまでの最短の歩数
    dist={}; dq=deque()
    for q,r in starts: 
        if r not in dist: dist[r]=1; dq.append(r)
    while dq:
        u=dq.popleft()
        for v in ADJ.get(u,()):
            if v in wing[m] and v not in dist: dist[v]=dist[u]+1; dq.append(v)
    arr=[dist[u]+1 for u in dist for r in L2[u] if r in wing[t]]
    other=[r for u in dist for r in L2[u] if not (r in wing[t] or r in wing[w] or r in wing[m])]
    amp=[rep(lambda s:-s,sg,0) for sg in (1,-1)]
    lines.append((min(arr),len(arr),len(other)))
    print(f'D2 線{i}  主{w}→ずらし{m}→主{t}  最短の歩数 {min(arr)}  出口の数 {len(arr)}  他の扇へ漏れる口 {len(other)}  符号 +1→{amp[0]} −1→{amp[1]}')
print('D2 5本の歩数が等しい',len({l[0] for l in lines})==1)
# D3
def add5(X,Y,twist=True,passrule=True):
    c=0; out=0
    for i in range(5):
        v=digit(i,X[i],Y[i],c,twist,passrule)
        if v is None: return None
        s,c=v; out+=s*5**i
    return out+3125*c
def val(X): return sum(X[i]*5**i for i in range(5))
cases=[[2]*5,[-2]*5,[2,-2]*2+[2],[-2,2]*2+[-2]]
T=[(a,b) for a in cases for b in cases]
random.seed(0); T+=[([random.randint(-2,2) for _ in range(5)],[random.randint(-2,2) for _ in range(5)]) for _ in range(20000)]
for tw,pr,lab in ((True,True,'設計どおり'),(False,True,'対照 ねじれなし'),(True,False,'対照 通過の規則なし')):
    ok=sum(add5(a,b,tw,pr)==val(a)+val(b) for a,b in T)
    print(f'D3 {lab}: {ok}/{len(T)}')
