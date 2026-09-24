# 道具2（その三）：下端の規則（桁番号の引き算 i−a）を外し、中心 C の側へ内向きに描いた型紙（φ^−3k、小数部）で受ける
#  事前の基準
#   T1 内向きの型紙 k=1..40 で一桁の表が基準と一致（40/40）
#   T2 整数部だけ（小数部の入力0）の 10/20/40 桁：貫く入力と無作為100組で和が一致、段数は その二 と同じ
#   T3 整数部5桁＋小数部10桁の足し算 100組で和が一致（5^10 倍した整数で照合）
#   規則として残すもの：相手の目印が描かれていなければ何も渡らない（恒等）。桁番号の引き算は使わない
#   対照：比を φ だけずらすと和が落ちる
import random
exec(open('st2.py').read().split('def fastadd')[0])
PHIm1=U.zsub(U.PHI,U.ONE)
_up=upow
def upow(n):
    if n>=0: return _up(n)
    p=U.ONE
    for _ in range(-n): p=U.zmul(p,PHI3i)
    return p

def ref2(i): return U.zadd(C2,U.zmul(upow(i),D0))
# 内向きは座標の係数が大きく、浮動小数の升目分けが効かない。隣の候補は桁0の型紙の番号の対から取り、
# 描いた座標で norm2 を整数で判定し直す（全部そろって一致すること、を数える）
_src=open('vp2.py').read().split('def draw2(i):')[1].split('\nexec(')[0]
_base=None
def draw_in(i):
    global _base
    u=upow(i); P2=lambda v:U.zadd(C2,U.zmul(u,U.zsub(U.zadd(v,v),C2)))
    Z0=tuple(z0); tpl=TPL
    cells=[P2(U.zadd(tuple(v),Z0)) for v in tpl['cells']]
    need=U.norm2(U.zmul(u,(2*U.PHI[0],2*U.PHI[1],2*U.PHI[2],2*U.PHI[3])))
    if _base is None:
        G0=draw2(0); c0=[P2b(v) for v in tpl['cells']] if False else None
    adj=defaultdict(list); bad=0
    for a,b in BASEPAIRS:
        if U.norm2(U.zsub(cells[a],cells[b]))==need: adj[cells[a]].append(cells[b])
        else: bad+=1
    assert bad==0, bad
    S=set(cells)
    sa=P2(U.zadd(tuple(tpl['sa']),Z0)); sb=P2(U.zadd(tuple(tpl['sb']),Z0)); scr=[P2(U.zadd(tuple(v),Z0)) for v in tpl['screen']]
    def ring(g):
        for kk in (0,1):
            r=[U.zadd(g,U.zmul(u,U.zmul((2,0,0,0),U.zmul(U.PHI,U.zt(kk+2*j))))) for j in range(5)]
            if all(x in S for x in r): return r
    rA,rB=ring(sa),ring(sb)
    def auto(src,init,op):
        st={q:init for q in src}; fr=list(src)
        while fr:
            nx=[]
            for a in fr:
                for b in adj[a]:
                    if b not in st: st[b]=op(st[a]); nx.append(b)
            fr=nx
        return st
    base=U.zsub(rA[0],sa)
    return dict(scr=scr,PA=auto(rA,base,ZR),PB=auto(rB,base,ZR),SA=auto(rA,(1,0),TW),SB=auto(rB,(1,0),TW),
                dA=auto(rA,0,lambda n:n+1),dB=auto(rB,0,lambda n:n+1))
_u0=U.ONE; _c0=[U.zadd(C2,U.zsub(U.zadd(U.zadd(tuple(v),tuple(z0)),U.zadd(tuple(v),tuple(z0))),C2)) for v in TPL['cells']]
_need0=U.norm2((2*U.PHI[0],2*U.PHI[1],2*U.PHI[2],2*U.PHI[3]))
BASEPAIRS=[(a,b) for a in range(len(_c0)) for b in range(len(_c0)) if a!=b and U.norm2(U.zsub(_c0[a],_c0[b]))==_need0]
print("桁0の隣の対（全対を整数で判定）",len(BASEPAIRS))
print("確認：桁0と桁3を draw_in で描いた表が基準と一致", table(draw_in(0))==tabs[0], table(draw_in(3))==tabs[0])
tin=[table(draw_in(-k)) for k in range(1,41)]

print(f"T1 内向きの型紙 k=1..40 で表が基準と一致 {sum(t==tabs[0] for t in tin)}/40")
LO=40
TAB={i:(tabs[i] if i>=0 else tin[-i-1]) for i in range(-LO,41)}
IDX2={ref2(i):i for i in range(-LO,41)}
def fadd(xd,yd,lo,hi,wrong=False):     # 桁 lo..hi−1 に入力、その下 −LO..lo−1 は 0
    rng=list(range(-LO,hi)); X=dict(zip(range(lo,hi),xd)); Y=dict(zip(range(lo,hi),yd))
    f={i:tuple(TAB[i][(X.get(i,0),Y.get(i,0),c)][1] for c in (-1,0,1)) for i in rng}
    A=dict(f); P=dict(f); a,b=1,1; ra=rb=PHI3i; steps=0
    while a<hi+LO:
        s=U.zmul(ra,PHIm1) if wrong else ra; nA={}
        for i in rng:
            j=IDX2.get(U.zadd(C2,U.zmul(s,U.zsub(ref2(i),C2))))
            q=P[j] if j is not None else (-1,0,1)
            nA[i]=tuple(A[i][q[c]+1] for c in range(3))
        P=A; A=nA; a,b=a+b,a; ra,rb=U.zmul(ra,rb),ra; steps+=1
    c={i+1:A[i][1] for i in rng}; c[-LO]=0
    s={i:TAB[i][(X.get(i,0),Y.get(i,0),c[i])][0] for i in range(lo,hi)}
    return steps+1, s, c[hi]
def val(s,lo,hi,cN): return sum(s[i]*5**(i-lo) for i in range(lo,hi))+cN*5**(hi-lo)
def num(d): return sum(v*5**i for i,v in enumerate(d))
random.seed(6)
for nd in (10,20,40):
    xd=[2]+[1]*(nd-1); k,s,cN=fadd(xd,xd,0,nd); kw,sw,cw=fadd(xd,xd,0,nd,wrong=True)
    lim=(5**nd-1)//2; ok=0
    for _ in range(100):
        X,Y=random.randint(-lim,lim),random.randint(-lim,lim); kk,ss,cc=fadd(todig_n(X,nd),todig_n(Y,nd),0,nd); ok+=val(ss,0,nd,cc)==X+Y
    print(f"T2 {nd}桁 貫く入力 和 {val(s,0,nd,cN)==2*num(xd)} 段数 {k}   無作為 {ok}/100   対照（比ずらし）和 {val(sw,0,nd,cw)==2*num(xd)}")
lim=(5**15-1)//2; ok=0; ks=set()
for _ in range(100):
    X,Y=random.randint(-lim,lim),random.randint(-lim,lim)          # 5^10 倍した整数：下10桁が小数部
    kk,ss,cc=fadd(todig_n(X,15),todig_n(Y,15),-10,5); ok+=val(ss,-10,5,cc)==X+Y; ks.add(kk)
print(f"T3 整数部5桁＋小数部10桁 和 {ok}/100 段数 {sorted(ks)}")
