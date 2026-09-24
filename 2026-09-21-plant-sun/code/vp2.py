# 検定VP2 次の桁の型紙の置き場を、ひし形の入れ子（同じ中心・比 φ³）が決める
#  NR2/NR3：細ひし形は、同じ中心をもつ円環の中心の黄金のひし形と比 φ³（回転なし）で相似（帯のひし形は次の行の円環が要る）
#  置き方：桁 i+1 の型紙 = 桁 i の型紙を、軸上の帯のひし形の中心 C のまわりに φ³ 倍（2x' = M + φ³(2x−M)、整数で持つ）
#   平行移動・回転の量は外から与えない。中心 C と比 φ³ はひし形の組から読む
# 事前の基準
#   T3 桁 0..40 の型紙すべてで一桁の表が基準と一致
#   T4 この置き方で描き足す加算器の 10/20/40 桁の足し算（各100組）が整数の和と一致
#   確認 桁 i の軸上の帯のひし形の像が、桁 i の担体の外に次の行を足したときの大きなひし形と一致する（NR3 の相手そのもの）
import sys,json,math,random
exec(open('/home/claude/vp1.py').read().split('# ---- 描く')[0])
d=pickle.load(open('/home/claude/prop/slit14.pkl','rb'))
F=U.gaps(d['cells']); rh=[cyc for a,cyc in F if min(U.GAP_NAME,key=lambda x:abs(x-a))==0.8123]
x0=U.xy(z0)[0]
ax=[c for c in rh if abs(sum(U.xy(v)[0] for v in c)/4-x0)<1e-6 and max(U.xy(v)[1] for v in c)>78]
C2=U.zadd(ax[0][0],ax[0][2])                      # 2C
PHI3=U.zmul(U.zmul(U.PHI,U.PHI),U.PHI)
ZR=lambda z:U.zmul(z,U.zt(1)); TW=lambda s:(-s[1],s[0])
print('軸上の帯のひし形',len(ax),' 中心の2倍 M',C2)
def upow(n):
    p=U.ONE
    for _ in range(n): p=U.zmul(p,PHI3)
    return p
def draw2(i):
    u=upow(i); P2=lambda v:U.zadd(C2,U.zmul(u,U.zsub(U.zadd(v,v),C2)))   # 2倍の座標
    Z0=tuple(z0); tpl=TPL
    cells=[P2(U.zadd(tuple(v),Z0)) for v in tpl['cells']]
    N2=4*U.norm2(u)[0], 4*U.norm2(u)[1]
    need=U.norm2(U.zmul(u,(2*U.PHI[0],2*U.PHI[1],2*U.PHI[2],2*U.PHI[3])))
    sc=(1+5**.5)/2**1; s=((1+5**.5)/2)**(3*i)*2
    by=defaultdict(list)
    for q in cells:
        x,y=U.xy(q); by[(math.floor(x/(2*s)),math.floor(y/(2*s)))].append(q)
    adj=defaultdict(list)
    for q in cells:
        x,y=U.xy(q); gx,gy=math.floor(x/(2*s)),math.floor(y/(2*s))
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for r in by.get((gx+dx,gy+dy),()):
                    if r!=q and U.norm2(U.zsub(q,r))==need: adj[q].append(r)
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
exec(open('/home/claude/vp1.py').read().split('# ---- 描く（ここから担体を使わない）----')[1].split('ref=table')[0].split("def rep")[1].join(["def rep",""]) if False else "")
def rep(f,z,n):
    for _ in range(n): z=f(z)
    return z
def digit(G,x,y,cin):
    a=-2*(x+cin); b=2*y
    L=[q for q in G['scr'] if rep(ZR,G['PA'][q],a%10)==rep(ZR,G['PB'][q],b%10)]
    lv={(G['dA'][q]-G['dB'][q])//2 for q in L}
    if len(lv)!=1: return None
    s=lv.pop(); same={rep(TW,G['SA'][q],a%4)==rep(TW,G['SB'][q],b%4) for q in L}
    if len(same)!=1: return None
    if same.pop(): return s,0
    if s!=0: return s,(-1 if s>0 else 1)
    return s,((1 if cin>0 else -1) if cin else 0)
R=range(-2,3)
table=lambda G:{(x,y,c):digit(G,x,y,c) for x in R for y in R for c in (-1,0,1)}
ref=table(draw2(0))
print('桁0の表が和と整合',all(v and x+y+c==5*v[1]+v[0] for (x,y,c),v in ref.items()))
tabs=[ref]
for i in range(1,41): tabs.append(table(draw2(i)))
print(f'T3 入れ子に置いた桁 1..40 の型紙で表が基準と一致 {sum(t==ref for t in tabs[1:])}/40')
def todig(X):
    out=[]
    while X: r=((X+2)%5)-2; out.append(r); X=(X-r)//5
    return out or [0]
def add_n(X,Y):
    xd,yd=todig(X),todig(Y); n=max(len(xd),len(yd)); c=0; out=[]; i=0
    while i<n or c:
        x=xd[i] if i<len(xd) else 0; y=yd[i] if i<len(yd) else 0
        s,c=tabs[i][(x,y,c)]; out.append(s); i+=1
    return sum(s*5**i for i,s in enumerate(out))
random.seed(6)
for nd in (10,20,40):
    lim=(5**nd-1)//2
    print(f'T4 {nd}桁 ',sum(add_n(X,Y)==X+Y for X,Y in [(random.randint(-lim,lim),random.randint(-lim,lim)) for _ in range(100)]),'/100')
