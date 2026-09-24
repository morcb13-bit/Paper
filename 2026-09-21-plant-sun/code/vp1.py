# 検定VP1 型紙だけで描く一桁の加算器と、繰り上がりで桁を描き足す加算器
#  型紙：翼0の五角形721枚・スリットの五芒星二つ・着地35枚を、頂点（z0）からの Z[ζ] 整数差分として保存（JSON）
#  描く：任意の平行移動 T（整数4成分）と回転 ζ^k で p' = T + ζ^k·(p−z0)。隣接は描いた後に norm2 = NCELL で判定し直す
#  担体（build_stack, gen10）は型紙を作るときに一度読むだけ。描いてからは使わない
# 事前の基準
#   T1 無作為な T と k（0..9）で描いた 12 個の型紙すべてで、一桁の表 75 通りが PC5 の表と一致
#   T2 繰り上がりが最上位を越えるたびに型紙を一枚描き足す加算器で、10桁・20桁・40桁の無作為な足し算（各100組）が整数の和と一致し、
#      描いた型紙の枚数が結果の桁数と一致
#   対照 スリットを鏡の対でない五芒星に替えた型紙 → 表が崩れる
import sys,json,pickle,random,math
sys.path.insert(0,'/home/claude/Paper/2026-09-21-plant-sun/code')
import b13_chain_units as U
from collections import defaultdict
# ---- 型紙を作る（担体を読むのはここだけ）----
d=pickle.load(open('/home/claude/prop/slit14.pkl','rb')); g10=pickle.load(open('/home/claude/prop/g10.pkl','rb'))
z0=g10['z0']; w0=[q for q in g10['Q'] if 0 in g10['cell_w'][q]]
rel=lambda p:list(U.zsub(p,z0))
def ring5_in(g,S):
    for kk in (0,1):
        r=[U.zadd(g,U.zmul(U.PHI,U.zt(kk+2*i))) for i in range(5)]
        if all(x in S for x in r): return r
S=set(w0); others=[g for g in d['stars'] if g not in (d['sa'],d['sb']) and abs(U.xy(g)[1]-U.xy(d['sa'])[1])>0.5]
TPL=dict(cells=[rel(q) for q in w0],sa=rel(d['sa']),sb=rel(d['sb']),screen=[rel(q) for q in d['screen']])
json.dump(TPL,open('/home/claude/template_digit.json','w'))
BAD=dict(TPL); BAD['sb']=rel(others[0])
# ---- 描く（ここから担体を使わない）----
def draw(tpl,T,k):
    P=lambda v:U.zadd(T,U.zmul(U.zt(k),tuple(v)))
    cells=[P(v) for v in tpl['cells']]; Sx=set(cells)
    by=defaultdict(list)
    for q in cells:
        x,y=U.xy(q); by[(int(x//2),int(y//2))].append(q)
    adj=defaultdict(list)
    for q in cells:
        x,y=U.xy(q); gx,gy=int(x//2),int(y//2)
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for r in by.get((gx+dx,gy+dy),()):
                    if r!=q and U.norm2(U.zsub(q,r))==U.NCELL: adj[q].append(r)
    sa,sb=P(tpl['sa']),P(tpl['sb']); scr=[P(v) for v in tpl['screen']]
    rA,rB=ring5_in(sa,Sx),ring5_in(sb,Sx)
    def auto(src,init,op):
        st={q:init for q in src}; fr=list(src)
        while fr:
            nx=[]
            for u in fr:
                for v in adj[u]:
                    if v not in st: st[v]=op(st[u]); nx.append(v)
            fr=nx
        return st
    ZR=lambda z:U.zmul(z,U.zt(1)); TW=lambda s:(-s[1],s[0])
    base=U.zsub(rA[0],sa)
    G=dict(scr=scr,PA=auto(rA,base,ZR),PB=auto(rB,base,ZR),SA=auto(rA,(1,0),TW),SB=auto(rB,(1,0),TW),
           dA=auto(rA,0,lambda n:n+1),dB=auto(rB,0,lambda n:n+1))
    return G
def rep(f,z,n):
    for _ in range(n): z=f(z)
    return z
ZR=lambda z:U.zmul(z,U.zt(1)); TW=lambda s:(-s[1],s[0])
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
def table(G): return {(x,y,c):digit(G,x,y,c) for x in R for y in R for c in (-1,0,1)}
ref=table(draw(TPL,(0,0,0,0),0))
okref=all(v and x+y+c==5*v[1]+v[0] for (x,y,c),v in ref.items())
print('基準の表（原点・回転0）75通りが和と整合',okref)
random.seed(5); same=0; placements=[]
for i in range(12):
    T=tuple(random.randint(-500,500) for _ in range(4)); k=random.randrange(10); placements.append((T,k))
    same+= table(draw(TPL,T,k))==ref
print(f'T1 無作為な平行移動と回転で描いた型紙 {same}/12 が同じ表')
print('対照 鏡の対でないスリットの型紙：表が同じ',table(draw(BAD,(0,0,0,0),0))==ref)
# ---- T2 桁を描き足す加算器 ----
drawn=[]
def get(i):
    while len(drawn)<=i:
        j=len(drawn); T=U.zmul(U.zt(j%10),(60*j,0,0,0)); drawn.append(table(draw(TPL,T,j%10)))
    return drawn[i]
def todig(X):
    out=[]
    while X: r=((X+2)%5)-2; out.append(r); X=(X-r)//5
    return out or [0]
def add_grow(X,Y):
    xd,yd=todig(X),todig(Y); n=max(len(xd),len(yd)); c=0; out=[]; i=0
    while i<n or c:
        x=xd[i] if i<len(xd) else 0; y=yd[i] if i<len(yd) else 0
        s,c=get(i)[(x,y,c)]; out.append(s); i+=1
    while len(out)>1 and out[-1]==0: out.pop()
    return sum(s*5**i for i,s in enumerate(out)), len(out)
for nd in (10,20,40):
    ok=okn=0
    for _ in range(100):
        lim=(5**nd-1)//2; X=random.randint(-lim,lim); Y=random.randint(-lim,lim)
        v,L=add_grow(X,Y); ok+= v==X+Y; okn+= L==len(todig(X+Y))
    print(f'T2 {nd}桁: 和が一致 {ok}/100  桁数と描いた型紙の数が一致 {okn}/100  描いた型紙の総数 {len(drawn)}')
