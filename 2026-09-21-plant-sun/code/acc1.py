# 仮想のアキュムレータ（担体なし）：レジスタを位相の組（指し ζ^k・表裏 TW^m、20で一周）で持ち、
# 書き戻しは着地した五角形での二つの位相場の比（PB/PA と SB/SA）を読んで次の周の源の回しにする。s を整数で渡さない
#  事前の基準
#   OK  10桁のレジスタに k=1..20 回足し込み、各回の中身が整数の累積和と一致（無作為100通り）／貫く入力でも一致
#   対照1 書き戻しの向きを逆（比の逆数）にすると落ちる   対照2 書き戻しをしないと2回目以降が落ちる
#  繰り上がりは前回の一括の渡し（φ^(−3a) の相似、二つの記憶）。繰り上がりの規則（通過の規則）はまだ値で持つ（§4 の残り）
#  照合（観察者側）だけでレジスタの位相を整数に直す
import random, pickle
src=open('vp2.py').read().split('ref=table(draw2(0))')[0]
exec(src)
PHIm1=U.zsub(U.PHI,U.ONE); PHI3i=U.zmul(U.zmul(PHIm1,PHIm1),PHIm1)
def upw(n):
    p=U.ONE
    for _ in range(abs(n)): p=U.zmul(p,PHI3 if n>0 else PHI3i)
    return p
_up=upow
def upow(n): return upw(n)
# 内向きの描画（隣の候補は桁0の全対判定から、判定は整数）
_c0=[U.zadd(U.zadd(tuple(v),tuple(z0)),U.zadd(tuple(v),tuple(z0))) for v in TPL['cells']]
_need0=U.norm2((2*U.PHI[0],2*U.PHI[1],2*U.PHI[2],2*U.PHI[3]))
BP=[(a,b) for a in range(len(_c0)) for b in range(len(_c0)) if a!=b and U.norm2(U.zsub(_c0[a],_c0[b]))==_need0]
def draw_any(i):
    if i>=0: return draw2(i)
    u=upw(i); P2=lambda v:U.zadd(C2,U.zmul(u,U.zsub(U.zadd(v,v),C2))); Z0=tuple(z0)
    cells=[P2(U.zadd(tuple(v),Z0)) for v in TPL['cells']]
    need=U.norm2(U.zmul(u,(2*U.PHI[0],2*U.PHI[1],2*U.PHI[2],2*U.PHI[3])))
    adj=defaultdict(list)
    for a,b in BP:
        assert U.norm2(U.zsub(cells[a],cells[b]))==need; adj[cells[a]].append(cells[b])
    S=set(cells); sa=P2(U.zadd(tuple(TPL['sa']),Z0)); sb=P2(U.zadd(tuple(TPL['sb']),Z0)); scr=[P2(U.zadd(tuple(v),Z0)) for v in TPL['screen']]
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
    return dict(scr=scr,PA=auto(rA,base,ZR),PB=auto(rB,base,ZR),SA=auto(rA,(1,0),TW),SB=auto(rB,(1,0),TW),dA=auto(rA,0,lambda n:n+1),dB=auto(rB,0,lambda n:n+1))
LO,HI=3,10
G={i:draw_any(i) for i in range(-LO,HI)}
Z0=tuple(z0); D0=U.zsub(U.zadd(U.zadd(tuple(TPL['sa']),Z0),U.zadd(tuple(TPL['sa']),Z0)),C2)
REF={i:U.zadd(C2,U.zmul(upw(i),D0)) for i in range(-LO,HI)}; IDX={v:i for i,v in REF.items()}
print("描いた型紙",len(G))
# 一桁：レジスタ (k,m) と入力 b で着地を読み、繰り上がりと、書き戻す位相の組を返す
_C={}
def dig(i,k,m,b,cin):
    key=(i,k,m,b,cin)
    if key in _C: return _C[key]
    g=G[i]; ka=(k-2*cin)%10; ma=(m-2*cin)%4
    L=[q for q in g['scr'] if rep(ZR,g['PA'][q],ka)==rep(ZR,g['PB'][q],b%10)]
    lv={(g['dA'][q]-g['dB'][q])//2 for q in L}; assert len(lv)==1
    s=lv.pop(); same={rep(TW,g['SA'][q],ma)==rep(TW,g['SB'][q],b%4) for q in L}; assert len(same)==1
    cout=0 if same.pop() else ((-1 if s>0 else 1) if s!=0 else ((1 if cin>0 else -1) if cin else 0))
    q=L[0]                                        # 書き戻し：着地した五角形での二つの位相場の比
    kk=[e for e in range(10) if rep(ZR,g['PA'][q],e)==g['PB'][q]][0]
    mm=[e for e in range(4) if rep(TW,g['SA'][q],e)==g['SB'][q]][0]
    _C[key]=(cout,(kk,mm)); return _C[key]
def add_step(reg,yb,mode):
    rng=range(-LO,HI); R=dict(reg); Y=dict(yb)
    for i in range(-LO,0): R[i]=(0,0); Y[i]=0
    f={i:tuple(dig(i,*R[i],2*Y[i],c)[0] for c in (-1,0,1)) for i in rng}
    A=dict(f); P=dict(f); ra=rb=PHI3i; st=0
    while not all(len(set(v))==1 for v in A.values()):
        nA={}
        for i in rng:
            j=IDX.get(U.zadd(C2,U.zmul(ra,U.zsub(REF[i],C2))))
            q=P[j] if j is not None else (-1,0,1)
            nA[i]=tuple(A[i][q[c]+1] for c in range(3))
        P=A; A=nA; ra,rb=U.zmul(ra,rb),ra; st+=1
    cin={i:(A[i-1][1] if i>0 else A[-1][1]) for i in range(HI)}
    out={}
    for i in range(HI):
        kk,mm=dig(i,*R[i],2*Y[i],cin[i])[1]
        out[i]={'wb':(kk,mm),'rev':((-kk)%10,(-mm)%4),'none':R[i]}[mode]
    return out, A[HI-1][1], st+1
# 観察者側の読み：位相の組 → 桁の値（a≡−2r を 10 と 4 で）
PH={((-2*r)%10,(-2*r)%4):r for r in range(-2,3)}
def val(reg): return sum(PH[reg[i]]*5**i for i in range(HI))
def todig(X,n):
    out=[]
    while X: r=((X+2)%5)-2; out.append(r); X=(X-r)//5
    return (out+[0]*n)[:n]
M=5**HI
def wrap(X): return ((X+(M-1)//2)%M)-(M-1)//2        # 10桁のレジスタは 5^10 で一周（あふれは捨てる）
random.seed(7)
for mode in ('wb','rev','none'):
    ok=0; allk=0; steps=set()
    for t in range(100):
        reg={i:(0,0) for i in range(HI)}; acc=0; good=True
        for k in range(20):
            Yv=random.randint(-(M-1)//2,(M-1)//2); yd=todig(Yv,HI)
            reg,co,st=add_step(reg,{i:yd[i] for i in range(HI)},mode); steps.add(st)
            acc=wrap(acc+Yv)
            if val(reg)!=acc: good=False; break
        ok+=good
    print(f"{ {'wb':'書き戻し','rev':'対照1 逆向き','none':'対照2 書き戻しなし'}[mode]}: 20回の足し込みがすべて一致 {ok}/100  段数 {sorted(steps)}")
# 貫く入力：レジスタ 2,1,1,...（=x）に 2,1,1,... を足す
reg={i:(0,0) for i in range(HI)}; xd=[2]+[1]*(HI-1); X=sum(v*5**i for i,v in enumerate(xd))
reg,_,_=add_step(reg,{i:xd[i] for i in range(HI)},'wb'); reg,_,st=add_step(reg,{i:xd[i] for i in range(HI)},'wb')
print("貫く入力（2回目で全桁を繰り上がりが貫く）:", val(reg)==wrap(2*X), " 段数",st)
