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

# 一桁：A＝アキュムレータの位相の組、B＝語から鏡を通して渡された位相の組
_C={}
def dig(i,k,m,kb,mb,cin):
    key=(i,k,m,kb,mb,cin)
    if key in _C: return _C[key]
    g=G[i]; ka=(k-2*cin)%10; ma=(m-2*cin)%4
    L=[q for q in g['scr'] if rep(ZR,g['PA'][q],ka)==rep(ZR,g['PB'][q],kb)]
    lv={(g['dA'][q]-g['dB'][q])//2 for q in L}
    if len(lv)!=1: _C[key]=None; return None
    s=lv.pop(); same={rep(TW,g['SA'][q],ma)==rep(TW,g['SB'][q],mb) for q in L}
    if len(same)!=1: _C[key]=None; return None
    cout=0 if same.pop() else ((-1 if s>0 else 1) if s!=0 else ((1 if cin>0 else -1) if cin else 0))
    q=L[0]
    kk=[e for e in range(10) if rep(ZR,g['PA'][q],e)==g['PB'][q]][0]
    mm=[e for e in range(4) if rep(TW,g['SA'][q],e)==g['SB'][q]][0]
    _C[key]=(cout,(kk,mm)); return _C[key]
class Bad(Exception): pass
def D(*a):
    r=dig(*a)
    if r is None: raise Bad()
    return r
def add_step(acc,opB):
    rng=range(-LO,HI); R=dict(acc); B=dict(opB)
    for i in range(-LO,0): R[i]=(0,0); B[i]=(0,0)
    f={i:tuple(D(i,*R[i],*B[i],c)[0] for c in (-1,0,1)) for i in rng}
    A=dict(f); P=dict(f); ra=rb=PHI3i; st=0
    while not all(len(set(v))==1 for v in A.values()):
        nA={}
        for i in rng:
            j=IDX.get(U.zadd(C2,U.zmul(ra,U.zsub(REF[i],C2))))
            q=P[j] if j is not None else (-1,0,1)
            nA[i]=tuple(A[i][q[c]+1] for c in range(3))
        P=A; A=nA; ra,rb=U.zmul(ra,rb),ra; st+=1
        if st>30: raise Bad()
    cin={i:(A[i-1][1] if i>0 else A[-1][1]) for i in range(HI)}
    return {i:D(i,*R[i],*B[i],cin[i])[1] for i in range(HI)}
PH={((-2*r)%10,(-2*r)%4):r for r in range(-2,3)}
def phase(r): return ((-2*r)%10,(-2*r)%4)
def mirror(p): return ((-p[0])%10,(-p[1])%4)
def val(reg): return sum(PH[reg[i]]*5**i for i in range(HI))
def todig(X,n):
    out=[]
    while X: r=((X+2)%5)-2; out.append(r); X=(X-r)//5
    return (out+[0]*n)[:n]
M=5**HI
def wrap(X): return ((X+(M-1)//2)%M)-(M-1)//2
# ---- メモリ：語 (u1,u0) の10桁の置き場＝各桁の目印 REF[i] を u0·E0+u1·E1 だけ動かした点 ----
E0=U.zmul(upw(20),D0); E1=U.zmul(U.zt(1),U.zmul(upw(21),D0))
def smul(n,z):
    r=U.ZERO
    for _ in range(abs(n)): r=U.zadd(r,z) if n>0 else U.zsub(r,z)
    return r
def addr_pts(u1,u0,shift=0):
    T=U.zadd(smul(u0+shift,E0),smul(u1,E1))
    return [U.zadd(REF[i],T) for i in range(HI)]
ADDR=[(u1,u0) for u1 in range(-2,3) for u0 in range(-2,3)]
MEMPTS=set(p for a in ADDR for p in addr_pts(*a)); assert len(MEMPTS)==25*HI
def run(mode):
    ok=0
    for t in range(100):
        MEM={}; model={}
        for a in ADDR:
            v=random.randint(-(M-1)//2,(M-1)//2); d=todig(v,HI); model[a]=v
            for p,r in zip(addr_pts(*a),d): MEM[p]=phase(r)
        acc={i:(0,0) for i in range(HI)}; macc=0; good=True
        try:
            for step in range(20):
                op=random.choice(('LOAD','STORE','ADD')); a=random.choice(ADDR)
                pts=addr_pts(*a,shift=(1 if mode=='addr' else 0))
                if any(p not in MEM for p in pts): raise Bad()
                if op=='LOAD':  acc={i:MEM[pts[i]] for i in range(HI)}; macc=model[a]
                if op=='STORE':
                    for i in range(HI): MEM[pts[i]]=acc[i]
                    model[a]=macc
                if op=='ADD':
                    opB={i:(MEM[pts[i]] if mode=='nomirror' else mirror(MEM[pts[i]])) for i in range(HI)}
                    acc=add_step(acc,opB); macc=wrap(macc+model[a])
                if val(acc)!=macc or any(val({i:MEM[p] for i,p in enumerate(addr_pts(*b))})!=model[b] for b in ADDR): good=False; break
        except Bad: good=False
        ok+=good
    return ok
random.seed(11)
print("メモリ25語（番地＝平衡5進二桁）、20手×100通り")
print("  本番                     毎手 アキュムレータと25語が照合と一致", run('ok'),"/100")
print("  対照1 番地の点を E0 一つずらす", run('addr'),"/100")
print("  対照2 読み出しで鏡を通さない", run('nomirror'),"/100")
