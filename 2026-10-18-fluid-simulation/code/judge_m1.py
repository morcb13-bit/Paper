import json, os
from judge_d0 import zm, za, zs, sgn
phi=(1+5**0.5)/2
def load(b,k):
    f='mv_%s_%d.json'%(b,k)
    return json.load(open(f)) if os.path.exists(f) else None
def comp(rec,a): return [(r[2*a],r[2*a+1]) for r in rec]
def windows(T,n=10):
    L=len(T); out=[]
    for i in range(n):
        s=(0,0)
        for t in T[i*L//n:(i+1)*L//n]: s=za(s,t)
        out.append(s)
    return out
def big(W):
    S=(0,0);Q=(0,0)
    for w in W: S=za(S,w); Q=za(Q,zm(w,w))
    return S, sgn(za(zm(S,S),zs(Q,-9)))>0
def f(z): return z[0]+z[1]*phi
rate={}
print('物体 速さ 帳簿 | 抵抗 ΣΔP_x（手番あたり≈） 雑音越え | 上下 ΣΔP_y（手番あたり≈） 雑音越え')
for b,k in [('ball',0),('ball',8),('ball',4),('ball',2),('ball',1),('wing+',2),('wing+',1),('wing-',2),('wing-',1)]:
    r=load(b,k)
    if r is None: print(b,k,'未完'); continue
    X=comp(r['rec'],0); Y=comp(r['rec'],1); n=len(X)
    Sx,bx=big(windows(X)); Sy,by=big(windows(Y)); rate[(b,k)]=(Sx,n)
    sp='止め' if k==0 else '1/%d'%k
    print(f"{b:5s} {sp:4s} {r['ok']} | {f(Sx)/n:+8.1f} {'越え' if bx else '内'} | {f(Sy)/n:+8.1f} {'越え' if by else '内'}")
# 抵抗の単調な増え方（手番あたりで比べる：Sx_a·n_b と Sx_b·n_a を整数で比べる）
order=[('ball',0),('ball',8),('ball',4),('ball',2),('ball',1)]
ok=True
for a,b in zip(order,order[1:]):
    if a in rate and b in rate:
        Sa,na=rate[a]; Sb,nb=rate[b]
        gt=sgn(za(zs(Sb,na),zs(Sa,-nb)))>0
        print('  抵抗',a[1],'→',b[1],'で増える:',gt); ok&=gt
print('抵抗が速さとともに増える（止め<8<4<2<1）:',ok)
