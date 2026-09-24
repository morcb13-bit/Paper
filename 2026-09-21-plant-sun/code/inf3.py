# 検定：円環の中心の頂点の型が、置き換えで育てたひし形のペンローズの頂点の型に入るか
# 型＝連続接続（辺）で繋がる隣の向き（36°刻み）の間隔の列。回転・反転は同一視
import b13_chain_units as U, inf1 as I, collections
def canon(ks):
    ks=sorted(ks); g=[(ks[(i+1)%len(ks)]-ks[i])%10 or 10 for i in range(len(ks))]
    c=[tuple(g[i:]+g[:i]) for i in range(len(g))]; c+= [tuple(reversed(x)) for x in c]
    return min(c)
def dirs(vec, units):                          # vec が units[k] に一致する k
    return units.index(vec) if vec in units else None
# 1) 置き換えで育てたひし形の頂点の型（辺＝長さ1の辺、向き ζ^k）
T=I.wheel
for _ in range(9): T=I.step(T)
nb=collections.defaultdict(set); UN=[U.zt(k) for k in range(10)]
for _,A,B,C in T:
    for x,y in ((A,B),(A,C),(B,C)):
        v=U.zsub(y,x)
        if v in UN: nb[x].add(UN.index(v)); nb[y].add(UN.index(U.zsub(x,y)))
def r2(z): x,y=U.xy(z); return x*x+y*y
Rmax=max(r2(v) for v in nb)
typesP=collections.Counter(canon(ks) for v,ks in nb.items() if 0.02*Rmax< r2(v) <0.36*Rmax)   # 中心と縁を除く
ctr=canon(nb[U.ZERO])
print("ひし形の頂点の型（中心と縁を除く、9段）",len(typesP),"種")
for t,n in typesP.most_common(): print("  ",t,n)
print("  （種の中心の型",ctr,"）")
TP=set(typesP)
# 2) 円環の中心
exec(open('gen10.py').read().split("pid={}")[0])
def judge(pts, units, name):
    S=set(pts); z0=(2,-2,0,-3)
    d=lambda z:r2(U.zsub(z,z0)); M=max(d(p) for p in pts)
    inn=[p for p in pts if d(p)<0.36*M]
    ty=collections.Counter(); bad=collections.Counter()
    for p in inn:
        ks=[k for k,u in enumerate(units) if U.zadd(p,u) in S]
        if len(ks)<2: bad[('孤立',len(ks))]+=1; continue
        t=canon(ks); ty[t]+=1
        if t not in TP: bad[t]+=1
    ok=sum(n for t,n in ty.items() if t in TP)
    print(f"{name}: 内側 {len(inn)} 点のうち ひし形の頂点の型に入る {ok}")
    for t,n in bad.most_common(8): print("   外れ",t,n)
    return ok,len(inn)
RC=list({place(c,k) for (w,r,c,k) in rings})
judge(RC, U.CONT, "円環の中心（辺＝CONT）")
cells=U.fits([place(c,k) for (w,r,c,k) in rings])
NB=[U.zadd(U.ZERO,u) for u in [U.zmul(U.PHI,U.zt(k)) for k in range(10)]]
judge(list(cells), NB, "負の対照：五角形の中心（辺＝隣の五角形 φ）")
