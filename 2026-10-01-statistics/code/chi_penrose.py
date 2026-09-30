# ペンローズの床で「χ二乗検定」を離散の操作で書き直す
# 問い：丘の高さは、番地の景色（隣の数）で決まるか
# 着地：源から T 刻の「留まる／隣へ移る」道を一本ずつ引いた電子の着地番地
# 型紙（H0）：同じ距離の輪の中では、どの番地にも同じだけ着地する（景色は効かない）
# 残渣：組ごとに O×n_輪 − N_輪×n_組（整数）。統計量は残渣の絶対値の和
# 判定：着地はそのままに、ラベルを同じ輪の中で K 回入れ替え、実測以上の統計量が出た回数 k を数える
# 負の対照：隣の数のラベルを同じ輪の中で入れ替えたもの（組の大きさは同じ、景色とは無関係）で同じことをする → k が小さくならなければ OK
import math, random, sys
PHI=(1+5**0.5)/2
def penrose(n):
    tris=[]
    for i in range(10):
        b=(math.cos((2*i-1)*math.pi/10),math.sin((2*i-1)*math.pi/10)); c=(math.cos((2*i+1)*math.pi/10),math.sin((2*i+1)*math.pi/10))
        if i%2==0: b,c=c,b
        tris.append((0,(0,0),b,c))
    for _ in range(n):
        out=[]
        for col,A,B,C in tris:
            if col==0:
                P=(A[0]+(B[0]-A[0])/PHI,A[1]+(B[1]-A[1])/PHI); out+= [(0,C,P,B),(1,P,C,A)]
            else:
                Q=(B[0]+(A[0]-B[0])/PHI,B[1]+(A[1]-B[1])/PHI); R=(B[0]+(C[0]-B[0])/PHI,B[1]+(C[1]-B[1])/PHI)
                out+= [(1,R,C,A),(1,Q,R,B),(0,R,Q,A)]
        tris=out
    return tris
tris=penrose(8); e=tris[0]; s=1/math.hypot(e[2][0]-e[1][0],e[2][1]-e[1][1])
WH=min(34,int(s*0.72)); HH=round(WH*0.6)
idx={}; X=[]; Y=[]; edges=set()
def vid(p):
    k=(round(p[0]*s*1e3),round(p[1]*s*1e3))
    if k not in idx: idx[k]=len(X); X.append(p[0]*s); Y.append(p[1]*s)
    return idx[k]
inR=lambda p: abs(p[0]*s)<=WH and abs(p[1]*s)<=HH
for _,A,B,C in tris:
    for Q in (B,C):
        if inR(A) and inR(Q):
            a,b=vid(A),vid(Q); edges.add((min(a,b),max(a,b)))
n=len(X); adj=[[] for _ in range(n)]
for a,b in edges: adj[a].append(b); adj[b].append(a)
src=min(range(n),key=lambda v:(X[v]-5.3)**2+(Y[v]-2.1)**2)
T=int(sys.argv[1]) if len(sys.argv)>1 else 30
NE=int(sys.argv[2]) if len(sys.argv)>2 else 20000
K=200
# 残りの刻 k で番地 u から始まる道の数 R[k][u]（整数の足し算）
R=[[1]*n]
for k in range(1,T+1):
    prev=R[-1]; R.append([prev[u]+sum(prev[w] for w in adj[u]) for u in range(n)])
rng=random.Random(13)
def one():
    u=src
    for k in range(T,0,-1):
        r=rng.randrange(R[k][u]); prev=R[k-1]
        for w in [u]+adj[u]:
            if r<prev[w]: u=w; break
            r-=prev[w]
    return u
land=[0]*n
for _ in range(NE): land[one()]+=1
dist=[math.hypot(X[v]-X[src],Y[v]-Y[src]) for v in range(n)]
shell=[int(d) for d in dist]; RMAX=10
def stat(landv,label):
    tot=0; cells={}
    for v in range(n):
        if shell[v]>RMAX: continue
        c=cells.setdefault((shell[v],label[v]),[0,0]); c[0]+=landv[v]; c[1]+=1
    sh={}
    for (a,l),(o,m) in cells.items():
        t=sh.setdefault(a,[0,0]); t[0]+=o; t[1]+=m
    res={}
    for (a,l),(o,m) in cells.items():
        No,ns=sh[a]; r=o*ns-No*m; tot+=abs(r); res[(a,l)]=r
    return tot,res,cells,sh
inner=[v for v in range(n) if shell[v]<=RMAX]
byshell={}
for v in inner: byshell.setdefault(shell[v],[]).append(v)
def null_land():
    L=[0]*n
    for a,vs in byshell.items():
        No=sum(land[v] for v in vs)
        for _ in range(No): L[rng.choice(vs)]+=1
    return L
deg=[len(adj[v]) for v in range(n)]
shuf=deg[:]; r2=random.Random(7)
for a,vs in byshell.items():
    ls=[deg[v] for v in vs]; r2.shuffle(ls)
    for v,l in zip(vs,ls): shuf[v]=l
labels={'隣の数':deg, '隣の数を輪の中で入れ替えたラベル（負の対照）':shuf}

print(f'刻T={T} 電子={NE} 輪の半径≤{RMAX} 番地={len(inner)} 対照K={K}')
for name,lab in labels.items():
    S,res,cells,sh=stat(land,lab)
    k=0; rr=random.Random(99)
    for _ in range(K):
        pl=lab[:]
        for a,vs in byshell.items():
            ls=[lab[v] for v in vs]; rr.shuffle(ls)
            for v,l in zip(vs,ls): pl[v]=l
        if stat(land,pl)[0]>=S: k+=1
    # 観察者側：古典のχ二乗
    chi=0.0; dof=0
    for (a,l),(o,m) in cells.items():
        No,ns=sh[a]; E=No*m/ns
        if E>0: chi+=(o-E)**2/E; dof+=1
    dof-=len(sh)
    print(f'\n[{name}] 統計量（残渣の絶対値の和）={S}  対照で実測以上={k}/{K}  （観察者側 χ²={chi:.1f}, 自由度={dof}）')
    if name=='隣の数':
        agg={}
        for (a,l),r in res.items():
            if a>=1: agg[l]=agg.get(l,0)+r
        # 組ごとの O と E（輪をまたいで合計、E は整数比で表示用）
        OE={}
        for (a,l),(o,m) in cells.items():
            if a<1: continue
            No,ns=sh[a]; t=OE.setdefault(l,[0,0.0]); t[0]+=o; t[1]+=No*m/ns
        for l in sorted(OE): print(f'  隣{l}: 着地={OE[l][0]}  型紙の個数={OE[l][1]:.0f}  残渣の符号={"+" if agg[l]>0 else "−"}')
