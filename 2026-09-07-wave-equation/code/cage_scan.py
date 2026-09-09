#  cage_scan.py ── 檻の 2 の冪と 24 の可除性を、剰余を伸ばして数える
#
#  検定S1 担体
#      OK なら：RT面30(次数4)・笠10(次数2と4)・600胞体120(次数12)が組める
#  検定S2 一歩
#      OK なら：T^T T = t^2 I
#      NG なら：以降の位数はすべて無効
#  検定S3 v247 の再現（負の対照の代わり）
#      OK なら：600胞体の T^k=I が m=5,7,11,13 で 4 / 600 / 60 / 840 に一致
#      NG なら：一歩の書き方か担体が v247 と違う。数値を出さない
#  検定S4 24 の可除性と 5 の身分
#      成り立たない例が出れば、その並びは規則ではない

import math, random
exec(open('rt_walk.py').read().split('# ================================================================ 実行')[0])

# ---------------------------------------------------------------- 担体
def rt_face_graph():
    verts,n_ico,edges,_ = rt_graph()
    faces,_ = rt_faces(verts,n_ico,edges)
    fedges={}
    for i,f in enumerate(faces):
        for k in range(4): fedges.setdefault(frozenset((f[k],f[(k+1)%4])),[]).append(i)
    nbr={i:[] for i in range(len(faces))}
    for e,fs in fedges.items():
        if len(fs)==2: nbr[fs[0]].append(fs[1]); nbr[fs[1]].append(fs[0])
    for i in nbr: nbr[i].sort()
    def direction(u,v):
        d=vsub(verts[u],verts[v])
        for c in d:
            s=psign(c)
            if s:
                if s<0: d=tuple(pneg(x) for x in d)
                break
        return d
    dirs={}
    for a,b in edges: dirs.setdefault(direction(a,b),[]).append(frozenset((a,b)))
    zs=[]
    for dv,es in sorted(dirs.items()):
        fs=set()
        for e in es: fs.update(fedges[e])
        zs.append(sorted(fs))
    belt=set(zs[0]); rest=[f for f in range(30) if f not in belt]
    seen=set(); caps=[]
    for f in rest:
        if f in seen: continue
        comp=[f]; seen.add(f); st=[f]
        while st:
            x=st.pop()
            for y in nbr[x]:
                if y in rest and y not in seen: seen.add(y); comp.append(y); st.append(y)
        caps.append(sorted(comp))
    return nbr, caps[0]

def cell600():
    import itertools
    P=[]
    for i in range(4):
        for s in (2,-2):
            v=[PZ]*4; v[i]=(s,0); P.append(tuple(v))
    for a in (1,-1):
        for b in (1,-1):
            for c in (1,-1):
                for d in (1,-1):
                    P.append(((a,0),(b,0),(c,0),(d,0)))
    EVEN=[p for p in itertools.permutations(range(4))
          if sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])%2==0]
    for p in EVEN:
        for s1 in (1,-1):
            for s2 in (1,-1):
                for s3 in (1,-1):
                    b=[PZ,(s1,0),(-s2,s2),(0,s3)]     # 0, ±1, ±(phi-1), ±phi
                    v=[None]*4
                    for i in range(4): v[p[i]]=b[i]
                    P.append(tuple(v))
    P=sorted(set(P))
    def n2(u,v):
        s=PZ
        for a,b in zip(u,v): s=padd(s,pmul(psub(a,b),psub(a,b)))
        return s
    best=None
    for i in range(1,len(P)):
        d=n2(P[0],P[i])
        if best is None or pcmp(d,best)<0: best=d
    nbr={i:[] for i in range(len(P))}
    for i in range(len(P)):
        for j in range(i+1,len(P)):
            if n2(P[i],P[j])==best: nbr[i].append(j); nbr[j].append(i)
    for i in nbr: nbr[i].sort()
    return P, nbr

# ---------------------------------------------------------------- 一歩
def build(nbr, t):
    arcs=[(u,v) for u in sorted(nbr) for v in nbr[u]]
    idx={a:i for i,a in enumerate(arcs)}
    rows=[]
    for (u,v) in arcs:
        d=len(nbr[u]); assert (2*t)%d==0, f"c が整数にならない d={d} t={t}"
        c=(2*t)//d
        r={}
        for w in nbr[u]: r[idx[(w,u)]]=r.get(idx[(w,u)],0)+c
        r[idx[(v,u)]]=r.get(idx[(v,u)],0)-t
        rows.append([(j,val) for j,val in sorted(r.items()) if val])
    return arcs, rows

def check_norm(arcs, rows, t):
    n=len(arcs)
    cols={}
    for i,r in enumerate(rows):
        for j,v in r: cols.setdefault(j,{})[i]=v
    for a in range(n):
        for b in range(a,n):
            ca=cols.get(a,{}); cb=cols.get(b,{})
            s=sum(v*cb[i] for i,v in ca.items() if i in cb)
            if s != (t*t if a==b else 0): return False
    return True

# ---------------------------------------------------------------- 多項式（mod m）
def pdiv_mod(a,b,m):
    a=[x%m for x in a]; db=len(b)-1; inv=pow(b[-1],-1,m)
    q=[0]*max(1,len(a)-db)
    for i in range(len(a)-1,db-1,-1):
        if a[i]%m:
            f=a[i]*inv%m; q[i-db]=f
            for j,bj in enumerate(b): a[i-db+j]=(a[i-db+j]-f*bj)%m
    r=a[:db]
    while len(r)>1 and r[-1]%m==0: r.pop()
    return q,r

def ptrim(a,m):
    a=[x%m for x in a]
    while len(a)>1 and a[-1]==0: a.pop()
    return a if a else [0]

def pgcd_mod(a,b,m):
    a=ptrim(a,m); b=ptrim(b,m)
    if not a or (len(a)==1 and a[0]==0): a=[1] if (not b or (len(b)==1 and b[0]==0)) else b
    if not b or (len(b)==1 and b[0]==0): return [x*pow(a[-1],-1,m)%m for x in a]
    while True:
        _,r=pdiv_mod(a,b,m); r=ptrim(r,m)
        a,b=b,r
        if not b or (len(b)==1 and b[0]%m==0): break
    inv=pow(a[-1],-1,m); return [x*inv%m for x in a]

def pmul_mod(a,b,m):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b): r[i+j]=(r[i+j]+x*y)%m
    return r

def poly_lcm_mod(a,b,m):
    g=pgcd_mod(a,b,m); q,_=pdiv_mod(pmul_mod(a,b,m),g,m)
    inv=pow(q[-1],-1,m); return [x*inv%m for x in q]

def pmulmod(a,b,f,m):
    _,r=pdiv_mod(pmul_mod(a,b,m),f,m)
    r=list(r)
    while len(r)<len(f)-1: r.append(0)
    return r[:len(f)-1]

def ppowmod(a,e,f,m):
    r=[1]+[0]*(len(f)-2)
    a=list(a)
    while len(a)<len(f)-1: a.append(0)
    while e:
        if e&1: r=pmulmod(r,a,f,m)
        a=pmulmod(a,a,f,m); e>>=1
    return r

def xpoly(f):
    x=[0,1]
    while len(x)<len(f)-1: x.append(0)
    return x[:max(1,len(f)-1)]

# ---------------------------------------------------------------- 最小多項式
def apply_rows(rows,x,m):
    return [sum(c*x[j] for j,c in r)%m for r in rows]

def lin_dep(K,m):
    r=len(K); n=len(K[0])
    M=[list(K[i])+[1 if j==i else 0 for j in range(r)] for i in range(r)]
    piv=0
    for c in range(n):
        p=None
        for i in range(piv,r):
            if M[i][c]%m: p=i; break
        if p is None: continue
        M[piv],M[p]=M[p],M[piv]
        inv=pow(M[piv][c],-1,m)
        M[piv]=[x*inv%m for x in M[piv]]
        for i in range(r):
            if i!=piv and M[i][c]%m:
                fq=M[i][c]
                M[i]=[(M[i][j]-fq*M[piv][j])%m for j in range(len(M[i]))]
        piv+=1
        if piv==r: break
    for i in range(r):
        if all(x%m==0 for x in M[i][:n]):
            co=M[i][n:]
            while len(co)>1 and co[-1]%m==0: co.pop()
            if any(x%m for x in co): return [x%m for x in co]
    return None

def minpoly_mod(rows,m,tries=3,cap=90):
    best=[1]
    for _ in range(tries):
        v=[random.randrange(m) for _ in range(len(rows))]
        K=[v]; f=None
        for _ in range(cap):
            v=apply_rows(rows,v,m); K.append(v)
            f=lin_dep(K,m)
            if f: break
        if f is None: raise RuntimeError("最小多項式が出ない")
        best=poly_lcm_mod(best,f,m) if best!=[1] else f
    return best

# ---------------------------------------------------------------- 位数
def ddf_degrees(f,m):
    degs=set(); cur=list(f)
    d=0; h=xpoly(cur)
    while len(cur)>1 and d<len(f):
        d+=1
        h=ppowmod(xpoly(cur),m**d,cur,m)
        t=list(h)
        while len(t)<2: t.append(0)
        t[1]=(t[1]-1)%m
        t=ptrim(t,m)
        if len(t)==1 and t[0]==0:
            degs.add(d); cur=[1]; break
        g=pgcd_mod(cur,t,m)
        if len(g)>1:
            degs.add(d)
            while len(cur)>1:
                q,r=pdiv_mod(cur,g,m)
                if len(r)==1 and r[0]%m==0: cur=q
                else: break
    return degs

def is_prime(n):
    if n<2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n%p==0: return n==p
    d=n-1; s=0
    while d%2==0: d//=2; s+=1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x=pow(a,d,n)
        if x in (1,n-1): continue
        for _ in range(s-1):
            x=x*x%n
            if x==n-1: break
        else: return False
    return True

def rho(n):
    if n%2==0: return 2
    while True:
        x=random.randrange(2,n); y=x; c=random.randrange(1,n); d=1
        while d==1:
            x=(x*x+c)%n; y=(y*y+c)%n; y=(y*y+c)%n
            d=math.gcd(abs(x-y),n)
        if d!=n: return d

def fact(n,out=None):
    if out is None: out={}
    if n==1: return out
    if is_prime(n): out[n]=out.get(n,0)+1; return out
    d=rho(n); fact(d,out); fact(n//d,out); return out

def order_mod(f,m):
    degs=ddf_degrees(f,m) or {len(f)-1}
    L=1
    for d in degs:
        v=m**d-1; L=L*v//math.gcd(L,v)
    e=0
    while m**e < len(f)-1: e+=1
    L*=m**e
    one=[1]+[0]*(len(f)-2)
    if ppowmod(xpoly(f),L,f,m)!=one: return None
    k=L
    for p in sorted(fact(L)):
        while k%p==0 and ppowmod(xpoly(f),k//p,f,m)==one: k//=p
    return k

# ================================================================ 実行
random.seed(13)
nbr_rt, cap = rt_face_graph()
S=set(cap); nbr_cap={f:[g for g in nbr_rt[f] if g in S] for f in cap}
P600, nbr600 = cell600()

print("検定S1 担体")
print(f"  RT面    {len(nbr_rt)}面   次数{sorted(set(len(v) for v in nbr_rt.values()))} 辺{sum(len(v) for v in nbr_rt.values())//2}")
print(f"  笠      {len(nbr_cap)}面   次数{sorted(set(len(v) for v in nbr_cap.values()))} 辺{sum(len(v) for v in nbr_cap.values())//2}")
print(f"  600胞体 {len(nbr600)}頂点 次数{sorted(set(len(v) for v in nbr600.values()))} 辺{sum(len(v) for v in nbr600.values())//2}")
ok1=(len(nbr_rt)==30 and len(nbr_cap)==10 and len(nbr600)==120
     and set(len(v) for v in nbr600.values())=={12})
print(f"  → 検定S1 {'OK' if ok1 else 'NG'}")

CAR={}
for name,nb,t in (("RT面30",nbr_rt,2),("笠10",nbr_cap,2),("600胞体",nbr600,6)):
    CAR[name]=build(nb,t)+(t,)

print("\n検定S2 一歩（T^T T = t² I）")
ok2=True
for name in ("RT面30","笠10"):
    arcs,rows,t=CAR[name]; r=check_norm(arcs,rows,t); ok2&=r
    print(f"  {name:<8} 有向辺{len(arcs):>5} t={t}  {r}")
arcs,rows,t=CAR["600胞体"]
print(f"  {'600胞体':<8} 有向辺{len(arcs):>5} t={t}  （1440次元。検定S3で代える）")
print(f"  → 検定S2 {'OK' if ok2 else 'NG'}")

print("\n検定S3 v247 の再現（600胞体・t=6）")
exp={5:4,7:600,11:60,13:840}; ok3=True
for m,e in exp.items():
    k=order_mod(minpoly_mod(rows,m),m); ok3&=(k==e)
    print(f"  m={m:<3} k={k}   v247は{e}   {'一致' if k==e else '不一致'}")
print(f"  → 検定S3 {'OK' if ok3 else 'NG'}")

def cls(m): return "分岐" if m==5 else ("分解" if m%5 in (1,4) else "惰性")

print("\n検定S4 檻の2の冪と 24")
MS=[5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61]
verdict={}
for name in ("RT面30","笠10","600胞体"):
    arcs,rows,t=CAR[name]
    print(f"\n  【{name}】 t={t}")
    print(f"  {'m':>3} {'5':<4} {'檻':>28} {'2冪':>5} {'3冪':>4} {'24':>4}")
    for m in MS:
        k=order_mod(minpoly_mod(rows,m),m)
        if k is None:
            print(f"  {m:>3} {cls(m):<4} {'—':>28}"); continue
        v2=0; kk=k
        while kk%2==0: v2+=1; kk//=2
        v3=0; kk=k
        while kk%3==0: v3+=1; kk//=3
        verdict.setdefault(cls(m),[]).append((name,m,k%24==0))
        print(f"  {m:>3} {cls(m):<4} {k:>28} {v2:>5} {v3:>4} {'割る' if k%24==0 else '—':>4}")

print("\n  5の身分ごとの 24 の可除性")
for c in ("分岐","分解","惰性"):
    rows_=verdict.get(c,[])
    yes=sum(1 for _,_,b in rows_ if b)
    print(f"    {c}: {yes}/{len(rows_)} で 24 が割る")
