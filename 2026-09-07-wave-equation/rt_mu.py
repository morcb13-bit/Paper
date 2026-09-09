#  歩きの行列 D^{-1}A の側（v247 §1-3 と同じ量）を担体ごとに出す
exec(open('rt_walk.py').read().split('# ================================================================ 実行')[0])
from fractions import Fraction
from collections import Counter
verts,n_ico,edges,elen2=rt_graph(); faces,vadj=rt_faces(verts,n_ico,edges)
fedges={}
for i,f in enumerate(faces):
    for k in range(4): fedges.setdefault(frozenset((f[k],f[(k+1)%4])),[]).append(i)
nbr={i:[] for i in range(len(faces))}
for e,fs in fedges.items():
    if len(fs)==2: nbr[fs[0]].append(fs[1]); nbr[fs[1]].append(fs[0])
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
zones=[]
for dv,es in dirs.items():
    fs=set()
    for e in es: fs.update(fedges[e])
    zones.append(sorted(fs))
z=set(zones[0]); restf=[f for f in range(30) if f not in z]
seen=set(); caps=[]
for f in restf:
    if f in seen: continue
    comp=[f]; seen.add(f); st=[f]
    while st:
        x=st.pop()
        for y in nbr[x]:
            if y in restf and y not in seen: seen.add(y); comp.append(y); st.append(y)
    caps.append(sorted(comp))

def charpoly_frac(A):
    n=len(A); I=[[Fraction(1 if i==j else 0) for j in range(n)] for i in range(n)]
    Mprev=I; cs=[Fraction(1)]
    for k in range(1,n+1):
        AM=[[sum(A[i][t]*Mprev[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        ck=Fraction(-sum(AM[i][i] for i in range(n)),k)
        Mprev=[[AM[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]
        cs.append(ck)
    return list(reversed(cs))

def clear(p):
    from math import lcm
    L=1
    for c in p: L=lcm(L,c.denominator)
    return [int(c*L) for c in p], L

def show(name, sub):
    S=set(sub); nb={f:[g for g in nbr[f] if g in S] for f in sub}
    idx={f:i for i,f in enumerate(sub)}
    n=len(sub)
    P=[[Fraction(0) for _ in range(n)] for _ in range(n)]
    for f in sub:
        d=len(nb[f])
        for g in nb[f]: P[idx[f]][idx[g]]=Fraction(1,d)
    cp=charpoly_frac(P)
    ip,L=clear(cp)
    print(f"\n【{name}】 面 {n}  次数 {sorted(Counter(len(nb[f]) for f in sub).items())}")
    # 因数分解（有理根＋2次）を分母つきで
    fs=[]; cur=ip[:]
    # 有理根 p/q
    for num in range(-8,9):
        for den in (1,2,3,4,6,8,12):
            from math import gcd
            if gcd(abs(num),den)!=1 and not (num==0 and den==1): continue
            g=[-num,den]
            while len(cur)>1:
                q,rem=poly_divmod(cur,g)
                if q is not None and len(rem)==1 and rem[0]==0: fs.append(g); cur=q
                else: break
    for a in range(1,40):
        for b in range(-40,41):
            for c in range(-40,41):
                g=[c,b,a]
                while len(cur)>2:
                    q,rem=poly_divmod(cur,g)
                    if q is not None and len(rem)==1 and rem[0]==0: fs.append(g); cur=q
                    else: break
    cnt=Counter(tuple(x) for x in fs)
    for f,k in sorted(cnt.items(), key=lambda kv:(len(kv[0]),kv[0])):
        f=list(f)
        s=""
        if len(f)==3:
            disc=f[1]*f[1]-4*f[2]*f[0]
            sq=1; d=disc
            i=2
            while i*i<=abs(d):
                while d% (i*i)==0: d//=i*i; sq*=i
                i+=1
            s=f"   判別式 {disc} = {sq*sq}·{d}" if sq>1 else f"   判別式 {disc}"
        print(f"  {pstr(f):<26} 重複 {k}{s}")
    if len(cur)>1: print(f"  残り {pstr(cur)}  （{len(cur)-1}次）")

show("RT 全体（30面・閉じている）", list(range(30)))
show("笠（10面・壁あり）", caps[0])
show("帯（10面・輪）", sorted(z))
