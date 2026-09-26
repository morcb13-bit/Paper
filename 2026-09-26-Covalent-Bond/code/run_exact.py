import helix_check as h, sympy as sp
R=lambda M: sp.Matrix([sp.radsimp(sp.expand(x)) for x in M])
def step(T,sg):
    p,q,r,t=T; return [r,t,R(h.rot(p,r,t,sg)),R(h.rot(q,r,t,sg))]
def trip(cs):
    a,b,c=cs[-3]-cs[-4],cs[-2]-cs[-3],cs[-1]-cs[-2]
    return sp.sign(sp.radsimp(sp.expand(a.dot(b.cross(c)))))
def chain(N,mode):
    T=h.V[:]; cs=[h.cen(T)]; signs=[]
    for k in range(N):
        if mode=="screw" or len(cs)<4: sg=1
        else:  # 物理的な向きを一歩ごとに逆へ（負の対照）
            want=-trip(cs); sg=None
            for s in (1,-1):
                T2=step(T,s)
                if trip(cs+[h.cen(T2)])==want: sg=s;break
        T=step(T,sg); cs.append(h.cen(T)); signs.append(sg)
    return cs,signs
def measure(cs,label):
    print("==",label); ok=True; res=[]
    for n in range(1,len(cs)):
        D=cs[n]-cs[0]; d2=sp.radsimp(sp.expand(D.dot(D)))
        if not d2.is_Rational: print(n,"距離²が有理数でない",d2); ok=False; continue
        E=25*2**n*d2
        if not E.is_Integer: print(n,"E非整数",E); ok=False; continue
        E=int(E); dd=E-60*n*n*2**n
        fam=5 if n%2 else 15
        good=dd%6==0 and dd//6>0 and (dd//6)%fam==0 and sp.integer_nthroot((dd//6)//fam,2)[1] and (dd//6)%2==1
        s=sp.integer_nthroot((dd//6)//fam,2)[0] if good else None
        ok&=bool(good); res.append((n,d2,E,dd//6 if dd%6==0 else dd/6,s))
        print(f"n={n:2d} 距離²={d2}  E={E}  d={res[-1][3]}  "+(f"{fam}系×{s}²" if good else "NG"))
    return ok,res
cs,_=chain(20,"screw"); ok,res=measure(cs,"本体（ねじ一定）")
pred=[5,15,5,15,125,135,5,735]
print("E1,E2=",res[0][2],res[1][2]," d予言一致:",[r[3] for r in res[:8]]==pred)
print("平方根:",[r[4] for r in res])
print("本体 合格:",ok)


