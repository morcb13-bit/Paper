# 一本道の検定（事前に決めた基準）
# 合格: E1=150, E2=1050, d_n の予言列と一致, 奇数歩=5系・偶数歩=15系, d_n 奇数 (n<=N)
# 負の対照: 回す向きを ±120° で交互 → 二系統の分岐が崩れること
import sympy as sp, sys
from sympy import sqrt, Rational, Matrix
N=int(sys.argv[1]) if len(sys.argv)>1 else 20
V=[Matrix(v) for v in [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]]  # 中心→辺の距離1
def rot(p,a,b,sgn):
    u=(b-a); u=u/sqrt(u.dot(u)); c=Rational(-1,2); s=sgn*sqrt(3)/2
    x=p-a
    return a + x*c + u.cross(x)*s + u*(u.dot(x))*(1-c)
def cen(T): return sum(T,Matrix([0,0,0]))/4
def issq(k):
    if k<=0: return None
    r=sp.integer_nthroot(k,2); return r[0] if r[1] else None
def run(signs):
    # T=[p,q,r,t]: (p,q)=前と共有した辺, (r,t)=向かい側の辺。向かい側の辺 r→t のまわりに回す
    T=V[:]; cs=[cen(T)]
    for sg in signs:
        p,q,r,t=T
        p2,q2=rot(p,r,t,sg),rot(q,r,t,sg)
        T=[r,t,p2,q2]   # 新しい共有辺は (r,t)、向かい側は回した (p,q)
        T=[sp.simplify(x) for x in T]; cs.append(cen(T))
    out=[]
    for n in range(1,len(cs)):
        D=cs[n]-cs[0]; d2=sp.nsimplify(sp.simplify(D.dot(D)))
        E=sp.simplify(25*2**n*d2)
        out.append((n,d2,E))
    return out
def judge(out,label):
    print("==",label)
    ok=True
    for n,d2,E in out:
        if not (E.is_Integer):
            print(n,"E 非整数",E); ok=False; continue
        E=int(E); dd=E-60*n*n*2**n
        if dd%6: print(n,"E",E,"(E-60n^2 2^n) が6で割れない"); ok=False; continue
        d=dd//6; fam=5 if n%2 else 15
        s=issq(d//fam) if d>0 and d%fam==0 else None
        flag = (s is not None) and d%2==1
        ok &= flag
        print(f"n={n:2d} 距離²={d2}  E={E}  d={d}  {fam}系×{s}²" if s else f"n={n:2d} 距離²={d2}  E={E}  d={d}  NG")
    return ok
if __name__=="__main__":
    for sg in (1,-1):
        o=run([sg]*3)
        print("向き",sg,"一歩の距離²",o[0][1],"二歩",o[1][1])
