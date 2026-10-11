# 正四面体 n 枚を三回軸のまわりに 120/n 度ずつずらして同じ中心に重ねる
# 型(i)：すべて同じ向き  型(ii)：一つおきに上下反転（n=2 で星型八面体）
# 40度などは整数格子に乗らないので、同一判定だけ50桁の数で行う（観察用）
from mpmath import mp, mpf, cos, sin, pi, matrix, lu_solve
from itertools import combinations
mp.dps=50; EPS=mpf(10)**-30
s3=mp.sqrt(3); s2=mp.sqrt(2)
# 三回軸を z に立てた正四面体（中心 0、上の頂点 (0,0,3)）
base=[(0,0,mpf(3))]+[(2*s2*cos(2*pi*j/3),2*s2*sin(2*pi*j/3),mpf(-1)) for j in range(3)]
def rot(p,th,flip):
    x,y,z=p; c,s=cos(th),sin(th)
    return (c*x-s*y,s*x+c*y,-z if flip else z)
def tets(n,kind):
    return [[rot(p,2*pi*k/(3*n),kind=='ii' and k%2==1) for p in base] for k in range(n)]
def sub(a,b): return tuple(a[i]-b[i] for i in range(3))
def dot(a,b): return sum(a[i]*b[i] for i in range(3))
def crs(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def planes(T):
    P=[]
    for f in combinations(range(4),3):
        a,b,c=(T[i] for i in f); o=T[[i for i in range(4) if i not in f][0]]
        n_=crs(sub(b,a),sub(c,a)); d=dot(n_,a)
        if dot(n_,o)>d: n_=tuple(-x for x in n_); d=-d
        P.append((n_,d,(a,b,c)))
    return P
class PS:
    def __init__(s): s.p=[]
    def key(s,q):
        for i,r in enumerate(s.p):
            if all(abs(q[k]-r[k])<EPS for k in range(3)): return i
        s.p.append(q); return len(s.p)-1
def in_tri(q,a,b,c):
    n_=crs(sub(b,a),sub(c,a))
    return all(dot(crs(sub(v,u),sub(q,u)),n_)>=-EPS for u,v in((a,b),(b,c),(c,a)))
def run(n,kind):
    Ts=tets(n,kind); ps=PS()
    for T in Ts:
        for p in T: ps.key(p)
    nv=len(ps.p)
    # 辺が他の正四面体の面を突き抜ける点（辺の内部）
    pieces=[]; cross=PS()
    for ti,T in enumerate(Ts):
        for a,b in combinations(T,2):
            cut=PS()
            for tj,U in enumerate(Ts):
                if tj==ti: continue
                for n_,d,(p,q,r) in planes(U):
                    den=dot(n_,sub(b,a))
                    if abs(den)<EPS: continue
                    t=(d-dot(n_,a))/den
                    if EPS<t<1-EPS:
                        x=tuple(a[i]+t*(b[i]-a[i]) for i in range(3))
                        if in_tri(x,p,q,r): cut.key(x); cross.key(x)
            pieces.append(len(cut.p)+1)
    # 芯（すべての正四面体の共通部分）
    H=[(n_,d) for T in Ts for n_,d,_ in planes(T)]
    U={}
    for n_,d in H:
        m=mp.sqrt(dot(n_,n_)); k=tuple(mp.nstr(x/m,20) for x in n_)+(mp.nstr(d/m,20),); U[k]=(n_,d)
    H=list(U.values()); cv=PS()
    for (n1,d1),(n2,d2),(n3,d3) in combinations(H,3):
        M=matrix([list(n1),list(n2),list(n3)])
        if abs(mp.det(M))<EPS: continue
        x=lu_solve(M,matrix([d1,d2,d3])); x=(x[0],x[1],x[2])
        if all(dot(n_,x)<=d+EPS for n_,d in H): cv.key(x)
    fc=[sum(abs(dot(n_,x)-d)<EPS for x in cv.p) for n_,d in H]
    fc=sorted(k for k in fc if k>=3)
    return nv,len(cross.p),sorted(set(pieces)),len(cv.p),len(fc),fc
for kind in('i','ii'):
    print('型',kind,'（枚数・角度・頂点・突き抜け点・一辺の切れ数・芯の頂点・芯の面・芯の面の角数）')
    for n in range(1,6):
        nv,nx,pc,cvn,cf,fc=run(n,kind)
        print(' ',n,f'{120/n:g}度',nv,nx,pc,cvn,cf,fc)
