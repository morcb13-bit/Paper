"""検定SC1：星型立方体（立方体を互い違いに二つ重ねた形）の性質
数は Z[√2] の整数組 (a,b)=a+b√2 だけで扱う。割り算は使わず、点は同次座標 (X,Y,Z,W)。
対照：同じ手順を星型八面体（正四面体＋点対称の鏡像）に通し、芯＝正八面体(頂点6・面8)、先端8 が出ること。
出なければ手順が壊れている（NG）。
調べる量：(1)二つの立体の頂点の合計と共有数 (2)芯（共通部分）の頂点・辺・面と面の角数
(3)先端（一方の外にはみ出す部分）の個数と、一つの先端の頂点の数（格子点での連結成分）
"""
from itertools import combinations
class Q:
    __slots__=('a','b')
    def __init__(s,a,b=0): s.a,s.b=a,b
    def __add__(s,o): o=o if isinstance(o,Q) else Q(o); return Q(s.a+o.a,s.b+o.b)
    def __sub__(s,o): o=o if isinstance(o,Q) else Q(o); return Q(s.a-o.a,s.b-o.b)
    def __mul__(s,o): o=o if isinstance(o,Q) else Q(o); return Q(s.a*o.a+2*s.b*o.b,s.a*o.b+s.b*o.a)
    def __neg__(s): return Q(-s.a,-s.b)
    def sign(s):
        a,b=s.a,s.b
        if a>=0 and b>=0: return 0 if a==0 and b==0 else 1
        if a<=0 and b<=0: return -1
        l,r=a*a,2*b*b      # a と b の符号が違う
        if l==r: return 0
        return (1 if a>0 else -1) if l>r else (1 if b>0 else -1)
    def iszero(s): return s.a==0 and s.b==0
def q(v): return v if isinstance(v,Q) else Q(v)
def det3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
           -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
           +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
def dot(n,p): return n[0]*p[0]+n[1]*p[1]+n[2]*p[2]

def solid(planes): return [([q(c) for c in n],q(d)) for n,d in planes]

def vertices(P):
    """面を三つずつ解き、すべての面の内側にある点（同次座標）"""
    out=[]
    for (n1,d1),(n2,d2),(n3,d3) in combinations(P,3):
        M=[n1,n2,n3]; W=det3(M)
        if W.iszero(): continue
        X=det3([[d1,n1[1],n1[2]],[d2,n2[1],n2[2]],[d3,n3[1],n3[2]]])
        Y=det3([[n1[0],d1,n1[2]],[n2[0],d2,n2[2]],[n3[0],d3,n3[2]]])
        Z=det3([[n1[0],n1[1],d1],[n2[0],n2[1],d2],[n3[0],n3[1],d3]])
        if W.sign()<0: X,Y,Z,W=-X,-Y,-Z,-W
        h=(X,Y,Z,W)
        if all((d*W-dot(n,h[:3])).sign()>=0 for n,d in P):
            if not any(same(h,g) for g in out): out.append(h)
    return out
def same(h,g): return all((h[i]*g[3]-g[i]*h[3]).iszero() for i in range(3))
def on(n,d,h): return (dot(n,h[:3])-d*h[3]).iszero()

def key(n,d): return tuple((c.a,c.b) for c in n)+((d.a,d.b),)
def core_stats(P):
    U={}
    for n,d in P: U.setdefault(key(n,d),(n,d))
    P=list(U.values())
    V=vertices(P)
    faces=[]
    for n,d in P:
        k=sum(on(n,d,h) for h in V)
        if k>=3 and not any(all(on(n,d,h)==on(m,e,h) for h in V) for m,e in [f[1] for f in faces]):
            faces.append((k,(n,d)))
    E=sum(1 for h,g in combinations(V,2) if sum(on(n,d,h) and on(n,d,g) for n,d in P)>=2)
    return len(V),E,sorted(f[0] for f in faces),V

NB=[(a,b,c) for a in(-1,0,1) for b in(-1,0,1) for c in(-1,0,1) if (a,b,c)!=(0,0,0)]
def tips(PA,PB,L):
    """A の内側（厳密）かつ B の外側（厳密）の格子点の連結成分（26近傍）"""
    pts=set()
    for i in range(-2*L,2*L+1):
        for j in range(-2*L,2*L+1):
            for k in range(-2*L,2*L+1):
                p=(Q(i),Q(j),Q(k)); W=Q(L)
                if all((d*W-dot(n,p)).sign()>0 for n,d in PA) and any((d*W-dot(n,p)).sign()<0 for n,d in PB):
                    pts.add((i,j,k))
    seen=set(); comps=[]
    for s in pts:
        if s in seen: continue
        st=[s]; seen.add(s); c=0
        while st:
            x,y,z=st.pop(); c+=1
            for dx,dy,dz in NB:
                t=(x+dx,y+dy,z+dz)
                if t in pts and t not in seen: seen.add(t); st.append(t)
        comps.append(c)
    return len(comps)

def run(name,PA,PB,L):
    PA,PB=solid(PA),solid(PB)
    VA,VB=vertices(PA),vertices(PB)
    shared=sum(any(same(h,g) for g in VB) for h in VA)
    nv,ne,fs,CV=core_stats(PA+PB)
    ta,tb=tips(PA,PB,L),tips(PB,PA,L)
    print(f'[{name}] 頂点 A{len(VA)} B{len(VB)} 共有{shared} 合計{len(VA)+len(VB)-shared}')
    print(f'   芯：頂点{nv} 辺{ne} 面{len(fs)} 面の角数{fs}')
    if name.startswith('三回'):
        for h in CV: print('     芯の頂点', [(c.a,c.b) for c in h])
    print(f'   先端：Aのはみ出し{ta} Bのはみ出し{tb} 計{ta+tb}')

r2=Q(0,1)
# 対照：星型八面体
TA=[((1,1,-1),1),((1,-1,1),1),((-1,1,1),1),((-1,-1,-1),1)]
TB=[(tuple(-c for c in n),d) for n,d in TA]
CA=[((1,0,0),1),((-1,0,0),1),((0,1,0),1),((0,-1,0),1),((0,0,1),1),((0,0,-1),1)]
# 三回軸のまわり180度：法線 (−1,2,2)/3 など → −x+2y+2z ≤ 3
CB3=[]
for i in range(3):
    n=[2,2,2]; n[i]=-1
    CB3+= [(tuple(n),3),(tuple(-c for c in n),3)]
# 四回軸のまわり45度：x±y ≤ √2, z ≤ 1
CB4=[((1,1,0),r2),((-1,-1,0),r2),((1,-1,0),r2),((-1,1,0),r2),((0,0,1),1),((0,0,-1),1)]
L=12
run('対照 星型八面体',TA,TB,L)
run('三回軸180度',CA,CB3,L)
run('四回軸45度',CA,CB4,L)
