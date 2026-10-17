# 教科書のダイヤモンド：単位格子（辺 12）。原子＝面心立方の点＋その (3,3,3) ずらし。
# 教科書の格子を x→12−x で鏡に映した向きにとる（内部の原子＝面心立方の点＋(−3,3,3)）。ダイヤモンドは鏡に映してもダイヤモンド
from itertools import product, combinations
B=list(product([0,1],repeat=3)); odd=lambda n:sum(n)%2==1
def tet(o,e,w): return [tuple(o[i]+b[i]*e for i in range(3)) for b in B if odd(b)==w]
def atoms(rng):
    S=set()
    for c in product(rng,repeat=3):
        base=[12*c[0],12*c[1],12*c[2]]
        for f in [(0,0,0),(0,6,6),(6,0,6),(6,6,0)]:
            p=tuple(base[i]+f[i] for i in range(3)); S.add(p)
            S.add((p[0]-3,p[1]+3,p[2]+3))
    return S
inside=lambda p:all(0<=x<=12 for x in p)
A=atoms(range(-2,3)); cell=[p for p in A if inside(p)]
d2=lambda p,q:sum((p[i]-q[i])**2 for i in range(3))
interior=[p for p in cell if all(0<x<12 for x in p) and all(x%6==3 for x in p)]
print("格子の中の原子（境界を含む）",len(cell),"内部",sorted(interior))
print("内部原子の小立方体の番地",[tuple(x//6 for x in p) for p in interior],"奇か",[odd(tuple(x//6 for x in p)) for p in interior])
# 内部原子の結合の先＝A の小さな正四面体の頂点か
SUB=[tet((6*n[0],6*n[1],6*n[2]),6,True) for n in B if odd(n)]
for p in interior:
    nb=sorted(q for q in A if d2(p,q)==27)
    T=[T for T in SUB if all(abs(v[i]-p[i])==3 for v in T for i in range(3))][0]
    print(" ",p,"結合の先",nb,"＝A の正四面体の頂点",sorted(nb)==sorted(T))
V=set(v for T in SUB for v in T); print("シェルピンスキー2段目の頂点",len(V),"のうち面の中心",sum(1 for v in V if sorted(x in (0,12) for x in v).count(True)==1),"角",sum(1 for v in V if all(x in(0,12) for x in v)))
print("角 8 のうちシェルピンスキーに入らない角",[c for c in product([0,12],repeat=3) if c not in V],"原子か",[c in A for c in product([0,12],repeat=3) if c not in V])
for c in [c for c in product([0,12],repeat=3) if c not in V]:
    print("  角",c,"の隣で格子の中",[q for q in A if d2(c,q)==27 and inside(q)])
# 結晶全体：どの原子も結合 4 本、最近接の長さ² 27、結合どうしの内積
core=[p for p in A if all(-12<=x<=24 for x in p)]
print("結合の本数（中央付近の原子）",set(sum(1 for q in A if d2(p,q)==27) for p in core))
p=interior[0]; nb=[q for q in A if d2(p,q)==27]; u=[tuple(q[i]-p[i] for i in range(3)) for q in nb]
print("結合の向き",u,"内積",sorted(set(sum(a[i]*b[i] for i in range(3)) for a,b in combinations(u,2))))
# 頂点側の原子（面の中心）の 4 隣は鏡像 B の正四面体
for v in sorted(V):
    if sorted(x in (0,12) for x in v).count(True)==1:
        nb=sorted(q for q in A if d2(v,q)==27); vec=[tuple(q[i]-v[i] for i in range(3)) for q in nb]
        par=[ (sum((x>0) for x in w))%2 for w in vec]
        print(" 面の中心",v,"の隣",nb,"格子の中",sum(inside(q) for q in nb)); break
# 単位格子あたり：原子 8、結合 16
cnt=sum(1 for p in A if all(0<=x<12 for x in p)); bonds=sum(1 for p in A if all(0<=x<12 for x in p) for q in A if d2(p,q)==27)//2
print("単位格子あたり 原子",cnt,"結合",bonds)
