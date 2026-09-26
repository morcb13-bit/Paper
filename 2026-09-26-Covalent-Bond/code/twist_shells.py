# ひねって2頂点で接する組（中心間＝外接球の半径）の殻と、二つの中心を見込む角。sympy で厳密に。
import sympy as sp
from itertools import product, combinations
phi=(1+sp.sqrt(5))/2
V=set()
for s in product([1,-1],repeat=3): V.add(s)
for a in [phi-1,-(phi-1)]:
    for b in [phi,-phi]:
        base=[0,a,b]
        for k in range(3): V.add(tuple(base[(i-k)%3] for i in range(3)))
V=[sp.Matrix(v) for v in sorted(V,key=str)]
assert len(V)==20
cases={}
for i,j in combinations(range(20),2):
    g=sp.nsimplify(sp.simplify(V[i].dot(V[j])))
    if g in cases or sp.N(9+6*g)<=0: continue
    cases[g]=(i,j)
target76=sp.simplify(12*(sp.sqrt(5)-2)**2)          # cos=1/φ³ → 76.35°
target77=sp.Rational(12*49,1024)                     # cos=7/32 → 77.36°
for g,(i,j) in cases.items():
    v1,v2=V[i],V[j]
    al=sp.Rational(3,2)/(3+g); c=v1.cross(v2); be=sp.sqrt(sp.simplify((3-al**2*(6+2*g))/c.dot(c)))
    u=sp.simplify(al*(v1+v2)+be*c)
    assert sp.simplify(u.dot(u)-3)==0
    B=[sp.simplify(x-2*((x.dot(u)-sp.Rational(3,2))/3)*u) for x in V]
    pts=[]; 
    for x in V: pts.append(('A',x))
    for x in B:
        if not any(sp.simplify((x-y).dot(x-y))==0 for y in V): pts.append(('B',x))
    sh={}
    for who,x in pts:
        own=sp.Matrix([0,0,0]) if who=='A' else u
        oth=u if who=='A' else sp.Matrix([0,0,0])
        key=sp.simplify((x-oth).dot(x-oth))        # もう一方の中心までの距離²（自分の中心までは常に R²=3）
        sh.setdefault(key,[]).append(who)
    print(f"\n=== 接する2頂点の内積 g={g}（{sp.N(g,4)}）  点 {len(pts)} 個  殻 {len(sh)} 枚")
    keys=sorted(sh,key=lambda k:sp.N(k))
    for k in keys:
        cosv=sp.N(sp.sqrt(k)/(2*sp.sqrt(3)))   # 近い方は常に R：cos = d遠/(2R)
        ang=sp.N(sp.acos(sp.sqrt(k)/(2*sp.sqrt(3)))*180/sp.pi,6)
        e76=sp.simplify(k-target76)==0; e77=sp.simplify(k-target77)==0
        print(f"  相手の中心までの距離² {sp.N(k,8):>12}  点 {len(sh[k]):2d}（A{sh[k].count('A')}/B{sh[k].count('B')}）  見込む角 {ang}°  76.35°:{e76} 77.4°:{e77}")
    s=[sp.N(k) for k in keys]
    print("  向かい合う組の和:",[round(float(s[t]+s[-1-t]),9) for t in range(len(s)//2)])
