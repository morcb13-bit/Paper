# 正四面体の一本道 と RT の管 の関係（厳密計算：sympy, Z[φ]）
# 事前の基準
#  Q1 管の一層の一歩（中心の変位 t）の横²:縦² が 1:4 か（一本道の中心の一歩と同じ円錐か）
#     対照：輪の中の面接触の一歩 30 本、管で試した他の一歩 → 1:4 以外が出れば検査は NG を返せる
#  Q2 一本道の支点辺と軸の cos²=1/5 が、RT の五回軸どうしの cos² と一致するか
#  Q3 円錐の上の一歩の回転角：一本道 cosΔ=-1/4、RT 72°
import sympy as sp, itertools
p=(1+sp.sqrt(5))/2
G=[sp.Matrix(v) for v in [(1,p,0),(-1,p,0),(0,1,p),(0,-1,p),(p,0,1),(p,0,-1)]]
S2=1+p**2
def phys(n): return sum((G[k]*n[k] for k in range(6)),sp.zeros(3,1))  # 長さ×√S2
ax=G[2]
def ratio(n):
    v=phys(n); h2=sp.nsimplify(sp.expand((v.dot(ax))**2/(ax.dot(ax))))
    l2=sp.expand(v.dot(v)-h2)
    return sp.radsimp(sp.simplify(l2/h2)) if h2!=0 else None, sp.simplify(v.dot(v)/S2)
# 一本道の中心の一歩：前進²12/5、中心間²3 → 横²=3/5
print("一本道 中心の一歩 横²:縦² =",sp.Rational(3,5)/sp.Rational(12,5))
r,L=ratio((2,2,2,0,1,-1))
print("Q1 管の一層 t=(2,2,2,0,1,-1): 横²/縦² =",r,"  長さ² =",sp.nsimplify(sp.radsimp(L)))
# 対照：面接触の一歩（±1が4つ）で面法線が…→ 30本すべての比を集計
steps=set()
for idx in itertools.combinations(range(6),4):
    for sg in itertools.product((1,-1),repeat=4):
        n=[0]*6
        for i,s in zip(idx,sg): n[i]=s
        steps.add(tuple(n))
from collections import Counter
C=Counter()
for n in steps:
    v=phys(n); 
    if v.dot(v)==0: continue
    rr,_=ratio(n); C[str(rr)]+=1
print("対照（±1が4つの一歩 全",len(steps),"本）横²/縦² の分布:")
for k,v in C.most_common(): print("   ",k,v)
# Q2 五回軸どうし
c2=sp.simplify((G[0].dot(G[2]))**2/(G[0].dot(G[0])*G[2].dot(G[2])))
print("Q2 隣り合う五回軸の cos² =",sp.radsimp(c2))
# Q3
print("Q3 円錐 cos²=1/5 の上で隣が直交する回転: cosΔ =",sp.solve(sp.Rational(1,5)+sp.Rational(4,5)*sp.Symbol('c'),'c'))
print("    RT 五回軸の隣: 1/5+4/5cos72 =",sp.radsimp(sp.Rational(1,5)+sp.Rational(4,5)*sp.cos(sp.pi*2/5)), " = 1/√5 ?",sp.simplify(sp.Rational(1,5)+sp.Rational(4,5)*sp.cos(sp.pi*2/5)-1/sp.sqrt(5))==0)
