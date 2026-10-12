# 検定2：手の向き（星型八面体の B 先端の空き席）の角を整数の組で書き、実測と照合する
# 判定は整数（内積, 長さ²）の交差掛け算だけ。浮動小数は外部の実測角を組に換算する一回だけ。
from itertools import combinations
import math, random

A = [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
B = [tuple(-x for x in v) for v in A]
SEATS = A + B                      # A から埋める順
dot = lambda u,v: sum(a*b for a,b in zip(u,v))

def hands(e):                      # L殻の電子数 e → 空き席（B側に残る）
    return SEATS[e:]

def pairs(hs):                     # 手どうしの (内積, 長さ²) の集合
    return sorted({(dot(u,v), dot(u,u)) for u,v in combinations(hs,2)})

# 外部の実測（[知識] NIST CCCBDB 気相）：この換算だけ浮動小数
MEAS = {"CH4":(4,109.47), "NH3":(5,106.67), "H2O":(6,104.48)}
def to_pair(deg, tol=0.1, qmax=12):  # 角 ±tol に入る最小分母の (p,q)：cos = p/q
    lo, hi = sorted([math.cos(math.radians(deg-tol)), math.cos(math.radians(deg+tol))])
    for q in range(1,qmax+1):
        for p in range(-q,q+1):
            if lo <= p/q <= hi: return (p,q)
    return None

same = lambda a,b: a[0]*b[1] == b[0]*a[1]      # p/q の一致を整数で

print("分子  手  模型(内積,長さ²)  実測の組  判定")
for name,(e,deg) in MEAS.items():
    hs = hands(e); mp = pairs(hs); m = to_pair(deg)
    ok = all(same(x,m) for x in mp)
    print(f"{name:4} {len(hs):2}  {mp}  {m}  {'OK' if ok else 'NG'}")

print("\n負の対照：8先端からばらばらに手を選んだときに出る組")
allp = sorted({(dot(u,v),dot(u,u)) for u,v in combinations(SEATS,2)})
print(" 出うる組", allp)
random.seed(1)
for k in (4,3,2):
    cnt = {}
    for _ in range(2000):
        key = tuple(pairs(random.sample(SEATS,k))); cnt[key]=cnt.get(key,0)+1
    print(f" 手{k}本:", sorted(cnt.items(), key=lambda x:-x[1])[:4])
