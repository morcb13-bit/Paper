# 正 n 単体を n+1 次元の整数座標で：頂点 v_i = (n+1)e_i - (1,…,1)、中心は原点
# 手どうしの組（内積, 長さ²）＝ (-(n+1), n(n+1)) → 約分して (-1, n)
from itertools import combinations
for n in range(1,6):
    d=n+1
    V=[[(d if k==i else 0)-1 for k in range(d)] for i in range(d)]
    s=[sum(c) for c in zip(*V)]
    prs={(sum(a*b for a,b in zip(u,v)), sum(a*a for a in u)) for u,v in combinations(V,2)}
    print(f"n={n} 頂点{d} 和{s} 組{prs}")
# 実測との照合（換算だけ浮動小数）
import math
for nm,deg,n in [("BeH2",180.0,1),("BH3",120.0,2),("CH4",109.4712,3),("H2O",104.4776,4)]:
    print(nm, deg, "模型", round(math.degrees(math.acos(-1/n)),4))
print("NH3 106.769  幾何平均 √(3·4):", round(math.degrees(math.acos(-1/math.sqrt(12))),4))
