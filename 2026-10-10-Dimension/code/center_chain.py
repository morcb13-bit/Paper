# 検定 F1：中心を通る中の立方体を同じ比でたどったとき、どの次元で 1（中心を通る線分一本）に着くか
# 基準（走らせる前に決める）
#  - 一歩：n次元立方体 {0,1}^n の中の k立方体のうち、中心が n次元立方体の中心と一致するもの
#    （頂点の和×2 ＝ 2^k ×(1,…,1) を整数で判定）。辺の二乗 m を毎段同じにする（同じ比＝フラクタル）
#  - 着地：k＝1（中心を通る線分）。その線分の中点が中心の点
#  - 予想：m＝2 なら n＝2^j だけが 1 に着く。外れる次元が一つでもあれば NG
#  - 対照：m＝3 では n＝3^j だけが着く、と別の答えが出ること（比を変えると答えが変わる＝NGを返せる）
from itertools import product

def centered_inner(n, m):
    # 中心を通る中の立方体が存在するか（総当たりの代わりに一つ作って判定。n<=12）
    if n % m:
        return None
    k = n // m
    # ブロック i は座標 [i*m,(i+1)*m)。原点 (0..0)、辺のベクトルはブロック全部を 1 にする
    verts = []
    for s in product((0, 1), repeat=k):
        v = []
        for b in s:
            v += [b] * m
        verts.append(v)
    tot = [sum(v[i] for v in verts) for i in range(n)]
    ok = all(2 * t == len(verts) for t in tot)  # 頂点の和×2 ＝ 2^k
    return k if ok else None

def chain(n, m):
    path = [n]
    while n > 1:
        k = centered_inner(n, m)
        if k is None:
            return path, False
        path.append(k)
        n = k
    return path, True

for m in (2, 3):
    print(f"比（辺の二乗）m={m}")
    for n in range(1, 17):
        p, ok = chain(n, m)
        print(f"  {n:2d}: {' → '.join(map(str, p))}  {'1 に着く' if ok and p[-1]==1 else '途中で止まる'}")
