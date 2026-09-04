# 検定UN（単数の星と回転中心の照合）
#
#   検定WD で、72°回転が遠くまで効く星は中心1個と φ⁸ の5個だった。
#   その5個は r² = φ¹⁶ = 610 + 987φ で N = 1、つまり単数。
#   では 411個のうち r² が単数になる星は、その5個だけか。
#
#   OK   単数の星＝回転中心の5個（一語で繋がる）
#   NG   別の集合（単数と回転は無関係）
#   必ず落ちる設定  N(r²) の計算が φ² 倍で 1 を保つこと（単数の定義の検算）
import json, math
from collections import Counter, defaultdict
import b13_chain_units as U

cells = {tuple(int(x) for x in k.split(",")): a
         for k, a in json.load(open("carrier_1245.json"))["cells"].items()}
adj = defaultdict(set)
for q, a in cells.items():
    vs = [U.zadd(q, U.zt(a + 2 * i)) for i in range(5)]
    for i in range(5):
        u, w = vs[i], vs[(i + 1) % 5]
        adj[u].add(w); adj[w].add(u)
CEN = []
for area, cyc in U.gaps(cells):
    if abs(area - 2.9389) < 0.01:
        c = (0, 0, 0, 0)
        for p in cyc: c = U.zadd(c, p)
        if any(len(adj[p]) < 2 for p in cyc): continue
        CEN.append(c)
phi = (1 + 5 ** .5) / 2

def N(p, q): return p * p + p * q - q * q      # N(p+qφ)

print("必ず落ちる設定：φ の冪のノルムが (−1)^k であること")
x = (1, 0)
for k in range(1, 9):
    x = (x[1], x[0] + x[1])                    # φ^k = F_{k-1} + F_k φ
    if k in (2, 4, 8, 16): print("   φ^%-2d = %4d + %4dφ   N = %3d" % (k, x[0], x[1], N(*x)))
print()

rows = []
for c in CEN:
    p, q = U.norm2(c)                          # 10倍座標の二乗長 = 100·r²
    assert p % 100 == 0 and q % 100 == 0
    p, q = p // 100, q // 100
    rows.append((math.hypot(*U.xy(c)) / 10.0, p, q, N(p, q), c))

print("五芒星 411 個の r² と N(r²)")
agg = defaultdict(int)
for r, p, q, nn, c in rows: agg[(round(r, 6), p, q, nn)] += 1
for (r, p, q, nn), m in sorted(agg.items()):
    mark = "  ← 単数" if abs(nn) == 1 else ""
    print("   r %9.6f   r² = %5d + %5dφ   N = %8d   %3d 個%s" % (r, p, q, nn, m, mark))

units = [(r, c) for r, p, q, nn, c in rows if abs(nn) == 1]
print()
print("r² が単数の星 %d 個  距離 %s" % (len(units), sorted(set(round(r, 6) for r, c in units))))

# 回転中心（検定WD の5個）と照合
V10 = set(tuple(10 * x for x in p) for p in adj)
XY10 = {p: U.xy(p) for p in V10}
Rcar = max(math.hypot(*v) for v in XY10.values()) / 10.0
PTS = sorted(V10, key=lambda p: math.hypot(*XY10[p]))
def closes(c, k=2):
    """同じ距離の頂点はまとめて判定する。個別に見ると同距離の並び順で R が変わる"""
    z = U.zt(k); cx, cy = U.xy(c); lim = Rcar - math.hypot(cx, cy) / 10.0
    g = defaultdict(list)
    for p in V10:
        d = math.hypot(XY10[p][0] - cx, XY10[p][1] - cy) / 10.0
        if d <= lim: g[round(d, 9)].append(p)
    R = 0.0
    for d in sorted(g):
        if all(U.zadd(c, U.zmul(U.zsub(p, c), z)) in V10 for p in g[d]): R = d
        else: return R, lim
    return R, lim

rot = set(c for c in CEN if closes(c)[0] > 20)
uni = set(c for r, c in units)
print()
print("回転中心（R>20）%d 個  /  単数の星 %d 個  /  共通 %d 個"
      % (len(rot), len(uni), len(rot & uni)))
print("回転中心だが単数でない %d 個   単数だが回転中心でない %d 個"
      % (len(rot - uni), len(uni - rot)))
