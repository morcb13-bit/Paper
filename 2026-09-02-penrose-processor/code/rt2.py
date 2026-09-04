# 検定RT2（広げた担体で梯子が続くか）
#
#   担体を半径260まで広げた。φ¹⁰ = 122.9919 の殻が入る。
#   φ⁴・φ⁶・φ⁸・φ¹⁰ の単数の星それぞれについて、72°回転が効く半径 R を測る。
#
#   OK   φ¹⁰ の5個が新しい R 類を作り、φ⁴→φ⁶→φ⁸→φ¹⁰ と R が上がり続ける
#   NG   φ¹⁰ が単数なのに R が上がらない（φ⁸ に担体と別の境界がある）
#   必ず落ちる設定  36°回転はどこでも 0
#   検算       図の中心の R が縁までの余裕と一致すること（担体が5回対称であること）
import json, math
from collections import Counter, defaultdict
import b13_chain_units_big as U

cells = {tuple(int(x) for x in k.split(",")): a
         for k, a in json.load(open("carrier_big.json"))["cells"].items()}
adj = defaultdict(set)
for q, a in cells.items():
    vs = [U.zadd(q, U.zt(a + 2 * i)) for i in range(5)]
    for i in range(5):
        u, w = vs[i], vs[(i + 1) % 5]
        adj[u].add(w); adj[w].add(u)
V10 = set(tuple(10 * x for x in p) for p in adj)
XY10 = {p: U.xy(p) for p in V10}
Rcar = max(math.hypot(*v) for v in XY10.values()) / 10.0
print("五角形 %d  頂点 %d  担体の半径 %.3f" % (len(cells), len(V10), Rcar))

CEN = []
for area, cyc in U.gaps(cells):
    if abs(area - 2.9389) < 0.01:
        c = (0, 0, 0, 0)
        for p in cyc: c = U.zadd(c, p)
        if any(len(adj[p]) < 2 for p in cyc): continue
        CEN.append(c)
print("五芒星 %d 個" % len(CEN))

def N(p, q): return p * p + p * q - q * q
phi = (1 + 5 ** .5) / 2

def closes(c, k=2, cap=None):
    z = U.zt(k); cx, cy = U.xy(c)
    lim = Rcar - math.hypot(cx, cy) / 10.0
    if cap: lim = min(lim, cap)
    g = defaultdict(list)
    for p in V10:
        d = math.hypot(XY10[p][0] - cx, XY10[p][1] - cy) / 10.0
        if d <= lim: g[round(d, 9)].append(p)
    R = 0.0
    for d in sorted(g):
        if all(U.zadd(c, U.zmul(U.zsub(p, c), z)) in V10 for p in g[d]): R = d
        else: return R, lim
    return R, lim

# 単数の星を拾う
units = defaultdict(list)
for c in CEN:
    p, q = U.norm2(c)
    if p % 100 or q % 100: continue
    p, q = p // 100, q // 100
    if abs(N(p, q)) == 1:
        units[round(math.hypot(*U.xy(c)) / 10.0, 6)].append(c)
print()
print("r² が単数の星")
for r in sorted(units):
    k = round(math.log(r) / math.log(phi))
    print("  r = %11.6f  = φ^%-2d   %2d 個" % (r, k, len(units[r])))

print()
print("検算：図の中心の R")
c0 = min(CEN, key=lambda c: math.hypot(*U.xy(c)))
R, lim = closes(c0)
print("  R = %.4f   縁までの余裕 = %.4f   一致 %s" % (R, lim, abs(R - lim) < 1e-9))

print()
print("必ず落ちる設定：36°回転")
print("  φ⁸ の星で R = %.4f" % closes(units[round(phi**8, 6)][0], k=1)[0])

print()
print("72°回転が効く半径")
for r in sorted(units):
    k = round(math.log(r) / math.log(phi))
    rs = []
    for c in units[r]:
        R, lim = closes(c)
        rs.append((R, lim))
    u = sorted(set(round(x[0], 4) for x in rs))
    print("  φ^%-2d (r=%10.6f)  R = %s   縁までの余裕 %.3f  個数 %d"
          % (k, r, u, rs[0][1], len(rs)))
