# 検定WD（窓の中で決まるか）
#
#   検定RT で、五芒星411個は「72°回転が破れる半径 R」で5つの類に分かれた。
#   R が星の窓の中の位置（垂直空間の像）で決まるなら、回転はカットアンドプロジェクトの
#   言葉に載る。決まらないなら、R は窓の外の量ということになる。
#
#   窓の座標  σ: ζ → ζ³ による共役像。Z[ζ₁₀] の厳密演算で出す
#
#   OK   R の類ごとに窓の中の領域が分かれる（別の類の点が混ざらない）
#   NG   混ざる（R は窓では決まらない）
#   必ず落ちる設定  窓の座標を無作為に振り直すと分離が消えること
#   併せて  5個の特別な星が窓のどこにいるか
import json, math, random
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
V10 = set(tuple(10 * x for x in p) for p in adj)
XY10 = {p: U.xy(p) for p in V10}
Rcar = max(math.hypot(*v) for v in XY10.values()) / 10.0
CEN = []
for area, cyc in U.gaps(cells):
    if abs(area - 2.9389) < 0.01:
        c = (0, 0, 0, 0)
        for p in cyc: c = U.zadd(c, p)
        if any(len(adj[p]) < 2 for p in cyc): continue
        CEN.append(c)
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


def conj(v):
    """σ: ζ → ζ³（窓の座標）"""
    out = (0, 0, 0, 0)
    for i, coef in enumerate(v):
        if coef: out = U.zadd(out, tuple(coef * x for x in U.zt((3 * i) % 10)))
    return out

R = [round(closes(c)[0], 4) for c in CEN]
W = [tuple(x / 10.0 for x in U.xy(conj(c))) for c in CEN]
print("五芒星 %d 個   R の類 %s" % (len(CEN), dict(sorted(Counter(R).items()))))
print()

# 類ごとの窓の中の広がり
print("窓の中での各類の位置")
for r in sorted(set(R)):
    pts = [W[i] for i in range(len(CEN)) if R[i] == r]
    rad = [math.hypot(*p) for p in pts]
    print("  R = %8.4f  %3d 個   窓半径 %.4f 〜 %.4f  （平均 %.4f）"
          % (r, len(pts), min(rad), max(rad), sum(rad) / len(rad)))
print()

# 分離しているか：別の類で窓座標が近い対があるか
def sep(Rlist):
    worst = None
    for i in range(len(CEN)):
        for j in range(i + 1, len(CEN)):
            if Rlist[i] == Rlist[j]: continue
            d = math.hypot(W[i][0] - W[j][0], W[i][1] - W[j][1])
            if worst is None or d < worst[0]: worst = (d, i, j)
    return worst
d, i, j = sep(R)
print("別の類どうしで窓の中が最も近い対： 距離 %.6f （R %.4f と R %.4f）" % (d, R[i], R[j]))
# 同じ類の中の最小距離と比べる
same = None
for a in range(len(CEN)):
    for b in range(a + 1, len(CEN)):
        if R[a] != R[b]: continue
        dd = math.hypot(W[a][0] - W[b][0], W[a][1] - W[b][1])
        if same is None or dd < same: same = dd
print("同じ類どうしの窓の中の最小距離： %.6f" % same)
print("→ 別の類が近づく距離が、同じ類の最小距離より大きければ、窓で分かれている")
print()

# 必ず落ちる設定：窓の座標を振り直す
rng = random.Random(7)
Wsave = W[:]
idx = list(range(len(CEN))); rng.shuffle(idx)
W = [Wsave[k] for k in idx]
d2, i2, j2 = sep(R)
print("必ず落ちる設定（窓の座標を振り直し）  別の類の最小距離 %.6f" % d2)
W = Wsave
print()

# 5個の特別な星は窓のどこにいるか
print("R が大きい上位6個の窓の座標")
order = sorted(range(len(CEN)), key=lambda i: -R[i])[:6]
for i in order:
    print("  R %9.4f   窓 (%8.4f, %8.4f)  |窓| %.4f   実空間の距離 %8.4f"
          % (R[i], W[i][0], W[i][1], math.hypot(*W[i]), math.hypot(*U.xy(CEN[i])) / 10.0))
