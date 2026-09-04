# 検定RT（回転の効く半径）
#
#   「5回対称に数を置いて72°回して重ねる」を盤面の操作として成立させるには、
#   その五芒星のまわりで 72°回転がタイリングを自分自身に写す必要がある。
#   どこまでの半径で写るかを、五芒星411個それぞれについて測る。
#
#   回転は Z[ζ₁₀] の厳密演算（ζ²倍 = 72°）。浮動小数は距離の並べ替えにしか使わない。
#
#   出すもの  R(i) = 中心 c(i) から見て、72°回転で頂点が頂点に写る最大半径
#   OK        R が φ³ より大きい星があり、R の分布が離散的（段をなす）
#   NG        中心以外は全部 φ³ 止まり（回転中心は図の中心だけ）
#   必ず落ちる設定  36°回転（ζ倍）ではどの星も閉じないこと
#   負の対照  五芒星でない点（円環中心）を中心にすると閉じないこと
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
# 10倍格子（五芒星の中心が整数になる）
V10 = set(tuple(10 * x for x in p) for p in adj)
XY10 = {p: U.xy(p) for p in V10}
Rcar = max(math.hypot(*v) for v in XY10.values()) / 10.0

CEN = []; 
for area, cyc in U.gaps(cells):
    if abs(area - 2.9389) < 0.01:
        c = (0, 0, 0, 0)
        for p in cyc: c = U.zadd(c, p)
        if any(len(adj[p]) < 2 for p in cyc): continue
        CEN.append(c)
RING = [tuple(10 * x for x in r) for r in cells]   # 円環中心（負の対照用）
print("頂点 %d  五芒星 %d  担体の半径 %.3f" % (len(V10), len(CEN), Rcar))

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


# --- 本番：72°（k=2） ---
res = []
for i, c in enumerate(CEN):
    R, lim = closes(c, 2)
    res.append((R, lim, i))
Rs = [r for r, l, i in res]
print()
print("72°回転で閉じる半径 R の分布（五芒星411個）")
h = Counter(round(r, 4) for r in Rs)
for v, cnt in sorted(h.items()):
    print("   R = %8.4f   %3d 個" % (v, cnt))
print("   最大 %.4f（φ³ = 4.2361、φ⁵ = 11.0902、φ⁷ = 29.0344）" % max(Rs))

# 縁で切れているものを除いた見方
cut = sum(1 for R, lim, i in res if abs(R - lim) < 1e-9)
print("   縁まで一度も破れなかった星 %d 個（R が余裕 lim と一致）" % cut)

# --- 必ず落ちる設定：36° ---
bad = [closes(c, 1)[0] for c in CEN[:60]]
print()
print("必ず落ちる設定 36°回転（先頭60個）  R の最大 %.4f" % max(bad))

# --- 負の対照：円環中心 ---
ctrl = [closes(c, 2)[0] for c in RING[:60]]
print("負の対照 円環中心を中心に72°（先頭60個）  R の最大 %.4f" % max(ctrl))

# --- R の大きい星はどこにいるか ---
print()
print("R の大きい順に10個（中心からの距離つき）")
res.sort(reverse=True)
for R, lim, i in res[:10]:
    d = math.hypot(*U.xy(CEN[i])) / 10.0
    print("   R %8.4f  （縁までの余裕 %7.4f）  中心からの距離 %8.4f" % (R, lim, d))
