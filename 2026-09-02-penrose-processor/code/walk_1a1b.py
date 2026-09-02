#  検定WK1  1a/1b の歩き手
#
#      監督の指定（2026-09-02）：
#          五芒星のところで 1a・1b と出てきて順繰りに埋めていく。
#          直進しか選べないときは 1a も 1b も同じ。3またのときは 1a は右、1b は左に進む。
#
#      測る量：中心の五芒星の各角から 1a・1b を出し、
#              (1) 4またに当たるまで／図の外へ出るまでの歩数
#              (2) 通った次数の列
#              (3) 二本が閉じるか（元の場所に戻るか）、鏡像になっているか
#      OK なら：二本が規則正しく分かれ、次数の列に周期か規則が出る
#      NG なら：すぐ4またに当たる／二本が同じ道になる／規則が読めない
#      未定義：4また（内部1330頂点のうち304個）。規則が何も言っていないので、
#              当たったらそこで止めて記録する。ここを決めるのは監督の側
#      必ず落ちる設定：右と左を入れ替えると、二本の役が入れ替わるだけであること

import json, math
from collections import Counter, defaultdict
import b13_chain_units as U

d = json.load(open("rings_integer.json"))
cells = {}
for k, v in d["cells"].items():
    q = tuple(json.loads(k)) if k.startswith("[") else tuple(eval(k))
    cells[q] = v

# ── 五角形の辺のグラフ ────────────────────────────────
adj = defaultdict(set)
for q, a in cells.items():
    vs = [U.zadd(q, U.zt(a + 2 * i)) for i in range(5)]
    for i in range(5):
        u, w = vs[i], vs[(i + 1) % 5]
        adj[u].add(w); adj[w].add(u)
XY = {p: U.xy(p) for p in adj}
deg = Counter(len(s) for s in adj.values())
print("頂点 %d 個／次数 %s" % (len(adj), dict(sorted(deg.items()))))

# ── 中心の五芒星の角 ──────────────────────────────────
faces = U.gaps(cells)
star = None
for area, cyc in faces:
    if abs(area - 2.9389) < 0.01:
        P = [U.xy(p) for p in cyc]
        g = (sum(a for a, _ in P) / len(P), sum(b for _, b in P) / len(P))
        if math.hypot(*g) < 0.01:
            star = cyc; break
print("中心の五芒星の角 %d 個\n" % (len(star) if star else 0))

# ── 歩き手 ────────────────────────────────────────────
def ang(a, b):
    return math.atan2(XY[b][1] - XY[a][1], XY[b][0] - XY[a][0])

def step(prev, cur, side):
    """cur に prev から入ったときの次の頂点。side=+1 が右、-1 が左"""
    outs = [w for w in adj[cur] if w != prev]
    if not outs:
        return None, "行き止まり"
    if len(outs) == 1:
        return outs[0], "直進"
    if len(outs) == 2:
        back = ang(cur, prev)
        rel = []
        for w in outs:
            t = (ang(cur, w) - back) % (2 * math.pi)      # 後ろ向きから測った角
            rel.append((t, w))
        rel.sort()
        # 後ろから時計回りに最初＝右、最後＝左（画面 y は下向きなので入れ替え）
        return (rel[-1][1] if side > 0 else rel[0][1]), "3また"
    return None, "4また"

def walk(start_prev, start, side, limit=400):
    prev, cur = start_prev, start
    path = [start_prev, start]
    kinds = []
    for _ in range(limit):
        nxt, kind = step(prev, cur, side)
        kinds.append(kind)
        if nxt is None:
            return path, kinds, kind
        prev, cur = cur, nxt
        path.append(cur)
        if cur == start and prev == start_prev:
            return path, kinds, "閉じた"
    return path, kinds, "打ち切り"

print("■ 中心の五芒星の10角から出す")
res = []
for i, v in enumerate(star):
    # 五芒星の角から外へ：角に接する辺のうち、原点から遠ざかる向きを最初の一歩にする
    outs = sorted(adj[v], key=lambda w: -math.hypot(*XY[w]))
    if not outs:
        continue
    first = outs[0]
    for side, nm in ((+1, "1a 右"), (-1, "1b 左")):
        path, kinds, end = walk(v, first, side)
        res.append((i, nm, len(path), Counter(kinds), end))

for i, nm, n, c, end in res[:8]:
    print("  角%2d %s … 歩数 %3d ／ %s ／ 終わり方 %s"
          % (i, nm, n, dict(c), end))

print("\n■ まとめ")
ends = Counter(e for _, _, _, _, e in res)
lens = [n for _, _, n, _, _ in res]
print("  終わり方 %s" % dict(ends))
print("  歩数 最小%d 最大%d 中央%d" % (min(lens), max(lens), sorted(lens)[len(lens)//2]))

# 1a と 1b が同じ道になっていないか
same = 0
for i in range(0, len(res), 2):
    if res[i][2] == res[i + 1][2] and res[i][3] == res[i + 1][3]:
        same += 1
print("  1a と 1b が同じ形になった角 %d / %d" % (same, len(res) // 2))
