#  検定WK4  異なる列は異なる軌跡になるか
#
#      監督の指定（2026-09-02）：
#          1,0,0,-1,0,0,1,1,-1 といった異なる情報を流していくと、決して同じ軌跡にはならない。
#
#      前の検定WK3 は記号を {連, 飛}（辺の種類）に取ったため走らなかった。誤り。
#      入力は辺の種類ではなく**頂点での選び方**であり、次数がそのまま選択肢の数になる。
#          次数2 … 出口1本 → 0 のみ
#          次数3 … 出口2本 → −1 と +1
#          次数4 … 出口3本 → −1, 0, +1 の三つとも使える
#      よって4または詰まる場所ではなく、桁が全部立つ場所である。
#
#      測る量：長さ k の全ての列 3^k について
#          (1) 最後まで流れるか（その頂点にその桁が無ければ止まる）
#          (2) 終わりの状態（前の頂点, 今の頂点）が列ごとに異なるか
#          (3) 軌跡（通った頂点の並び）が列ごとに異なるか
#      OK なら：流れた列の終わりの状態が列を一意に決める（単射）＝入力が状態を決める
#      NG なら：異なる列が同じ状態に落ちる（衝突）＝情報が消える
#      必ず落ちる設定：桁を全部 0 にした列と、全部 +1 の列が同じ軌跡になったら物差しが壊れている
#      対照：同じ長さの列を無作為に選んでも衝突が出ないこと（衝突が構造由来か偶然かを分ける）
#
#      桁の割り当て：入ってきた向きを基準に、出口を折れ角の順に並べる。
#          右から順に +1, 0, −1（出口が2本なら +1 と −1、1本なら 0）

import json, math, itertools
from collections import defaultdict, Counter
import b13_chain_units as U

d = json.load(open("rings_integer.json"))
cells = {}
for k, v in d["cells"].items():
    cells[tuple(json.loads(k)) if k.startswith("[") else tuple(eval(k))] = v
adj = defaultdict(set)
for q, a in cells.items():
    vs = [U.zadd(q, U.zt(a + 2 * i)) for i in range(5)]
    for i in range(5):
        u, w = vs[i], vs[(i + 1) % 5]
        adj[u].add(w); adj[w].add(u)
XY = {p: U.xy(p) for p in adj}
print("頂点 %d／次数 %s" % (len(adj), dict(sorted(Counter(len(v) for v in adj.values()).items()))))

def options(prev, cur):
    """入ってきた向きを基準に、出口を右から順に並べて桁を割り当てる"""
    outs = [w for w in adj[cur] if w != prev]
    if not outs:
        return {}
    back = math.atan2(XY[prev][1] - XY[cur][1], XY[prev][0] - XY[cur][0])
    rel = sorted(((math.atan2(XY[w][1] - XY[cur][1], XY[w][0] - XY[cur][0]) - back) % (2 * math.pi), w)
                 for w in outs)
    if len(rel) == 1:
        return {0: rel[0][1]}
    if len(rel) == 2:
        return {1: rel[0][1], -1: rel[-1][1]}
    return {1: rel[0][1], 0: rel[len(rel) // 2][1], -1: rel[-1][1]}

def flow(word, start):
    prev, cur = start
    trail = [prev, cur]
    for s in word:
        o = options(prev, cur)
        if s not in o:
            return None, None
        prev, cur = cur, o[s]
        trail.append(cur)
    return (prev, cur), tuple(trail)

# 出発：中心の五芒星の角から外向きに一歩
faces = U.gaps(cells)
star = None
for area, cyc in faces:
    if abs(area - 2.9389) < 0.01:
        P = [U.xy(p) for p in cyc]
        g = (sum(a for a, _ in P) / len(P), sum(b for _, b in P) / len(P))
        if math.hypot(*g) < 0.01:
            star = cyc; break
v0 = star[0]
v1 = max(adj[v0], key=lambda w: math.hypot(*XY[w]))
start = (v0, v1)
print("出発（中心の五芒星の角から外向き）\n")

print("■ 列の長さごとに")
print("  長さ   列の数   流れた   終わりの状態   軌跡   衝突なし")
for k in range(1, 11):
    ends = defaultdict(list); trails = set()
    n = 0
    for w in itertools.product((-1, 0, 1), repeat=k):
        e, t = flow(w, start)
        if e is None: continue
        n += 1; ends[e].append(w); trails.add(t)
    uniq = sum(1 for v in ends.values() if len(v) == 1)
    print("   %2d   %6d   %6d   %10d   %6d   %s"
          % (k, 3 ** k, n, len(ends), len(trails),
             ("%.3f" % (uniq / n)) if n else "-"))

print("\n■ 監督の列 1,0,0,-1,0,0,1,1,-1")
w = (1, 0, 0, -1, 0, 0, 1, 1, -1)
e, t = flow(w, start)
print("  %s" % ("流れない" if e is None else "流れた。頂点 %d 個の軌跡" % len(t)))
if e is not None:
    print("  終わりの位置 (%.3f, %.3f)" % XY[e[1]])

print("\n■ 必ず落ちる設定：全部0 と 全部+1 は別の軌跡か")
for k in (5, 8):
    a = flow((0,) * k, start)[1]; b = flow((1,) * k, start)[1]
    print("  長さ%d：全部0 %s／全部+1 %s／%s"
          % (k, "流れない" if a is None else "流れた",
             "流れない" if b is None else "流れた",
             "別の軌跡" if (a is None or b is None or a != b) else "同じ軌跡（物差しが壊れている）"))

print("\n■ 衝突の中身（長さ10で終わりの状態が同じになった列）")
ends = defaultdict(list)
for w in itertools.product((-1, 0, 1), repeat=10):
    e, t = flow(w, start)
    if e is not None: ends[e].append(w)
col = [v for v in ends.values() if len(v) > 1]
print("  衝突した状態 %d 個／最大 %d 本が同じ状態へ"
      % (len(col), max((len(v) for v in col), default=0)))
if col:
    v = max(col, key=len)
    for w in v[:3]:
        print("    %s" % (",".join("%+d" % x for x in w)))
