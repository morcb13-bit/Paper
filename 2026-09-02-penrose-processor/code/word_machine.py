#  検定WK3  語 → 位置。1a/1b の並びは計算機の入力になるか
#
#      監督の指定（2026-09-02）：1a,1a,1b,1b と変えるだけで計算機になるのでは。
#
#      前の測定で役が分かれた：
#          1a = 連（巡る）。108°で折り続けると5環で閉じる。輪21個・重心は五芒星と一致
#          1b = 飛（移る）。10本とも別の輪をつなぐ。輪の中には0本
#      よって語は {a,b} の列で、a＝108°の連を一歩、b＝飛を一歩。
#
#      計算機の骨格は「入力が状態を決めること」。それを測る。
#      測る量：長さ k の全ての語について
#          (1) 最後まで歩けるか
#          (2) 終わりの状態（前の環, 今の環）が語ごとに異なるか
#      OK なら：歩ける語が十分あり、終わりの状態が語をほぼ区別する
#      NG なら：ほとんどの語が歩けない／終わりが衝突して語が復元できない
#      必ず落ちる設定：語を逆順にすると別の終わりになること（順序が効いていること）
#      対照：a を「108°」ではなく「任意の連」にすると分岐して決まらなくなること

import json, math, itertools
from collections import defaultdict, Counter
import b13_chain_units as U

d = json.load(open("rings_integer.json"))
R = [tuple(r) for r in d["rings"]]; XY = {r: U.xy(r) for r in R}
NC = U.norm2(U.CONT[0]); NS = U.norm2(U.SKIP[0])
adj = defaultdict(dict)
for i, a in enumerate(R):
    for b in R[i + 1:]:
        n = U.norm2(U.zsub(a, b))
        if n == NC: adj[a][b] = "a"; adj[b][a] = "a"
        elif n == NS: adj[a][b] = "b"; adj[b][a] = "b"

def turn(p, c, n):
    ax, ay = XY[c][0] - XY[p][0], XY[c][1] - XY[p][1]
    bx, by = XY[n][0] - XY[c][0], XY[n][1] - XY[c][1]
    cr = ax * by - ay * bx; dt = ax * bx + ay * by
    return 180.0 - abs(math.degrees(math.atan2(cr, dt))), (1 if cr > 0 else -1)

SIDE = +1
def step(p, c, sym, free_a=False):
    """一歩。a=108°の連（free_a なら任意の連）、b=飛"""
    outs = [w for w in adj[c] if w != p and adj[c][w] == sym]
    if sym == "b":
        return outs[0] if len(outs) == 1 else None
    if free_a:
        return outs[0] if len(outs) == 1 else (None if len(outs) != 1 else outs[0])
    cand = [w for w in outs if abs(turn(p, c, w)[0] - 108) < 1e-6 and turn(p, c, w)[1] == SIDE]
    return cand[0] if len(cand) == 1 else None

def run(word, start, free_a=False):
    p, c = start
    for s in word:
        n = step(p, c, s, free_a)
        if n is None:
            return None
        p, c = c, n
    return (p, c)

# 出発は中心の五芒星の輪の上に固定する
inner = sorted(R, key=lambda r: math.hypot(*XY[r]))[:5]
start = None
for a in inner:
    for b in adj[a]:
        if adj[a][b] == "a":
            start = (a, b); break
    if start: break
print("出発 (前, 今) = %s\n" % (start,))

print("■ 語の長さごとに（出発は固定）")
print("  長さ  語の数  歩けた  終わりの状態の種類  衝突なしの割合")
for k in range(1, 9):
    ws = ["".join(w) for w in itertools.product("ab", repeat=k)]
    ends = {}
    for w in ws:
        e = run(w, start)
        if e is not None:
            ends.setdefault(e, []).append(w)
    walked = sum(len(v) for v in ends.values())
    uniq = sum(1 for v in ends.values() if len(v) == 1)
    print("   %2d   %5d   %5d      %5d          %s"
          % (k, len(ws), walked, len(ends),
             ("%.2f" % (uniq / walked)) if walked else "-"))

print("\n■ 監督の語 aabb を繰り返す")
for rep in range(1, 7):
    w = "aabb" * rep
    e = run(w, start)
    print("  %-24s → %s" % (w, "歩けない" if e is None else "環 %s" % (e[1],)))

print("\n■ 歩けた語だけを見る（長さ6）")
ws = ["".join(w) for w in itertools.product("ab", repeat=6)]
ok = [(w, run(w, start)) for w in ws]
ok = [(w, e) for w, e in ok if e is not None]
print("  歩けた %d / 64" % len(ok))
print("  そのうち b の本数の分布 %s" % dict(Counter(w.count("b") for w, _ in ok)))
print("  b が2つ続く語で歩けたもの %d 本" % sum(1 for w, _ in ok if "bb" in w))

print("\n■ 必ず落ちる設定：語を逆順にすると別の終わりか")
same = diff = dead = 0
for w, e in ok:
    e2 = run(w[::-1], start)
    if e2 is None: dead += 1
    elif e2 == e:  same += 1
    else:          diff += 1
print("  別の終わり %d ／ 同じ終わり %d ／ 逆順は歩けない %d" % (diff, same, dead))

print("\n■ 対照：a を任意の連にする（108°を外す）")
for k in (4, 6):
    ws = ["".join(w) for w in itertools.product("ab", repeat=k)]
    n = sum(1 for w in ws if run(w, start, free_a=True) is not None)
    print("  長さ%d：決まった一歩になった語 %d / %d" % (k, n, len(ws)))
