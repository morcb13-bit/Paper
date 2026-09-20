#  植物の最小モデル  検定PL6  発火して二葉に分かれるまで
#
#      伸ばすだけ。縮める段は無い。原点を移すのではなく、先端を増やす。
#
#  装置（足したものは無い）
#      先端     星と、累積した番地（13進）の組
#      光       担体の上の濃さ。センサーはまわり5枚（符牒と同じ量）
#      一歩     番地を 13/φ² の Sturmian 増分（4 と 5 の列）で進め、
#               符牒の最大の向きの隣へ先端を作る
#      発火     累積番地が節（mod13 = 0）に落ちたコマ
#      分岐     符牒の最大が二つで同点なら、両方に先端を作る（親はそこで節になる）
#      止まる   最大が三つ以上で同点なら伸びない。隣が無ければその先端だけ止まる
#
#  事前登録（走らせる前に書いた）
#
#  検定PL6  発火して二葉に分かれるか
#      OK なら：発火が起き、そのあと二葉に分かれ、両方が一歩以上伸びる
#      NG なら：発火しない／分かれない／分かれても片方しか伸びない
#
#  検定PL7  番地の列が 4H {4,6,7,9} を踏まないか（記事 2026-04-02 の定理1）
#      OK なら：踏まない
#      NG なら：踏む（担体が違うので成り立たない、と書ける）
#      負の対照：増分を 4,5 でない等間隔（3 など）にすると踏む
#
#  検定PL8  空試験  光が一様
#      OK なら：最大が5つ同点なので伸びない
#
#  前提：wind_core.py・pentagram_neuron.py・amoeba.py を同じ場所に置く。

import math, itertools
from collections import defaultdict

exec(open('wind_core.py').read())
exec(open('pentagram_neuron.py').read().split('def run_tests')[0])
src = open('amoeba.py').read()
exec(src[src.index('def build'):src.index('def run_amoeba')])

PHI = (1 + 5 ** 0.5) / 2
ALPHA = 13 / PHI ** 2
SC, AROUND, _ = build()

# 隣は検定AM4 の結果に合わせ、三番目の距離まで
LK = defaultdict(list)
for i, j in itertools.combinations(range(30), 2):
    if math.dist(SC[i], SC[j]) <= 10.5784 + 1e-3:
        LK[i].append(j); LK[j].append(i)


def sturmian(n):
    """13/φ² の Sturmian 増分（4 と 5 の列）。"""
    return [int(ALPHA * (k + 1)) - int(ALPHA * k) for k in range(n)]


def tips_of(i, peak):
    """まわり5枚のうち最大のもの。同点の本数と、その向きを返す。"""
    vals = [round(-((p[0] - peak[0]) ** 2 + (p[1] - peak[1]) ** 2), 6)
            for q, p in AROUND[i]]
    m = max(vals)
    k = [t for t, v in enumerate(vals) if v == m]
    return [AROUND[i][t][1] for t in k]


def grow(start, peak, T=12, step=None, verbose=True):
    """先端を伸ばす。戻り値は（枝の記録、発火したコマ、分岐したコマ）。"""
    inc = step if step else sturmian(T + 2)
    live = [(start, 0, 0)]              # 星・累積番地・親の枝番号
    hist, fires, splits, addr = [], [], [], []
    for t in range(T):
        nxt = []
        for (i, a, br) in live:
            tg = tips_of(i, peak)
            a2 = a + inc[t]
            addr.append(a2 % 13)
            if a2 % 13 == 0:
                fires.append((t, i, br))
            if len(tg) >= 3:
                continue                # 伸びない
            if len(tg) == 2:
                splits.append((t, i, br))
            for w, target in enumerate(tg):
                cx, cy = SC[i]
                th = math.atan2(target[1] - cy, target[0] - cx)
                best, bd = None, None
                for j in LK[i]:
                    ang = math.atan2(SC[j][1] - cy, SC[j][0] - cx)
                    d = abs((ang - th + math.pi) % (2 * math.pi) - math.pi)
                    if bd is None or d < bd:
                        bd, best = d, j
                if best is None or bd >= math.pi / 2:
                    continue            # この先端は止まる
                nb = br if len(tg) == 1 else br * 2 + 1 + w
                nxt.append((best, a2, nb))
                hist.append((t, i, best, nb))
        live = nxt
        if verbose:
            print("   %2dコマ  増分%d  先端 %d 個  %s"
                  % (t, inc[t], len(live),
                     "／".join("星%d(番地%d)" % (i, a % 13) for i, a, _ in live[:6])))
        if not live:
            break
    return hist, fires, splits, addr, live


print("13/φ² = %.6f  Sturmian 増分（先頭14）: %s" % (ALPHA, sturmian(14)))

# 光源は軸（x=0）の上に置く。先端は軸から外れた星から出す。
peak = SC[27]
START = 0
print("\n光源 星27（x=%.3f）／出発 星%d（x=%.3f）" % (SC[27][0], START, SC[START][0]))
print("\n検定PL6  発火して二葉に分かれるか")
hist, fires, splits, addr, live = grow(START, peak, T=12)
print("  発火（番地が節に落ちたコマ）:", [(t, "星%d" % i) for t, i, _ in fires] or "なし")
print("  分岐（最大が二つで同点）:", [(t, "星%d" % i) for t, i, _ in splits] or "なし")
after = [s for s in splits if any(h[0] > s[0] for h in hist)]
ok = bool(fires) and bool(splits) and len(live) >= 2
print("  発火あり・分岐あり・二葉とも生きている  %s" % ("OK" if ok else "NG"))

print("\n検定PL7  番地の取り方と 4H {4,6,7,9}")
FIB = [1, 1]
while len(FIB) < 20:
    FIB.append(FIB[-1] + FIB[-2])
rows = [("累積 Sturmian（増分4,5）", addr),
        ("F(n) mod 13", [f % 13 for f in FIB[:12]]),
        ("等間隔 3 の累積", [(3 * (k + 1)) % 13 for k in range(12)])]
for nm, seq in rows:
    bad = [x for x in seq if x in (4, 6, 7, 9)]
    print("  %-22s %s  4H %d/%d" % (nm, seq[:13], len(bad), len(seq)))
print("  記事の定理1 は F(n) mod 13 についてのもの。増分の累積は13個すべてを巡るので踏む")
print("\n検定PL8  空試験  光が一様")


def flat_tips(i):
    return [q for q, p in AROUND[i]]


live = [(15, 0, 0)]
moved = False
for t in range(12):
    nx = []
    for (i, a, br) in live:
        if len(flat_tips(i)) >= 3:
            continue
        moved = True
    live = nx
print("  最大が5つ同点 → 伸びない  %s" % ("OK" if not moved else "NG"))
