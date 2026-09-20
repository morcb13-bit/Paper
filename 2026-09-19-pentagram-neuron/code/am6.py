#  検定AM6  黄金整数の段で寄る（逆フィボナッチの細かい梯子）
#
#      星どうしの距離は d² = p + qφ（検定IF1・IF2 の側）。φ の冪は φ⁴ と φ⁶ の二段しか
#      無いが、黄金整数で取れば段はもっと細かい。大きい段から使い、動けなくなったら
#      一段落とす。段は戻さない。
#
#      予測A  段が細かくなるので、AM5（φの冪二段）より到達が増える
#      予測B  軸の上の星0 は変わらず同点で止まる
#      負の対照  段を昇順にすると寄りきらない
#
#      センサーと方向選択は一切変えない。

import math, itertools
from collections import defaultdict

exec(open('wind_core.py').read())
exec(open('pentagram_neuron.py').read().split('def run_tests')[0])
src = open('amoeba.py').read()
exec(src[src.index('def build'):src.index('def run_amoeba')])

PHI = (1 + 5 ** 0.5) / 2
SC, AROUND, _ = build()

pairs = defaultdict(list)
for i, j in itertools.combinations(range(30), 2):
    d = math.dist(SC[i], SC[j])
    pairs[round(d, 6)].append((i, j))
DS = sorted(pairs)
gold = {d: golden(d * d, tol=1e-4) for d in DS}
ng = [d for d in DS if gold[d] is None]
print("星どうしの距離 %d 種／d² が黄金整数でないもの %d 種" % (len(DS), len(ng)))
print("  近い順に8段：")
for d in DS[:8]:
    g = gold[d]
    print("    %8.4f   d² = %s" % (d, "%d+%dφ" % g if g else "—"))

RUNG = sorted(DS, reverse=True)
LK = {}
for d in RUNG:
    lk = defaultdict(list)
    for i, j in pairs[d]:
        lk[i].append(j); lk[j].append(i)
    LK[d] = lk


def ladder(start, peak, order, T=60):
    i, k, log = start, 0, []
    while len(log) < T and k < len(order):
        t = sense(AROUND, i, peak)
        if t is None:
            return i, log, "センサー同点"
        j = move(SC, LK[order[k]], i, t)
        if j == i:
            k += 1
            continue
        log.append(order[k])
        i = j
        if math.dist(SC[i], peak) < 1e-9:
            return i, log, "到達"
    return i, log, "行き先なし" if k >= len(order) else "止まらない"


peak = SC[27]
print("\n  降順（大きい段から）")
print("  出発   コマ数   結果        終わりの距離   使った段（最初の5つ）")
for s in (0, 5, 10, 15, 20, 25):
    i, log, why = ladder(s, peak, RUNG)
    print("   星%2d    %3d    %-10s %8.3f    %s"
          % (s, len(log), why, math.dist(SC[i], peak),
             ["%.3f" % x for x in log[:5]] or "—"))

print("\n  負の対照  昇順（小さい段から）")
print("  出発   コマ数   結果        終わりの距離")
for s in (0, 5, 10, 15, 20, 25):
    i, log, why = ladder(s, peak, RUNG[::-1])
    print("   星%2d    %3d    %-10s %8.3f" % (s, len(log), why, math.dist(SC[i], peak)))
