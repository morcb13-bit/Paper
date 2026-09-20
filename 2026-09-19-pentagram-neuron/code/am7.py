#  検定AM7  逆フィボナッチで寄る（濃さが上がらなければ段を落とす）
#
#      検定AM6 で、降りる条件が「動けない」だけだと降りないことが出た。大きい段でも
#      向きに合う隣が必ずあるので跳び続ける。落とす条件をセンサー自身の量で与える。
#
#      一歩   いまの段で、符牒の指す向きにいちばん近い隣へ跳ぶ。
#             跳んだ先の濃さ（まわり5枚の最大）が上がっていなければ、戻って段を落とす。
#             段は戻さない。
#      不変   センサー（まわり5枚の最大がただ一つ）と方向選択（向きに近い隣・90°未満）
#
#      予測A  降順で到達する。コマ数は昇順（検定AM6）より減る
#      予測B  軸の上の星0 は同点のまま
#      負の対照  濃さを見ずに落とす＝検定AM6 降順（振動する。既に出ている）

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
    pairs[round(math.dist(SC[i], SC[j]), 6)].append((i, j))
RUNG = sorted(pairs, reverse=True)
LK = {}
for d in RUNG:
    lk = defaultdict(list)
    for i, j in pairs[d]:
        lk[i].append(j); lk[j].append(i)
    LK[d] = lk


def conc(i, peak):
    """星 i の濃さ＝まわり5枚のうち最大の値（−距離²）。"""
    return max(-((p[0] - peak[0]) ** 2 + (p[1] - peak[1]) ** 2) for q, p in AROUND[i])


def climb(start, peak, order, T=60):
    i, k, steps, drops = start, 0, 0, 0
    while steps < T and k < len(order):
        t = sense(AROUND, i, peak)
        if t is None:
            return i, steps, drops, "センサー同点"
        j = move(SC, LK[order[k]], i, t)
        if j == i or conc(j, peak) <= conc(i, peak):
            k += 1; drops += 1
            continue
        i = j; steps += 1
        if math.dist(SC[i], peak) < 1e-9:
            return i, steps, drops, "到達"
    return i, steps, drops, "行き先なし" if k >= len(order) else "止まらない"


peak = SC[27]
print("段 %d 段（星どうしの距離の種類）／いちばん大きい段 %.3f／小さい段 %.3f"
      % (len(RUNG), RUNG[0], RUNG[-1]))
print("\n  降順＋濃さで落とす")
print("  出発   跳んだ回数   落とした回数   結果        終わりの距離")
for s in (0, 5, 10, 15, 20, 25):
    i, n, dr, why = climb(s, peak, RUNG)
    print("   星%2d      %3d          %3d      %-10s %8.3f"
          % (s, n, dr, why, math.dist(SC[i], peak)))

print("\n  くらべる先  昇順（検定AM6・小さい段から・濃さは見ない）")
print("  出発   コマ数   結果        終わりの距離")


def ladder(start, peak, order, T=60):
    i, k, n = start, 0, 0
    while n < T and k < len(order):
        t = sense(AROUND, i, peak)
        if t is None:
            return i, n, "センサー同点"
        j = move(SC, LK[order[k]], i, t)
        if j == i:
            k += 1
            continue
        i = j; n += 1
        if math.dist(SC[i], peak) < 1e-9:
            return i, n, "到達"
    return i, n, "行き先なし" if k >= len(order) else "止まらない"


for s in (0, 5, 10, 15, 20, 25):
    i, n, why = ladder(s, peak, RUNG[::-1])
    print("   星%2d    %3d    %-10s %8.3f" % (s, n, why, math.dist(SC[i], peak)))

print("\n  山を別の星に置いて確かめる（降順＋濃さ／出発は星10）")
print("   山     跳んだ回数   結果")
for pk in (0, 5, 15, 20, 25, 27, 29):
    i, n, dr, why = climb(10, SC[pk], RUNG)
    print("   星%2d      %3d       %s" % (pk, n, why))
