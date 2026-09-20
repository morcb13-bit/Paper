#  検定AM5  逆フィボナッチで寄る
#
#      担体に入っている φ の冪の距離は二段だけ（φ⁴ = 6.8541 が35本、φ⁶ = 17.9443 が6本）。
#      大きい段から使い、動けなくなったら段を落とす。段は戻さない。
#
#      予測A  AM4 の固定の隣より、到達までのコマ数が減るか同じ
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

RUNG = [("φ⁶", PHI ** 6), ("φ⁴", PHI ** 4)]


def links_at(L):
    lk = defaultdict(list)
    for i, j in itertools.combinations(range(30), 2):
        if abs(math.dist(SC[i], SC[j]) - L) < 2e-3:
            lk[i].append(j); lk[j].append(i)
    return lk


LK = {nm: links_at(L) for nm, L in RUNG}
for nm, L in RUNG:
    print("段 %s (%.4f)  %d 本" % (nm, L, sum(len(v) for v in LK[nm].values()) // 2))


def ladder(start, peak, order, T=40):
    """order の順に段を使う。動けなくなったら次の段へ落とす。戻さない。"""
    i, k, log = start, 0, []
    for _ in range(T):
        t = sense(AROUND, i, peak)
        if t is None:
            return i, log, "センサー同点"
        j = move(SC, LK[order[k]], i, t)
        if j == i:
            if k + 1 < len(order):
                k += 1
                continue
            return i, log, "行き先なし"
        log.append(order[k])
        i = j
        if math.dist(SC[i], peak) < 1e-9:
            return i, log, "到達"
    return i, log, "止まらない"


peak = SC[27]
print("\n  出発   使った段の列            コマ数   結果      終わりの距離")
for s in (0, 5, 10, 15, 20, 25):
    i, log, why = ladder(s, peak, ["φ⁶", "φ⁴"])
    print("   星%2d   %-22s %3d     %-10s %8.3f"
          % (s, "".join("6" if x == "φ⁶" else "4" for x in log) or "—",
             len(log), why, math.dist(SC[i], peak)))

print("\n  負の対照  段を昇順（φ⁴ → φ⁶）にする")
print("  出発   使った段の列            コマ数   結果      終わりの距離")
for s in (0, 5, 10, 15, 20, 25):
    i, log, why = ladder(s, peak, ["φ⁴", "φ⁶"])
    print("   星%2d   %-22s %3d     %-10s %8.3f"
          % (s, "".join("6" if x == "φ⁶" else "4" for x in log) or "—",
             len(log), why, math.dist(SC[i], peak)))

print("\n  くらべる先（検定AM4・隣を三番目の距離まで固定）")
lk3 = defaultdict(list)
for a, b in itertools.combinations(range(30), 2):
    if math.dist(SC[a], SC[b]) <= 10.5784 + 1e-3:
        lk3[a].append(b); lk3[b].append(a)
for s in (0, 5, 10, 15, 20, 25):
    i, n = s, 0
    for _ in range(40):
        t = sense(AROUND, i, peak)
        if t is None:
            break
        j = move(SC, lk3, i, t)
        if j == i:
            break
        i = j; n += 1
        if math.dist(SC[i], peak) < 1e-9:
            break
    print("   星%2d   コマ数 %3d   終わりの距離 %8.3f" % (s, n, math.dist(SC[i], peak)))
