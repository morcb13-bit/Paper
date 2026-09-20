#  検定AM4  隣を次の距離まで広げる（センサーと方向選択規則は変えない）
#
#      予測A  次数1の葉の停止が消える
#      予測B  軸の上の星の停止は残る（左右同値は隣を増やしても消えない）
#
#  変えるのは候補となる隣の集合だけ。符牒の読み方（まわり5枚の最大がただ一つ）と、
#  選び方（その向きにいちばん近い隣・90°以上外れていれば動かない）はそのまま。

import math, itertools
from collections import defaultdict

exec(open('wind_core.py').read())
exec(open('pentagram_neuron.py').read().split('def run_tests')[0])
src = open('amoeba.py').read()
exec(src[src.index('def build'):src.index('def run_amoeba')])

PHI = (1 + 5 ** 0.5) / 2
SC, AROUND, _ = build()

ds = sorted(set(round(math.dist(SC[i], SC[j]), 4)
                for i, j in itertools.combinations(range(30), 2)))
print("星どうしの距離（近い順に6つ）:", ds[:6], " φ⁴ = %.4f" % PHI ** 4)


def links_upto(n):
    lk = defaultdict(list)
    lim = ds[n - 1] + 1e-3
    for i, j in itertools.combinations(range(30), 2):
        if math.dist(SC[i], SC[j]) <= lim:
            lk[i].append(j); lk[j].append(i)
    return lk


def why(lk, i, peak):
    """止まった理由。"""
    t = sense(AROUND, i, peak)
    if t is None:
        return "センサー同点"
    return "行き先なし" if move(SC, lk, i, t) == i else "動ける"


peak = SC[27]
print("\n  隣の範囲      本数  次数1の星   到達した出発   止まった理由")
for n in (1, 2, 3):
    lk = links_upto(n)
    nb = sum(len(v) for v in lk.values()) // 2
    leaf = sum(1 for i in range(30) if len(lk[i]) == 1)
    arrive, reasons = 0, []
    for s in (0, 5, 10, 15, 20, 25):
        i, seen, cyc = s, [], None
        for _ in range(60):
            if i in seen:
                cyc = len(seen) - seen.index(i)
                break
            seen.append(i)
            t = sense(AROUND, i, peak)
            if t is None:
                break
            j = move(SC, lk, i, t)
            if j == i:
                break
            i = j
        if math.dist(SC[i], peak) < 1e-9:
            arrive += 1
            reasons.append("星%d 到達（%d コマ）" % (s, len(seen) - 1))
        elif cyc:
            reasons.append("星%d 振動（長さ%d・星%d のあたり）" % (s, cyc, i))
        else:
            reasons.append("星%d %s(星%d)" % (s, why(lk, i, peak), i))
    print("  %s(%.3f)  %3d 本   %2d 個      %d/6" % (("最小", "次", "三番目")[n - 1],
                                                     ds[n - 1], nb, leaf, arrive))
    for r in reasons:
        print("       ", r)
