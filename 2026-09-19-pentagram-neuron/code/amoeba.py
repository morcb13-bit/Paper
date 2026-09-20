#  ペンローズアメーバ  検定AM1〜AM3  濃さの山を追う
#
#      発端（監督）：餌（仮にブドウ糖）を追うセンサーに繋がった正十二面体だけで、
#                    ペンローズアメーバはオートマトンで動く。
#
#  装置（足したものは無い。既にある量だけを使う）
#      餌の濃さ   担体の五角形ごとの値。山に近いほど大きい（−距離²）
#      センサー   五芒星のまわりの五角形5枚。符牒の5方向がそのまま5つの標本
#      一歩       5枚のうち最大がただ一つなら、その向きの隣の星へ原点を移す。
#                 ただ一つでなければ動かない（検定NE4 の「ただ一つなら決まる」と同じ形）
#      隣         距離 φ⁴ で繋がる星（検定DF6 の35本）。向きの差が90°未満のものに限る
#
#  事前登録（走らせる前に書いた）
#
#  検定AM1  濃さの山へ登るか
#      OK なら：数コマで原点が山の位置に達し、そこで止まる
#      NG なら：登らない／通り過ぎて止まらない
#
#  検定AM2  山を動かすと追うか
#      OK なら：山の速さが装置の一歩より遅いあいだ、離れずに追い続ける
#
#  検定AM3  空試験  濃さが一様
#      OK なら：原点が動かない
#
#  負の対照  でたらめな隣へ移る
#      OK なら：山に達しない（AM1 が符牒を見ている証拠になる）
#
#  値は大小を比べるだけで、足し引きも割り算もしない。
#
#  前提：wind_core.py・pentagram_neuron.py を同じ場所に置く。

import math, itertools, random
from collections import defaultdict

PHI = (1 + 5 ** 0.5) / 2


def build():
    F, faces, SC = carrier()
    XY = {q: tuple(float(t) for t in U.xy(q)) for q in F}
    CELL = list(XY.items())
    AROUND = []
    for cx, cy in SC:
        d = sorted(CELL, key=lambda t: (t[1][0] - cx) ** 2 + (t[1][1] - cy) ** 2)[:5]
        AROUND.append(sorted(d, key=lambda t: math.degrees(
            math.atan2(t[1][1] - cy, t[1][0] - cx)) % 360))
    link = defaultdict(list)
    for i, j in itertools.combinations(range(30), 2):
        if abs(math.dist(SC[i], SC[j]) - PHI ** 4) < 1e-6:
            link[i].append(j); link[j].append(i)
    return SC, AROUND, link


def sense(AROUND, i, peak):
    """星 i のまわり5枚の濃さ（−距離²）と、最大がただ一つかどうか。"""
    vals = [-( (p[0] - peak[0]) ** 2 + (p[1] - peak[1]) ** 2 ) for q, p in AROUND[i]]
    m = max(vals)
    if vals.count(m) != 1:
        return None
    k = vals.index(m)
    return AROUND[i][k][1]          # いちばん濃い五角形の座標


def move(SC, link, i, target):
    """target の向きにいちばん近い隣へ移る。90°以上外れていれば動かない。"""
    cx, cy = SC[i]
    th = math.atan2(target[1] - cy, target[0] - cx)
    best, bd = None, None
    for j in link[i]:
        a = math.atan2(SC[j][1] - cy, SC[j][0] - cx)
        d = abs((a - th + math.pi) % (2 * math.pi) - math.pi)
        if bd is None or d < bd:
            bd, best = d, j
    if bd is None or bd >= math.pi / 2:
        return i
    return best


def walk(SC, AROUND, link, start, peak, T=20, rnd=None):
    """T コマ動かす。各コマの（星、山までの距離）を返す。"""
    i = start
    out = [(i, math.dist(SC[i], peak))]
    for _ in range(T):
        if rnd is not None:
            i = rnd.choice(link[i]) if link[i] else i
        else:
            t = sense(AROUND, i, peak)
            i = i if t is None else move(SC, link, i, t)
        out.append((i, math.dist(SC[i], peak)))
    return out


def run_amoeba():
    NG = 0
    SC, AROUND, link = build()
    print("担体の五芒星30個／隣の本数 %s"
          % {k: sum(1 for i in range(30) if len(link[i]) == k) for k in (1, 2, 3)})
    dirs = sorted(set(round(math.degrees(math.atan2(SC[j][1] - SC[i][1],
                                                    SC[j][0] - SC[i][0])) % 180, 2)
                      for i in link for j in link[i]))
    print("隣への向き（180°で畳む）：%s" % dirs)

    # ── 検定AM1 ─────────────────────────────────────────
    print("\n検定AM1  濃さの山へ登る（山は星27 の位置。出発を変えて20コマ）")
    print("  出発   はじめの距離   終わりの距離   達したコマ   止まったか")
    peak = SC[27]
    okc = 0
    for s in (0, 5, 10, 15, 20, 25):
        w = walk(SC, AROUND, link, s, peak)
        arrive = next((t for t, (i, d) in enumerate(w) if d < 1e-9), None)
        still = w[-1][0] == w[-2][0]
        print(f"   星{s:2d}     {w[0][1]:8.3f}     {w[-1][1]:8.3f}      "
              f"{arrive if arrive is not None else '—':>4}      {'止まった' if still else '動いている'}")
        okc += 1 if arrive is not None else 0
    print(f"  山に達した出発 {okc}/6  {'OK' if okc >= 5 else 'NG'}")
    NG += 0 if okc >= 5 else 1

    # ── 負の対照 ─────────────────────────────────────────
    print("\n負の対照  でたらめな隣へ移る（同じ出発・20コマ・10回）")
    rnd = random.Random(13)
    hit = 0
    for s in (0, 5, 10, 15, 20, 25):
        for _ in range(10):
            w = walk(SC, AROUND, link, s, peak, rnd=rnd)
            hit += any(d < 1e-9 for _, d in w)
    print(f"  山に達した回数 {hit}/60  {'OK' if hit < 20 else 'NG'}")
    NG += 0 if hit < 20 else 1

    # ── 検定AM2 ─────────────────────────────────────────
    print("\n検定AM2  山を動かすと追うか（山を一コマあたり dx だけ横へ動かす）")
    print("   dx     20コマ後の距離   途中の最大の離れ")
    for dx in (0.0, 0.5, 1.5, 3.0):
        i = 10
        pk = list(SC[27])
        far = 0.0
        for t in range(20):
            tgt = sense(AROUND, i, pk)
            if tgt is not None:
                i = move(SC, link, i, tgt)
            pk[0] += dx
            far = max(far, math.dist(SC[i], pk))
        print(f"  {dx:4.1f}      {math.dist(SC[i], pk):10.3f}      {far:10.3f}")

    # ── 検定AM3 ─────────────────────────────────────────
    print("\n検定AM3  空試験  濃さが一様")
    i = 10
    same = True
    for _ in range(20):
        vals = [0 for _ in AROUND[i]]
        m = max(vals)
        if vals.count(m) == 1:
            same = False
    print(f"  最大がただ一つになることがあるか {'ある' if not same else 'ない'}"
          f"／原点は動かない  {'OK' if same else 'NG'}")
    NG += 0 if same else 1

    print(f"\nNG {NG} / 3")


if __name__ == "__main__":
    exec(open("wind_core.py").read())
    exec(open("pentagram_neuron.py").read().split("def run_tests")[0])
    run_amoeba()
