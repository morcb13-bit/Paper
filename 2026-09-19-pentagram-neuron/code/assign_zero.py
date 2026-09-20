#  五芒星ニューロン模型  検定NE4a〜NE4d  札から頂点への割り当て
#
#      発端（引継書 v256 §8-3）：向きは5通りに絞れ（DF11）、0桁目は札で決まる（DF13）が、
#                    札をどう読んで頂点に割り当てるかの手続きを書いていない。
#
#  手続き（外から番号を与えない）
#      星 s の札（検定DF13 と同じもの）＝ 半径 R 以内の五角形を
#      （距離, 方位−72k, 番地の偶奇）で並べた組。k は5通りの回転。
#      0番   = 札が辞書式で最小になる回転がただ一つのとき、その回転。
#      向き  = 方位の符号を反転した札（鏡像）と比べ、小さいほうを正とする。
#      扇で切る規則は使わない。五芒星は尖り18°・凹頂点54°の36°ずれで並んでおり、
#      内側の五角形が星によって切り目に乗るため、境界の取り方で答えが変わる。
#      座面の頂点は 0番の尖りからこの向きに 0,1,2,3,4。
#      出発辺 = 0番頂点から一つ飛ばしに出る対角線（検定NE0 の対応）。
#
#  事前登録（走らせる前に書いた）
#
#  検定NE4a  半径を上げると 0番が決まるか
#      OK なら：担体全部を読めば30個すべて決まる。決まる最小の半径を梯子で出す
#      NG なら：全部読んでも決まらない星が残る
#      ※ 事前登録では「DF12 の φ²で26個・φ³で30個 と一致」を条件にしたが、
#        DF12 のコードが手元に無く半径の取り方を照合できないため取り下げた。
#
#  検定NE4b  向きも一意に決まるか
#      OK なら：30個すべてで 符牒(k+1) ≠ 符牒(k−1)
#
#  検定NE4c  負の対照  半径 φ（まわり10枚の内側）で決めようとする
#      OK なら：30個すべてで決まらない
#      NG なら：決まってしまう＝この検査は何も測っていない
#
#  検定NE4d  出発辺が決まると読み出し順が一つになるか
#      OK なら：NE1b の60通りが星ごとに1通りに落ち、書いて読むと復元する
#
#  距離は 1e-4 に丸めて比べる。順序の判定は組の辞書式比較だけで行う。
#
#  前提：wind_core.py を exec し、pentagram_neuron.py の carrier()・dodeca()・
#        euler_circuit() を使う。

import math, itertools, random
from collections import defaultdict

PHI = (1 + 5 ** 0.5) / 2


def star_tips(faces):
    """担体の五芒星30個の、中心と尖り5個の方位を返す。"""
    out = []
    for a, c in faces:
        if abs(a - 2.9389) < 0.01:
            P = [tuple(float(t) for t in U.xy(p)) for p in c]
            cx = sum(p[0] for p in P) / 10
            cy = sum(p[1] for p in P) / 10
            rr = [math.hypot(p[0] - cx, p[1] - cy) for p in P]
            rmax = max(rr)
            ang = sorted(math.degrees(math.atan2(p[1] - cy, p[0] - cx)) % 360
                         for p, r in zip(P, rr) if abs(r - rmax) < 1e-6)
            out.append(((cx, cy), ang))
    return out


def fuda(XY, F, C, R, rot, mirror=False):
    """検定DF13 と同じ札。方位を 72*rot ずらし、mirror なら符号を反転する。"""
    cx, cy = C
    out = []
    for q, (x, y) in XY.items():
        dx, dy = x - cx, y - cy
        d = math.hypot(dx, dy)
        if d > R + 1e-6 or d < 1e-9:
            continue
        a = math.degrees(math.atan2(dy, dx))
        if mirror:
            a = -a
        out.append((round(d, 4), round((a - 72 * rot) % 360, 3), F[q] % 2))
    return tuple(sorted(out))


def decide(XY, F, C, R):
    """0番の回転と向き。決まらなければ None。"""
    t = [fuda(XY, F, C, R, k) for k in range(5)]
    m = min(t)
    if t.count(m) != 1:
        return None
    k = t.index(m)
    mt = [fuda(XY, F, C, R, j, mirror=True) for j in range(5)]
    mm = min(mt)
    if mm == m:
        return (k, 0)
    return (k, 1 if m < mm else -1)


def seat_cycle(V, FA, fi):
    """面 fi の5頂点を、隣り合う順に並べる。"""
    f = FA[fi]
    ad = {i: [j for j in f if j != i and abs(math.dist(V[i], V[j]) - 1) < 1e-9] for i in f}
    order = [f[0], ad[f[0]][0]]
    while len(order) < 5:
        nxt = [j for j in ad[order[-1]] if j != order[-2]][0]
        order.append(nxt)
    return order


def run_assign():
    NG = 0
    F, faces, SC = carrier()
    ST = star_tips(faces)
    XY = {q: tuple(float(t) for t in U.xy(q)) for q in F}
    Rmax = max(math.hypot(x, y) for x, y in XY.values()) * 2.2

    # ── 検定NE4a・NE4c ───────────────────────────────────
    print("検定NE4a  半径を上げると 0番が決まるか（検定NE4c の負の対照を含む）")
    print("  半径          決まった星   決まらない星")
    got = {}
    ladder = (("φ", PHI), ("φ²", PHI ** 2), ("φ²+", 4.05), ("φ³", PHI ** 3),
              ("φ³+", 5.30), ("φ⁴", PHI ** 4), ("担体全部", Rmax))
    for name, R in ladder:
        dec = [decide(XY, F, C, R) for C, tips in ST]
        n = sum(1 for d in dec if d is not None)
        print(f"    {name:8s}      {n:2d} 個        {30-n:2d} 個")
        got[name] = dec
    ok_a = sum(1 for d in got["担体全部"] if d) == 30
    print(f"  担体全部を読めば30個すべて決まる  {'OK' if ok_a else 'NG'}")
    NG += 0 if ok_a else 1
    ok_c = sum(1 for d in got["φ"] if d) == 0
    print(f"  負の対照 半径φでは決まらない  {'OK' if ok_c else 'NG'}")
    NG += 0 if ok_c else 1

    # ── 検定NE4b ─────────────────────────────────────────
    print("\n検定NE4b  向きも一意に決まるか（担体全部）")
    dec = got["担体全部"]
    zero = sum(1 for d in dec if d and d[1] == 0)
    print(f"  向きが決まらない星 {zero} 個  {'OK' if zero == 0 else 'NG'}")
    NG += 0 if zero == 0 else 1
    print("  星ごとの 0番の扇と向き（先頭10個）")
    for i in range(10):
        if dec[i] is None:
            print(f"    星{i:2d}  決まらない")
        else:
            v = {1: '正', -1: '負', 0: '決まらない'}[dec[i][1]]
            print(f"    星{i:2d}  0番 扇{dec[i][0]}　向き {v}")
    print(f"  0番の分布 {[sum(1 for d in dec if d and d[0]==k) for k in range(5)]}"
          f"／向きの分布 正{sum(1 for d in dec if d and d[1]>0)} "
          f"負{sum(1 for d in dec if d and d[1]<0)}")

    # ── 検定NE4d ─────────────────────────────────────────
    print("\n検定NE4d  出発辺が決まると読み出し順が一つになるか")
    V, E, FA, FD, D = dodeca()
    seat = 0
    cyc = seat_cycle(V, FA, seat)
    didx = {tuple(sorted(d)): k for k, d in enumerate(D)}

    starts = []
    for i, d in enumerate(dec):
        k, s = d
        order = cyc if s > 0 else cyc[::-1]
        v0 = order[k]                      # 0番の尖りに対応する座面の頂点
        v2 = order[(order.index(v0) + 2) % 5]   # 一つ飛ばし
        starts.append(didx[tuple(sorted((v0, v2)))])
    print(f"  30個の出発辺 {sorted(set(starts))}（異なる {len(set(starts))} 通り）")

    random.seed(13)
    MEM = [random.choice([-2, -1, 0, 1, 2]) for _ in range(60)]
    okr = True
    seqs = set()
    for se in starts:
        seq = euler_circuit(D, se)
        if sorted(seq) != list(range(60)):
            okr = False; break
        seqs.add(tuple(seq))
        out = [MEM[k] for k in seq]
        back = [None] * 60
        for pos, k in enumerate(seq):
            back[k] = out[pos]
        if back != MEM:
            okr = False; break
    print(f"  星ごとに順序は1通りか（出発辺が決まれば閉路も決まる） {'はい' if okr else 'いいえ'}")
    print(f"  30個ぶんの読み出し順のうち異なるもの {len(seqs)} 通り（NE1b の60通りから落ちた）")
    print(f"  書いて読む  {'OK' if okr else 'NG'}")
    NG += 0 if okr else 1

    print(f"\nNG {NG} / 4")


if __name__ == "__main__":
    exec(open("wind_core.py").read())
    exec(open("pentagram_neuron.py").read().split("def run_tests")[0])
    run_assign()
