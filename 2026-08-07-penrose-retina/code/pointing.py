#  ペンタゴン網膜：同じ場所を指せるか
#  監督の指定：これは段数を上げる装置ではなく、同じ画像を見て会話するための
#              インターフェース。中心に寄せて見るのは人間と同じ。
#
#  対象を一行で（判別法9）：ある場所のまわりを半径 r だけ見たとき、
#  その見え（景色）で場所が一意に決まるか。決まれば「そこ」が二人で共有できる。
#
#  検定P1 局所の見えで場所が分かれるか
#      OK なら：見える範囲を広げるほど類が増える＝指せる
#      NG なら：どこも同じ景色＝場所を指す言葉が作れない
#
#  検定P2 一意になる半径はいくつか
#      有限で止まる：その半径が「指すのに要る視野」
#
#  照合には鏡映を入れる（判別法10）。担体の縁で切れた景色は除く（判別法7）。

import json
from fractions import Fraction as F
from qphi import Qp, zmul, zconj, zre, ZERO

Z = (0, 1, 0, 0)
def zneg(a): return tuple(-x for x in a)
def zsub(a, b): return tuple(a[i] - b[i] for i in range(4))
def zrot(a, k):
    b = a
    for _ in range(k % 5): b = zmul(b, Z)
    return b if (k % 10) < 5 else zneg(b)
def norm2(a): return zre(zmul(a, zconj(a)))

d = json.load(open("expanded.json"))
rings = [tuple(F(x) for x in v) for v in d["rings"]]
R_FIG = max(norm2(v).val() ** 0.5 for v in rings)
print(f"担体：環{len(rings)}  実半径 {R_FIG:.4f}")

def canon(nb):
    """回転10通り＋鏡映で畳んだ標準形（厳密）"""
    best = None
    for k in range(10):
        for mirror in (False, True):
            t = tuple(sorted(zconj(zrot(v, k)) if mirror else zrot(v, k)
                             for v in nb))
            if best is None or t < best: best = t
    return best

print()
print("=" * 70)
print(f"{'見る半径':>10s} {'縁で切れない中心':>16s} {'相異なる景色':>14s} {'最大の同一群':>14s}")
rows = []
for r in [4.3, 7.0, 9.2, 11.1, 15.0]:
    ok_centers = [c for c in rings if norm2(c).val() ** 0.5 <= R_FIG - r]
    if len(ok_centers) < 2: continue
    cls = {}
    for c in ok_centers:
        nb = [zsub(v, c) for v in rings
              if v != c and norm2(zsub(v, c)).val() <= r * r]
        cls.setdefault(canon(nb), []).append(c)
    big = max(len(v) for v in cls.values())
    rows.append((r, len(ok_centers), len(cls), big))
    print(f"{r:10.1f} {len(ok_centers):16d} {len(cls):14d} {big:14d}")

grew = rows[-1][2] / max(rows[0][2], 1)
print("検定P1:", "OK  半径を広げると景色が分かれる" if grew > 1 else "NG  分かれない")
uniq = [r for r, n, c, b in rows if b == 1]
print("検定P2:", f"一意になる半径 {uniq[0]}" if uniq
      else f"60環の中では一意にならない（最大の同一群 {rows[-1][3]}）")
print("        担体が小さいので、これは幾何の禁止ではない（判別法7）")
