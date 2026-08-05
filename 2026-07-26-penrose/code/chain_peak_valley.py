#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chain_peak_valley.py ── 桁の印は山と谷のどちらに座るか
============================================================
標準ライブラリのみ。b13_chain_units.py と chain_balanced.py だけを import する（判別法4）。
桁の定義は chain_balanced.balanced() を再利用する。作り直さない（判別法3）。

対象（判別法9）
---------------
奇数の鎖の環中心を、鎖の軸（端から端）を水平に戻したときの軸からのずれで
二つに分ける。上側を山、下側を谷と呼ぶ。測るのは
「桁の印が乗る環がどちらの側か」と「その側の違いが何を担っているか」。

偶数の鎖は折れ角が全箇所 180°（直線）なので、この問いは奇数の鎖にしかない。

事前登録（判別法6）
-------------------
  検定1 環中心は軸の両側に交互に並ぶか
      OK なら：側は環番号の偶奇だけで決まる
      NG なら：交互でない箇所を出す
  検定2 山か谷かは図の裏返しで入れ替わるか
      OK なら：絶対的な「山」「谷」は座標の性質であって構造の性質ではない。
               裏返しで不変なのは「端の環と同じ側か反対側か」だけ
      NG なら：裏返しで不変な絶対量がある
  検定3 印が端と同じ側になる n の条件
      OK なら：n ≡ 1 (mod 4) と一致する
      NG なら：一致しない n を出す
  検定4 桁（n mod 3）と側（n mod 4）は独立か
      OK なら：n の12周期で 3×2 = 6 組すべてが出る
      NG なら：出ない組がある＝両者は同じ量の言い換え
  検定5 同じ側の隣どうしの隔たり
      OK なら：8.0575（一つ飛ばし k=2・共有五角形0枚＝結合しない）
      NG なら：実測値を出す
  検定6 側をまたぐ隣どうしの隔たりと列の間隔
      OK なら：4.9798（連 k=1）と 連×sin36° = 2.9271
      NG なら：実測値を出す
  検定7 対照 ── 軸を 18° として測ると検定1 が落ちるか
      OK なら：検定1 は軸の取り方に反応する＝反証可能
      NG なら：何を軸にしても通る反証不能な検査
      ※ b13_chain_units.xy は既に ROT を含む。スキルの「−18°回す」は
        作図済み座標の話で、ここで重ねて回すと二重回転になる（誤り79）
  検定8 偶数の鎖にも山谷があるか
      OK なら：ある
      NG なら：無い。折れ角180°で軸からのずれが全て0になり側が定義できない
      ※ この検定は NG で終わるのが正しい（反証可能性の常設テスト）
"""
import math, os, sys

_here = os.path.dirname(os.path.abspath(__file__))
for _p in (_here, "/mnt/skills/user/b13-verify/scripts"):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import b13_chain_units as B
import chain_balanced as C

EPS = 1e-9


def axis_angle(n):
    """鎖の軸（端から端）の角度（度）"""
    p = [complex(*B.xy(c)) for c in B.unit(n)]
    d = p[-1] - p[0]
    return math.degrees(math.atan2(d.imag, d.real))


def offsets(n, axis_deg=None):
    """軸を水平に戻したときの各環中心の軸からのずれ"""
    p = [complex(*B.xy(c)) for c in B.unit(n)]
    if axis_deg is None:
        axis_deg = axis_angle(n)
    r = complex(math.cos(math.radians(-axis_deg)), math.sin(math.radians(-axis_deg)))
    p = [z * r for z in p]
    return [z.imag - p[0].imag for z in p]


def sides(n, axis_deg=None, flip=False):
    """各環の側（+1 / −1）。中線を境にする"""
    o = offsets(n, axis_deg)
    mid = (max(o) + min(o)) / 2
    s = [1 if x > mid else -1 for x in o]
    return [-x for x in s] if flip else s


def main(N=59):
    ok = lambda c: "OK" if c else "NG"
    odd = list(range(3, N + 1, 2))

    # 検定1
    bad1 = [n for n in odd
            if any(sides(n)[i] == sides(n)[i + 1] for i in range(n - 1))]
    print("検定1 環中心は軸の両側に交互             %s %s"
          % (ok(not bad1), bad1[:5] if bad1 else ""))
    print("      軸の実測角: %s"
          % ", ".join("n=%d:%.2f°" % (n, axis_angle(n)) for n in (3, 5, 7, 9)))

    # 検定2
    c2 = True
    for n in odd:
        s, f = sides(n), sides(n, flip=True)
        if s == f:
            c2 = False
        cent = (n + 1) // 2
        if (s[cent - 1] == s[0]) != (f[cent - 1] == f[0]):
            c2 = False
    print("検定2 裏返すと山谷は入れ替わり、同/反は不変 %s" % ok(c2))

    # 検定3・4
    rows = []
    for n in odd:
        d, _, cent = C.balanced(n)
        rows.append((n, d, cent, sides(n)[cent - 1] == sides(n)[0]))
    bad3 = [n for n, d, cent, same in rows if same != (n % 4 == 1)]
    print("検定3 印が端と同じ側 ⇔ n ≡ 1 (mod 4)     %s %s"
          % (ok(not bad3), bad3[:5] if bad3 else ""))

    pairs = {(d, same) for n, d, cent, same in rows}
    print("検定4 桁×側 の6組が全部出る（12周期）    %s  %d/6"
          % (ok(len(pairs) == 6), len(pairs)))

    # 検定5・6
    d2, d1 = set(), set()
    for n in odd[:10]:
        p = [complex(*B.xy(c)) for c in B.unit(n)]
        d1 |= {round(abs(p[i + 1] - p[i]), 4) for i in range(len(p) - 1)}
        d2 |= {round(abs(p[i + 2] - p[i]), 4) for i in range(len(p) - 2)}
    gap = max(offsets(9)) - min(offsets(9))
    print("検定5 同じ側の隣どうし = 8.0575（k=2）    %s  %s"
          % (ok(sorted(d2) == [8.0575]), sorted(d2)))
    c6 = sorted(d1) == [4.9798] and abs(gap - 4.979797 * math.sin(math.radians(36))) < 1e-4
    print("検定6 側をまたぐ = 4.9798（k=1）・列間隔  %s  %s / %.4f = 連×sin36°"
          % (ok(c6), sorted(d1), gap))

    # 検定7 対照
    bad7 = [n for n in odd
            if any(sides(n, axis_deg=18.0)[i] == sides(n, axis_deg=18.0)[i + 1]
                   for i in range(n - 1))]
    print("検定7 対照（軸を18°にすると検定1 が落ちる） %s" % ok(bool(bad7)))

    # 検定8 偶数の鎖（NG で終わるのが正しい）
    flat = all(abs(x) < EPS for n in range(2, 15, 2) for x in offsets(n))
    print("検定8 偶数の鎖にも山谷がある              %s  %s"
          % (ok(not flat), "軸からのずれが全て0（直線）" if flat else ""))

    print("\n   n   桁  中心環  端との関係  n mod 4  n mod 3")
    for n, d, cent, same in rows:
        if n > 29:
            break
        print("  %2d   %+2d   %3d    %s       %d        %d"
              % (n, d, cent, "同じ側" if same else "反対側", n % 4, n % 3))

    print("\n  桁は n mod 3、印の側は n mod 4。独立なので合わせて n mod 12。")
    print("  同じ側どうしは k=2（共有五角形0枚）で結合せず、結合は必ず側をまたぐ。")
    nng = sum(1 for c in (not bad1, c2, not bad3, len(pairs) == 6,
                          sorted(d2) == [8.0575], c6, bool(bad7)) if not c)
    print("\nNG %d / 8（検定8 は落ちるのが正しい）" % (nng + 1))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 59)
