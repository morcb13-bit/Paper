#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""#346 の反例を、比が振動する列の中に探す。

なぜ振動列か。既出の Woett の定理により、strongly complete な列は
a[n+1]/a[n] < φ が無限回起きる。したがって**比が一定で φ を超える列は
最初から候補にならない**（測るまでもない）。

残る形はこうなる ── 比が振動し、低い側で φ を下回りながら、
高い側で φ を超え、しかも三条件を満たす列があるか。
そういう列があれば比は収束しないので、#346 は反例を持つ。

  検定1 幾何平均を φ に固定した振動列は strongly complete にならない
      OK なら：Woett の定理と整合。走査の出発点が正しい
      NG なら：判定器か列の作り方が壊れている
  検定2 走査が「両方 True」を一件も返さない
      OK なら：掃いた範囲に反例は無い（**反例が無いことの証明ではない**）
      NG なら：候補が出た。疎な抜きを含めて再検査してから報告すること
  検定3 判定器が「両方 True」を返し得ることを、作りものの列で確かめる
      OK なら：検定2 の OK は反証不能な検査の産物ではない
      NG なら：検定2 の結果に意味が無い。判定器を直すこと

出力の読み方：**掃いた範囲だけが結果である**（判別法D）。
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from completeness import (PHI, is_strongly_complete, breaks_on_infinite_removal,
                          is_lacunary, graham)

R_LOW = (1.20, 1.30, 1.40, 1.45, 1.50, 1.55, 1.60)
R_HIGH = (1.30, 1.45, 1.60, 1.70, 1.80, 2.00)


def periodic(muls, a0, n):
    """倍率を周期的に掛けて作る整数列。"""
    s, i = [a0], 0
    while len(s) < n:
        v = int(round(s[-1] * muls[i % len(muls)]))
        s.append(max(v, s[-1] + 1))
        i += 1
    return s


def check(seq, depth=4, terms=18):
    lac, lo = is_lacunary(seq)
    sc, bad = is_strongly_complete(seq, depth, terms)
    inf, surv = breaks_on_infinite_removal(seq, terms)
    return lac, lo, sc, bad, inf, surv


def main():
    res = []

    print("検定1 幾何平均を φ に固定した振動列")
    ok1 = True
    for r1 in (1.20, 1.40, 1.55, PHI):
        s = periodic([r1, PHI * PHI / r1], 12, 24)
        _, _, sc, _, _, _ = check(s)
        print(f"    r1={r1:.3f}  r2={PHI * PHI / r1:.3f}   strongly {sc}")
        ok1 &= not sc
    res.append(("1", "幾何平均 φ の振動列は strongly complete にならない", ok1))

    print()
    print("検定2 (r1, r2) 平面の走査")
    print(f"    {'r1':<7}{'r2':<7}{'lac':<7}{'SC':<7}{'INF':<7}")
    hits = []
    for r1 in R_LOW:
        for r2 in R_HIGH:
            if r2 <= r1:
                continue
            s = periodic([r1, r2], 12, 26)
            lac, lo, sc, _, inf, surv = check(s)
            if sc:
                flag = "  <-- 両方 True" if (inf and lac) else ""
                print(f"    {r1:<7.2f}{r2:<7.2f}{lo:<7.2f}{str(sc):<7}{str(inf):<7}{flag}")
                if inf and lac:
                    hits.append((r1, r2, surv))
    res.append(("2", "走査が両方 True を返さない", len(hits) == 0))

    print()
    print("検定3 判定器が両方 True を返し得るか（対照）")
    gr = graham(22)
    lac, lo, sc, bad, inf, surv = check(gr)
    both = lac and inf
    print(f"    Graham 列: lac {lac} / SC {sc}（尾 {bad}）/ INF {inf}")
    res.append(("3", "対照が lac と INF を同時に返す", both))

    print()
    for no, name, ok in res:
        print(f"  検定{no} {name:<42} {'OK' if ok else 'NG'}")
    print()
    print(f"  {sum(o for _, _, o in res)}/{len(res)} OK")
    print()
    print("  掃いた範囲:")
    print(f"    r1 ∈ {R_LOW}")
    print(f"    r2 ∈ {R_HIGH}")
    print("    周期2・初項12・26項・尾は n=0..3・各尾18項")
    print("  この範囲の外については何も言っていない。")


if __name__ == "__main__":
    main()
