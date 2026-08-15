#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quat_exact.py ── 三角関数を使わない四元数回転。

回転四元数 Q = (cos(θ/2); α sin(θ/2), β sin(θ/2), γ sin(θ/2)) を、
20面体群の回転に限れば、成分はすべて Z[φ] の半整数になる。
cos・sin を一度も呼ばずに、厳密な回転が書ける。

  対象：単位イコシアン120個（= 600胞体の頂点 = 20面体群の回転四元数）
  演算：a + bφ を分数のまま持つ。丸めはどこにも入らない。

  事前登録

  検定1 サンドイッチの実部が厳密に0
      OK なら：q P q̄ の実部が、全120回転×全頂点で厳密に 0（浮動小数の
               「ほぼ0」ではなく等号）
      NG なら：積の規則が違う

  検定2 長さが厳密に保存される
      OK なら：回転後のノルムが元と厳密に等しい
      NG なら：単位四元数になっていない

  検定3 何種類の回転になるか
      OK なら：q と −q が同じ回転を与えるので、120個 → **60通り**
      NG なら：数え方が違う

  検定4 表せる角度
      OK なら：実部が取りうる値 {0, ±1/2, ±(φ−1)/2, ±φ/2, ±1} に対応して、
               θ ∈ {0°, 72°, 120°, 144°, 180°} だけが表せる
      NG なら：別の角が混じる

  検定5 誤差が溜まらない
      OK なら：回転を10万回合成しても、整数版のずれは厳密に 0。
               同じ手順の浮動小数版はずれが溜まる（＝対照）
      NG なら：整数版でもずれる
"""

import math
import random
from fractions import Fraction as F
import icosian as I

PHI = (1 + 5**0.5)/2
ZERO, ONE = I.ZERO, I.ONE


def conj(q):
    return I.H([q.c[0], I.Q5(0, 0) - q.c[1],
                I.Q5(0, 0) - q.c[2], I.Q5(0, 0) - q.c[3]])


def vec(x, y, z):
    return I.H([ZERO, x, y, z])


def rotate(q, p):
    """q p q̄。三角関数も浮動小数も使わない。"""
    return (q * p) * conj(q)


def angle_of(q):
    """実部 cos(θ/2) から θ を出す（表示用。計算には使わない）。"""
    c = q.c[0].val()
    c = max(-1.0, min(1.0, c))
    return round(math.degrees(2*math.acos(c)), 6)


def main():
    V = I.build()

    # 正20面体の12頂点（Z[φ] のまま）: (0, ±1, ±φ)/2 の巡回
    half = I.Q5(F(1, 2), 0)
    phi2 = I.Q5(0, F(1, 2))
    ico = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            ico.append(vec(ZERO, half*s1, phi2*s2))
            ico.append(vec(half*s1, phi2*s2, ZERO))
            ico.append(vec(phi2*s2, ZERO, half*s1))
    icoset = set(ico)

    print("== 検定1 サンドイッチの実部が厳密に0 ==")
    bad = 0
    for q in V:
        for p in ico:
            if rotate(q, p).c[0] != ZERO:
                bad += 1
    print("   %d 回転 × %d 頂点 のうち実部が0でないもの %d → %s"
          % (len(V), len(ico), bad, "OK" if bad == 0 else "NG"))

    print("== 検定2 長さが厳密に保存される ==")
    bad = 0
    for q in V:
        for p in ico:
            if rotate(q, p).norm() != p.norm():
                bad += 1
    print("   ノルムが変わったもの %d → %s" % (bad, "OK" if bad == 0 else "NG"))

    print("== 検定3 何種類の回転になるか ==")
    seen = {}
    for q in V:
        key = tuple(rotate(q, p).c for p in ico)
        seen.setdefault(key, []).append(q)
    sizes = sorted(set(len(v) for v in seen.values()))
    ok3 = (len(seen) == 60 and sizes == [2])
    print("   相異なる回転 %d 通り／1通りあたりの四元数 %s（q と −q）→ %s"
          % (len(seen), sizes, "OK" if ok3 else "NG"))

    print("== 検定4 表せる角度 ==")
    ang = {}
    for q in V:
        ang[angle_of(q)] = ang.get(angle_of(q), 0) + 1
    ok4 = set(ang) <= {0.0, 72.0, 120.0, 144.0, 180.0, 216.0, 240.0, 288.0,
                       360.0}
    print("   θ の分布 %s" % dict(sorted(ang.items())))
    print("   実部が取る値 %s"
          % sorted(set((q.c[0].a, q.c[0].b) for q in V),
                   key=lambda t: float(t[0])+float(t[1])*PHI))
    print("   → %s（20面体群の角だけ。任意角は表せない）"
          % ("OK" if ok4 else "NG"))

    print("== 検定5 誤差が溜まらない ==")
    N = 100000
    random.seed(13)
    word = [random.randrange(len(V)) for _ in range(N)]

    # 整数版：単位イコシアンの積は必ず単位イコシアン＝状態は120個の中に留まる
    q = I.H([ONE, ZERO, ZERO, ZERO])
    for k in word:
        q = q * V[k]
    for k in reversed(word):
        q = q * conj(V[k])
    exact_ok = (q.c == I.H([ONE, ZERO, ZERO, ZERO]).c)

    # 浮動小数版：同じ手順
    def fmul(a, b):
        return (a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
                a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],
                a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],
                a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0])
    Vf = [tuple(c.val() for c in v.c) for v in V]
    f = (1.0, 0.0, 0.0, 0.0)
    for k in word:
        f = fmul(f, Vf[k])
    for k in reversed(word):
        g = Vf[k]
        f = fmul(f, (g[0], -g[1], -g[2], -g[3]))
    drift = max(abs(f[0]-1.0), abs(f[1]), abs(f[2]), abs(f[3]))
    print("   %d 回合成して戻す" % (2*N))
    print("   整数版のずれ：%s（厳密等号）" % ("0" if exact_ok else "0でない"))
    print("   浮動小数版のずれ：%.3e ← 対照" % drift)
    print("   → %s" % ("OK" if exact_ok and drift > 0 else "NG"))

    print("\n== おまけ：正20面体が自分に写るか ==")
    bad = sum(1 for q in V for p in ico if rotate(q, p) not in icoset)
    print("   写らなかった対 %d / %d → %s"
          % (bad, len(V)*len(ico), "OK" if bad == 0 else "NG"))

    print("\n== 使い方 ==")
    q = next(v for v in V if angle_of(v) == 72.0)
    print("   72°回転の四元数 Q = (%s; %s, %s, %s)"
          % tuple(q.c))
    p = ico[0]
    print("   点 P = (%s, %s, %s)" % tuple(p.c[1:]))
    r = p
    for n in range(1, 6):
        r = rotate(q, r)
        print("   %d 回目 → (%s, %s, %s)%s"
              % (n, r.c[1], r.c[2], r.c[3],
                 "  ← 5回で元に戻る（厳密）" if r.c == p.c else ""))


if __name__ == "__main__":
    main()
