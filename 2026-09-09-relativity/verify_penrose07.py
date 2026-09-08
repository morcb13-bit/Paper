#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_penrose07.py ── 引継書 v246 §2 の主張を整数だけで検査する

判定はすべて整数の等式・不等式で行う。割り算も浮動小数も判定には使わない
（√5 は (p + q√5)/2 の形の係数として整数対で扱う）。
表示のためだけに実数へ落とす箇所には「表示専用」と書いてある。

各検定は (名前, 合否, 一行の所見) を返す。負の対照つきの検定は
「NG を返せること」もその場で確かめる。

    $ python3 verify_penrose07.py
"""

from itertools import product

RANGE_SMALL = range(-40, 41)     # 符号の全数確認に使う範囲
RANGE_ZERO  = range(-200, 201)   # N = 0 の全数探索に使う範囲
RANGE_MAT   = range(-30, 31)     # 行列の全数確認に使う範囲

results = []


def check(name, ok, note):
    results.append((name, ok, note))


# ---------------------------------------------------------------- 基本の量

def N(a, b):
    """段の印。観察者側の整理では絶対ノルム／世界間隔にあたる。"""
    return a * a + a * b - b * b


def half(ab):
    """半歩（1回ぶん）  (a, b) -> (b, a+b)"""
    a, b = ab
    return (b, a + b)


def full(ab):
    """一歩（2回ぶん）"""
    return half(half(ab))


def addr(k):
    """段 k（半歩単位）の番地 (F_{k-1}, F_k)"""
    ab = (1, 0)
    for _ in range(k):
        ab = half(ab)
    return ab


def pq(ab):
    """番地 -> 向き (p, q) = (2a + b, b)"""
    a, b = ab
    return (2 * a + b, b)


def half_pq(p, q):
    """半歩を (p, q) の側で書いた形。p, q は同じ偶奇なので整数のまま割れる。"""
    assert (p + 5 * q) % 2 == 0 and (p + q) % 2 == 0
    return ((p + 5 * q) // 2, (p + q) // 2)


def full_pq(p, q):
    """一歩を (p, q) の側で書いた形。"""
    assert (3 * p + 5 * q) % 2 == 0 and (p + 3 * q) % 2 == 0
    return ((3 * p + 5 * q) // 2, (p + 3 * q) // 2)


# ---------------------------------------------------------------- 検定1
# 一歩で N は不変、半歩で N は厳密に −N。

def t1():
    bad_full = bad_half = 0
    for a, b in product(RANGE_SMALL, RANGE_SMALL):
        if N(*full((a, b))) != N(a, b):
            bad_full += 1
        if N(*half((a, b))) != -N(a, b):
            bad_half += 1
    ok = (bad_full == 0 and bad_half == 0)
    check("検定1 一歩でN不変・半歩でN反転", ok,
          f"範囲 {RANGE_SMALL.start}〜{RANGE_SMALL.stop - 1} の全数。"
          f"一歩の例外 {bad_full}、半歩の例外 {bad_half}")

    # 負の対照：3回ぶんを「一歩」と誤って置くと不変にならないはず
    bad = sum(1 for a, b in product(range(-8, 9), range(-8, 9))
              if N(*half(half(half((a, b))))) != N(a, b))
    check("検定1-対照 3回ぶんではNが保たれない", bad > 0,
          f"NG を返せることの確認。食い違い {bad} 件（0 なら検定が壊れている）")


# ---------------------------------------------------------------- 検定2
# N = 0 を満たす番地は原点だけ。

def t2():
    zeros = [(a, b) for a, b in product(RANGE_ZERO, RANGE_ZERO) if N(a, b) == 0]
    ok = (zeros == [(0, 0)])
    check("検定2 光の線に番地が無い", ok,
          f"範囲 {RANGE_ZERO.start}〜{RANGE_ZERO.stop - 1} の全数。"
          f"N=0 の点 {zeros}")

    # 負の対照：N = a² + ab − b² を a² + ab − 6b² に差し替えると
    # (a+3b)(a−2b) と分解するので、原点以外の解が出る（判別式 25 が平方）
    other = [(a, b) for a, b in product(range(-30, 31), range(-30, 31))
             if a * a + a * b - 6 * b * b == 0 and (a, b) != (0, 0)]
    check("検定2-対照 別の形なら原点以外に解が出る", len(other) > 0,
          f"NG を返せることの確認。原点以外の解 {len(other)} 件")


# ---------------------------------------------------------------- 検定3
# 段の表：N = (−1)^k、p はリュカ数、q はフィボナッチ数、p² − 5q² = 4N。

def t3(kmax=12):
    rows, bad = [], 0
    L, F = [2, 1], [0, 1]
    for i in range(2, kmax + 2):
        L.append(L[-1] + L[-2])
        F.append(F[-1] + F[-2])
    for k in range(kmax + 1):
        ab = addr(k)
        p, q = pq(ab)
        n = N(*ab)
        if n != (-1) ** k:
            bad += 1
        if p != L[k] or q != F[k]:
            bad += 1
        if p * p - 5 * q * q != 4 * n:
            bad += 1
        rows.append((k, ab, n, p, q, p * p - 5 * q * q))
    check("検定3 段の表（N=(−1)^k、p=L_k、q=F_k、p²−5q²=4N）", bad == 0,
          f"段 0〜{kmax} で食い違い {bad}")
    return rows


# ---------------------------------------------------------------- 検定4
# 半歩を2回掛けると一歩になる。行列式は −1 と +1。

def t4():
    bad = 0
    for p, q in product(RANGE_MAT, RANGE_MAT):
        if (p + q) % 2:            # p, q は同じ偶奇の対だけを見る
            continue
        if half_pq(*half_pq(p, q)) != full_pq(p, q):
            bad += 1
    det_half = (1 * 1 - 5 * 1)     # (1/2)[[1,5],[1,1]] の 4·det
    det_full = (3 * 3 - 5 * 1)     # (1/2)[[3,5],[1,3]] の 4·det
    ok = (bad == 0 and det_half == -4 and det_full == 4)
    check("検定4 半歩×2 = 一歩、行列式 −1 と +1", ok,
          f"範囲 {RANGE_MAT.start}〜{RANGE_MAT.stop - 1} の全数で食い違い {bad}。"
          f"4·det = {det_half} と {det_full}（4 で割って −1 と +1）")


# ---------------------------------------------------------------- 検定5
# 光の間隔に掛かる係数（観察者側でいうボンディの k 係数）が φ^{2m} と一致する。
#
#   β = q√5 / p, γ = p/2 とすると
#   k² = (1+β)/(1−β) = (p + q√5)² / (p² − 5q²) = (p + q√5)² / 4
#   よって k = (p + q√5)/2。 φ^{2m} = (L_{2m} + F_{2m}√5)/2。
#   だから「p = L_{2m} かつ q = F_{2m}」が整数のままの判定になる。

def t5(mmax=6):
    L, F = [2, 1], [0, 1]
    for i in range(2, 2 * mmax + 2):
        L.append(L[-1] + L[-2])
        F.append(F[-1] + F[-2])
    bad = 0
    for m in range(mmax + 1):
        p, q = pq(addr(2 * m))
        if (p, q) != (L[2 * m], F[2 * m]):
            bad += 1
    check("検定5 光の間隔の係数 = φ^{2m}", bad == 0,
          f"一歩 0〜{mmax} で食い違い {bad}（判定は p=L_2m・q=F_2m の整数一致）")


# ---------------------------------------------------------------- 検定6
# ローレンツ収縮は一様。同時刻で測ると、前も後ろも同じ割合で縮む。
#
#   収縮後の位置は x/γ = 2x/p。割り算を使わずに
#   「どの2枚を取っても x_i' · x_j == x_j' · x_i」で一様さを判定する。

def t6(mmax=4):
    boards = [-3, -2, -1, 0, 1, 2, 3]
    bad = 0
    for m in range(1, mmax + 1):
        p, q = pq(addr(2 * m))
        # x' = 2x/p を通分して 2x のまま比べる（分母 p は共通）
        moved = [2 * x for x in boards]
        for (xi, mi), (xj, mj) in product(zip(boards, moved), repeat=2):
            if xi * mj != xj * mi:
                bad += 1
    check("検定6 ローレンツ収縮は一様（前後で違わない）", bad == 0,
          f"一歩 1〜{mmax}、板 {boards} の全対で食い違い {bad}")

    # 負の対照：図5がやっていた「前は割り後ろは掛ける」規則は一様でない
    p, q = pq(addr(4))                     # 一歩2 の係数
    skew = [(2 * x if x >= 0 else x * p * p) for x in boards]
    bad2 = sum(1 for (xi, mi), (xj, mj) in product(zip(boards, skew), repeat=2)
               if xi * mj != xj * mi)
    check("検定6-対照 前後で違う規則は一様にならない", bad2 > 0,
          f"NG を返せることの確認。食い違い {bad2} 件")


# ---------------------------------------------------------------- 検定7
# 列の向き (p, q) は半歩ごとに 45° を跨ぐ。
#
#   45° との比較は tan の比較で足りる：
#   時間側の成分 q√5/2 と空間側の成分 p/2 を比べる ⟺ 5q² と p² を比べる。
#   p² − 5q² = ±4 なので、符号がそのまま内外になる。

def t7(kmax=10):
    bad = 0
    prev = None
    for k in range(kmax + 1):
        p, q = pq(addr(k))
        inside = (p * p > 5 * q * q)       # 45° より下＝空間の並び
        if inside != (N(*addr(k)) > 0):
            bad += 1
        if prev is not None and inside == prev:
            bad += 1
        prev = inside
    check("検定7 半歩ごとに 45° を跨ぐ", bad == 0,
          f"段 0〜{kmax} で食い違い {bad}（判定は p² と 5q² の大小のみ）")


# ---------------------------------------------------------------- 検定8
# ものさしの幅。曲線と光の線のあいだを8等分したときの1本ぶんは 4/p²。
# 判定は「p² が段ごとに真に増える」ことだけで足りる（幅は 4/p² なので単調に縮む）。

def t8(mmax=5):
    ps = [pq(addr(2 * m))[0] for m in range(mmax + 1)]
    ok = all(ps[i] * ps[i] < ps[i + 1] * ps[i + 1] for i in range(len(ps) - 1))
    check("検定8 ものさしの幅 4/p² は単調に縮む", ok,
          f"p = {ps}、p² = {[x * x for x in ps]}")


# ---------------------------------------------------------------- 表示

def show_table(rows):
    print()
    print("段の表（検定3）")
    print(" k  (a, b)      N   p    q    p²−5q²   側")
    print(" " + "-" * 46)
    for k, ab, n, p, q, d in rows:
        side = "内" if n > 0 else "外"
        print(f"{k:2d}  ({ab[0]:3d},{ab[1]:4d})  {n:+d}  {p:4d} {q:4d}   {d:+d}    {side}")


def show_speeds(mmax=4):
    """表示専用。判定には使わない。"""
    from math import sqrt
    r5 = sqrt(5)
    print()
    print("表示専用（実数に落としたもの・判定には使っていない）")
    print(" 一歩  γ = p/2   β = q√5/p     光の間隔の比 (p+q√5)/2")
    print(" " + "-" * 52)
    for m in range(mmax + 1):
        p, q = pq(addr(2 * m))
        print(f" {m:3d}  {p/2:8.4f}  {q*r5/p:10.6f}  {(p + q*r5)/2:16.6f}")


def main():
    t1()
    t2()
    rows = t3()
    t4()
    t5()
    t6()
    t7()
    t8()

    print("=" * 62)
    print("verify_penrose07.py ── 引継書 v246 §2 の検査")
    print("=" * 62)
    for name, ok, note in results:
        print(f"[{'OK' if ok else 'NG'}] {name}")
        print(f"     {note}")
    show_table(rows)
    show_speeds()

    ng = [n for n, ok, _ in results if not ok]
    print()
    print("=" * 62)
    if ng:
        print(f"NG {len(ng)} 件： " + " / ".join(ng))
        return 1
    print(f"全 {len(results)} 件 OK（対照を含む）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
