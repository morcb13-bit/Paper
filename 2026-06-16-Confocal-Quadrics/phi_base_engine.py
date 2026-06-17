#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phi_base_engine.py  —  B13/Bp v147 homework (1)
計量層（共焦点 λ 根・無理数）を整数桁で読み下す装置。

部品:
  - Fibonacci 共焦点族 (a^2,b^2,c^2)=(21,13,8) の λ 三根ソルバ
  - 平衡φ進  digits {-1,0,1}      (黄金, place value = phi^k)
  - 平衡5進  digits {-2,-1,0,1,2} (place value = 5^k)
  - 符号反転 = 負号 の確認（balanced 対称）
  - 逆フィボナッチ重み = φ桁値 × √5 の確認
  - 誤差分散 = {n·alpha} 星状不一致（φ 最小 = 最均等）

Nature knows only addition. 有限の桁でも、無限に続く桁でも。
依存: mpmath のみ。
"""
from mpmath import mp, mpf, sqrt, polyroots, fabs, nint, e

mp.dps = 60
PHI = (1 + sqrt(5)) / 2
PSI = (1 - sqrt(5)) / 2            # -1/phi


# ---------------------------------------------------------------- 共焦点 λ 根
def confocal_lambda(x0, y0, z0, A=mpf(21), B=mpf(13), C=mpf(8)):
    """x^2/(A+L)+y^2/(B+L)+z^2/(C+L)=1 の三根 (降順).
    既定 (A,B,C)=(21,13,8)=(F8,F7,F6); 退化 L=-8,-13,-21; 軸差 8,5,13 (全Fib)."""
    x2, y2, z2 = mpf(x0) ** 2, mpf(y0) ** 2, mpf(z0) ** 2
    s1, s2, s3 = A + B + C, A * B + A * C + B * C, A * B * C
    P2 = x2 + y2 + z2
    P1 = x2 * (B + C) + y2 * (A + C) + z2 * (A + B)
    P0 = x2 * B * C + y2 * A * C + z2 * A * B
    # [P2 L^2+P1 L+P0] - [L^3+s1 L^2+s2 L+s3] = 0
    return sorted(polyroots([-1, P2 - s1, P1 - s2, P0 - s3]), reverse=True)


# ------------------------------------------------------- 平衡φ進 / 平衡5進
def phi_base(x, kmax=14, kmin=-50):
    """平衡φ進展開. 戻り: (digits dict {power:digit in -1..1}, residual)."""
    x = mpf(x)
    dg = {}
    for k in range(kmax, kmin - 1, -1):
        w = PHI ** k
        d = max(-1, min(1, int(nint(x / w))))
        dg[k] = d
        x = x - d * w
    return dg, x


def base5_balanced(x, kmax=4, kmin=-32):
    """平衡5進展開. digits in {-2,-1,0,1,2}."""
    x = mpf(x)
    dg = {}
    for k in range(kmax, kmin - 1, -1):
        w = mpf(5) ** k
        d = max(-2, min(2, int(nint(x / w))))
        dg[k] = d
        x = x - d * w
    return dg, x


def value(dg, base):
    return sum(d * mpf(base) ** k for k, d in dg.items())


def value_phi(dg):
    return sum(d * PHI ** k for k, d in dg.items())


# ---------------------------------------------------------------- 検査群
def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def star_discrepancy(alpha, N=4000):
    """{n*alpha} mod 1 の星状不一致 D*. 低いほど誤差が均等分散."""
    alpha = mpf(alpha)
    pts = sorted(float((n * alpha) % 1) for n in range(1, N + 1))
    D = 0.0
    for i, x in enumerate(pts):
        D = max(D, abs((i + 1) / N - x), abs(x - i / N))
    return D


def report(point=(2, 1, 1)):
    roots = confocal_lambda(*point)
    print(f"# Fibonacci confocal (a2,b2,c2)=(21,13,8), point {point}")
    for i, L in enumerate(roots, 1):
        print(f"  lambda_{i} = {mp.nstr(L, 30)}")

    print("\n# 平衡φ進 / 平衡5進 の収束（残差）")
    for i, L in enumerate(roots, 1):
        dgp, _ = phi_base(L, kmin=-50)
        dg5, _ = base5_balanced(L, kmin=-32)
        ep = fabs(L - value_phi(dgp))
        e5 = fabs(L - value(dg5, 5))
        print(f"  lambda_{i}: phi-base err≈{mp.nstr(ep,4):>11} (down phi^-50) | "
              f"base5 err≈{mp.nstr(e5,4):>11} (down 5^-32)")

    print("\n# 符号反転 = 負号 (balanced 対称)")
    dg, _ = phi_base(roots[0], kmin=-30)
    dgn, _ = phi_base(-roots[0], kmin=-30)
    print("  | (-1)*expand(L) - expand(-L) | =",
          mp.nstr(fabs(-value_phi(dg) - value_phi(dgn)), 4))

    print("\n# 逆フィボナッチ重み = φ桁値 × √5")
    for n in (8, 12, 16):
        r = (mpf(1) / fib(n)) / (PHI ** (-n))
        print(f"  n={n:>2}: (1/F_n)/(phi^-n) = {mp.nstr(r,10)}  -> sqrt5={mp.nstr(sqrt(5),10)}")

    print("\n# 誤差分散: {n*alpha} 星状不一致 D* (低 = 均等), N=4000")
    for name, a in [("phi", PHI), ("sqrt2", sqrt(2)),
                    ("e-2", e - 2), ("near-rational", mpf("0.142857"))]:
        print(f"  alpha={name:<14} D*≈ {star_discrepancy(a):.5f}")


if __name__ == "__main__":
    report()
