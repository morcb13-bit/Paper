"""
b13_constants_compare.py — 複数の数学定数を 3120 平衡展開で比較

目的:
    π の 3120 平衡展開で a_2 = -720 (C60 総剰余角と一致) が出た観測について、
    これが π 固有の性質か、それとも BASE=3120 自体の性質か(どんな定数でも
    各段に B13 関連の整数が立つ)を切り分ける。

    π, e, φ, √2, ln(2), γ(オイラー定数) を同じ手順で展開し、各定数の係数列
    を並べる。13 倍数や C60 関連整数の出現頻度を集計する。

判定基準:
    - π でのみ 720 が出るなら → π 固有 (円と C60 の幾何接続)
    - 他の定数でも頻繁に出るなら → BASE=3120 の性質
    - 13 倍数の出現率がランダム期待値 (1/13 ≈ 7.7%) と一致するなら偶然
"""

from decimal import Decimal, getcontext
from typing import Callable
import math

BASE = 3120
HALF_BASE = BASE // 2
N_TERMS = 30
PREC = 250


def _nearest_int(d: Decimal) -> int:
    """最近接整数 (b13_pi と同じ実装)."""
    floor_val = int(d)
    frac = d - Decimal(floor_val)
    if frac >= Decimal('0.5'):
        return floor_val + 1
    elif frac <= Decimal('-0.5'):
        return floor_val - 1
    return floor_val


def expand_3120(value: Decimal, n_terms: int = N_TERMS) -> list[int]:
    """任意の Decimal 値を 3120 平衡展開する.

    π の場合と完全に同じ手順を任意の定数に適用する.
    """
    getcontext().prec = PREC
    base = Decimal(BASE)
    coefficients = []
    residual = value
    for _ in range(n_terms):
        shifted = residual * base
        a_n = _nearest_int(shifted)
        coefficients.append(a_n)
        residual = shifted - Decimal(a_n)
    return coefficients


# ============================================================
# 高精度定数 (Decimal で 250 桁)
# ============================================================

def pi_decimal() -> Decimal:
    """Machin: π = 4(4 arctan(1/5) - arctan(1/239))"""
    getcontext().prec = PREC + 20
    def arctan_inv(x_inv):
        x = Decimal(1) / Decimal(x_inv)
        x2 = x * x
        result = Decimal(0)
        term = x
        n = 1
        while True:
            new_r = result + term / Decimal(n)
            if new_r == result:
                break
            result = new_r
            term = -term * x2
            n += 2
        return result
    pi = 4 * (4 * arctan_inv(5) - arctan_inv(239))
    getcontext().prec = PREC
    return +pi


def e_decimal() -> Decimal:
    """e = Σ 1/k!"""
    getcontext().prec = PREC + 20
    result = Decimal(0)
    fact = Decimal(1)
    for k in range(0, PREC + 50):
        if k > 0:
            fact *= k
        new_r = result + Decimal(1) / fact
        if new_r == result:
            break
        result = new_r
    getcontext().prec = PREC
    return +result


def sqrt_decimal(n: int) -> Decimal:
    """Newton 法で √n"""
    getcontext().prec = PREC + 20
    x = Decimal(n)
    s = x
    for _ in range(300):
        s_new = (s + x / s) / 2
        if s_new == s:
            break
        s = s_new
    getcontext().prec = PREC
    return +s


def phi_decimal() -> Decimal:
    """黄金比 φ = (1+√5)/2"""
    getcontext().prec = PREC + 20
    phi = (Decimal(1) + sqrt_decimal(5)) / Decimal(2)
    getcontext().prec = PREC
    return +phi


def ln2_decimal() -> Decimal:
    """ln(2) = Σ (-1)^(k+1) / k を高速化 (arctanh 系)
    ln(2) = 2 arctanh(1/3) = 2 Σ 1/((2k+1) 3^(2k+1))
    """
    getcontext().prec = PREC + 20
    result = Decimal(0)
    x_inv = Decimal(3)
    x = Decimal(1) / x_inv
    x2 = x * x
    term = x
    k = 0
    while True:
        denom = Decimal(2 * k + 1)
        new_r = result + term / denom
        if new_r == result:
            break
        result = new_r
        term *= x2
        k += 1
    result = 2 * result
    getcontext().prec = PREC
    return +result


# ============================================================
# 係数列の集計関数
# ============================================================

def factorize(n: int) -> str:
    if n == 0:
        return "0"
    sign = "-" if n < 0 else ""
    n = abs(n)
    if n == 1:
        return sign + "1"
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    parts = []
    for p, e in sorted(factors.items()):
        parts.append(str(p) if e == 1 else f"{p}^{e}")
    return sign + " × ".join(parts)


def has_13_factor(n: int) -> bool:
    return n != 0 and n % 13 == 0


def is_c60_related(n: int) -> bool:
    """C60 関連の整数 (角度総剰余 720, 頂点数 60, 五角形数 12 など)"""
    abs_n = abs(n)
    return abs_n in {12, 20, 30, 60, 90, 120, 180, 360, 720, 1560, 3120}


# ============================================================
# 実験本体
# ============================================================

def run_experiment():
    constants = {
        "π":   pi_decimal(),
        "e":   e_decimal(),
        "φ":   phi_decimal(),
        "√2":  sqrt_decimal(2),
        "√3":  sqrt_decimal(3),
        "√5":  sqrt_decimal(5),
        "ln2": ln2_decimal(),
    }

    print("=" * 80)
    print("複数定数の 3120 平衡展開 (二段目以降のみ表示, n=2..6)")
    print("=" * 80)
    print(f"{'定数':>5} | " + " | ".join(f"a_{n:>2}" for n in range(2, 7)))
    print("-" * 80)

    all_coefs = {}
    for name, val in constants.items():
        coefs = expand_3120(val)
        all_coefs[name] = coefs
        row = " | ".join(f"{coefs[n-1]:>6}" for n in range(2, 7))
        print(f"{name:>5} | {row}")

    print()
    print("=" * 80)
    print("『720』が現れる段 (絶対値で)")
    print("=" * 80)
    for name, coefs in all_coefs.items():
        hits = [(n + 1, c) for n, c in enumerate(coefs) if abs(c) == 720]
        if hits:
            for n, c in hits:
                print(f"  {name}: a_{n} = {c}")
        else:
            print(f"  {name}: なし")

    print()
    print("=" * 80)
    print("C60 関連整数 (720, 60, 360, 12, ...) が現れる段")
    print("=" * 80)
    for name, coefs in all_coefs.items():
        hits = [(n + 1, c) for n, c in enumerate(coefs[1:], start=1) if is_c60_related(c)]
        if hits:
            for n, c in hits:
                print(f"  {name}: a_{n} = {c}")
        else:
            print(f"  {name}: なし")

    print()
    print("=" * 80)
    print("13 倍数の出現率 (二段目以降 n=2..30, ランダム期待値 1/13 ≈ 7.7%)")
    print("=" * 80)
    print(f"{'定数':>5} | {'13倍数の数':>10} | {'全段数':>8} | {'出現率':>8}")
    print("-" * 80)
    for name, coefs in all_coefs.items():
        # 二段目以降のみカウント
        tail = coefs[1:]
        n_13 = sum(1 for c in tail if has_13_factor(c))
        n_total = len(tail)
        rate = n_13 / n_total * 100
        print(f"{name:>5} | {n_13:>10} | {n_total:>8} | {rate:>7.2f}%")

    print()
    print("=" * 80)
    print("π の係数列 (n=1..15) を再確認 (アンカー値)")
    print("=" * 80)
    pi_coefs = all_coefs["π"]
    for n in range(1, 16):
        c = pi_coefs[n - 1]
        fac = factorize(c)
        marker = ""
        if abs(c) == 720:
            marker = "  ← C60 総剰余角"
        elif c == 9802:
            marker = "  ← 26 × F_14 (F_14=377)"
        print(f"  a_{n:>2} = {c:>8}  ({fac}){marker}")

    print()
    print("=" * 80)
    print("二項確率による帰無仮説検定: 「13 倍数の出現はランダムである」")
    print("=" * 80)
    p_null = 1 / 13
    n_trials = N_TERMS - 1  # 二段目以降
    expected = n_trials * p_null
    print(f"試行数 n = {n_trials} (各定数の二段目〜{N_TERMS}段目)")
    print(f"期待値 = {n_trials} × (1/13) = {expected:.2f}")
    print()
    for name, coefs in all_coefs.items():
        tail = coefs[1:]
        observed = sum(1 for c in tail if has_13_factor(c))
        # 二項分布の両側 p 値 (近似)
        # P(X = k) = C(n,k) p^k (1-p)^(n-k)
        from math import comb
        def binom_pmf(k, n, p):
            return comb(n, k) * (p ** k) * ((1 - p) ** (n - k))
        # 両側 p 値: 期待値より離れた事象の確率の合計
        target_dev = abs(observed - expected)
        p_value = 0.0
        for k in range(n_trials + 1):
            if abs(k - expected) >= target_dev:
                p_value += binom_pmf(k, n_trials, p_null)
        signif = "★" if p_value < 0.05 else ""
        print(f"  {name:>5}: 観測 = {observed:>2}, p値 = {p_value:.3f} {signif}")


if __name__ == "__main__":
    run_experiment()
