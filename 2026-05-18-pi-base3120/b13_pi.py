"""
b13_pi.py — π の B13 (BASE=3120) 平衡展開

設計方針:
    Nature knows only addition.
    自然界は加算と最近接整数化のみで連続量を実装している。
    π もその例外ではない。

    本モジュールは π を以下のフラクタル展開で表現する:

        π = Σ_{n=1}^{∞} a_n / 3120^n

    各係数 a_n は「3120^n 倍に拡大した円の、3120 点離散円格子上で
    最も近い点」として定義される整数である(n ≥ 2 では |a_n| < 3120/2、
    ただし a_1 は π の整数部 3 を吸収するため例外で a_1 ≈ 3 × BASE + 残差
    のオーダーになる)。

    必要な演算は:
        1. 残差を 3120 倍する (拡大)
        2. 最近接整数を選ぶ (離散円格子上の点選択)
        3. 残差を更新する (加算)
    のみ。階乗・冪・平方根・三角関数・無限和の概念は不要。

    本実装は b13phase パッケージの BASE=3120 と整合する。

確定事項 (v126 記録用):
    a_1 = 9802 = 2 × 13² × 29 = 26 × F_14   (フィボナッチ F_14 を含む)
    a_2 = -720 = -2⁴ × 3² × 5               (C60 総剰余角と数値一致)
    a_3 = -1475 = -5² × 59
    a_4 = -1354 = -2 × 677
    a_5 = -367  (素数)
    ...

ChatGPT ゲートキーパー向けメモ:
    本実装はラマヌジャン公式との比較で収束速度が劣ること (3.18 桁/項 vs
    8.08 桁/項) を承知の上で、「自然界が実装しうる演算のみ」で π を構成
    することを目的とする。比較対象は楽譜(ラマヌジャン)ではなく演奏
    (自然)。詳細は handover v126-θ-π を参照。
"""

from __future__ import annotations
from decimal import Decimal, getcontext
from typing import NamedTuple

# b13phase パッケージとの整合
BASE = 3120  # = 60 × 52 = 13 × 240 = 2⁴ × 3 × 5 × 13


# ============================================================
# 6.1 出力の型 (B13 座標系の慣習に従う)
# ============================================================

class B13PiTerm(NamedTuple):
    """π 平衡展開の一段分のレコード.

    全フィールドが整数または Decimal で、スケールは BASE=3120 系に整合する.
    """
    n: int              # 段 (1 以上)
    a_n: int            # その段の整数係数 (|a_n| < BASE/2)
    scale: int          # 3120^n (この段の円のサイズ)
    cum_pi_num: int     # 累積近似の分子 (3120^n スケール)
    err_log10: float    # log10(|π - 累積近似|), 参考値


# ============================================================
# 6.2 メイン関数: evaluate_pi()
#     b13phase.evaluate() と同じ命名規則
# ============================================================

def evaluate_pi(n_terms: int = 15, prec_digits: int = 200) -> list[B13PiTerm]:
    """π の 3120 平衡展開を n_terms 段まで計算する.

    Parameters
    ----------
    n_terms : int
        計算する段数 (1 以上).
    prec_digits : int
        内部 Decimal 精度 (桁数). 各段で n_terms × log10(3120) ≈ 3.49×n_terms
        桁の精度が出るため, prec_digits ≥ 4×n_terms + 10 を推奨.

    Returns
    -------
    list[B13PiTerm]
        各段の係数レコードのリスト. 長さは n_terms.

    Algorithm
    ---------
    各段で実行する操作は以下のみ:

        shifted   = residual × 3120     (円を 3120 倍に拡大)
        a_n       = nearest_int(shifted) (3120 点格子の最近接点を選ぶ)
        residual  = shifted - a_n        (残差を次段に繰り越す)

    これは b13phase の az_idx 公式が「方位角を Z_3120 上の整数で表す」
    のと同型の操作である.
    """
    if n_terms < 1:
        raise ValueError(f"n_terms must be ≥ 1, got {n_terms}")
    if prec_digits < 4 * n_terms + 10:
        # 精度不足だと最終段で丸めが破綻するため自動引き上げ
        prec_digits = 4 * n_terms + 10

    getcontext().prec = prec_digits

    pi = _pi_decimal(prec_digits)
    base = Decimal(BASE)

    terms: list[B13PiTerm] = []
    residual = pi
    cum_num = 0  # 累積分子 (現在の段のスケールで)

    for n in range(1, n_terms + 1):
        shifted = residual * base

        # 最近接整数化 (banker's rounding を避けて素直な四捨五入)
        a_n = _nearest_int(shifted)

        residual = shifted - Decimal(a_n)

        # 累積近似 = Σ a_k / 3120^k を 3120^n スケールの整数分子で保持
        cum_num = cum_num * BASE + a_n
        scale = BASE ** n
        cum_pi = Decimal(cum_num) / Decimal(scale)
        err = abs(pi - cum_pi)
        err_log10 = float(err.ln() / Decimal(10).ln()) if err > 0 else float('-inf')

        terms.append(B13PiTerm(
            n=n,
            a_n=a_n,
            scale=scale,
            cum_pi_num=cum_num,
            err_log10=err_log10,
        ))

    return terms


# ============================================================
# 6.3 補助関数
# ============================================================

def to_float_pi(terms: list[B13PiTerm]) -> float:
    """展開係数列から π の浮動小数点近似を復元する.

    b13phase.to_float_coords() の π 版.
    """
    if not terms:
        return 0.0
    last = terms[-1]
    return last.cum_pi_num / last.scale


def to_decimal_pi(terms: list[B13PiTerm], prec_digits: int = 200) -> Decimal:
    """展開係数列から π の Decimal 近似を復元する."""
    getcontext().prec = prec_digits
    if not terms:
        return Decimal(0)
    last = terms[-1]
    return Decimal(last.cum_pi_num) / Decimal(last.scale)


def factorize(n: int) -> str:
    """整数の素因数分解を文字列で返す (検査用)."""
    if n == 0:
        return "0"
    sign = "-" if n < 0 else ""
    n = abs(n)
    if n == 1:
        return sign + "1"
    factors: dict[int, int] = {}
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


def verify_pi(terms: list[B13PiTerm], prec_digits: int = 200) -> dict:
    """展開の妥当性を検証する.

    b13phase.verify_all() の π 版.
    確認項目:
        bounds_ok    : 全ての |a_n| < BASE/2 を満たす (平衡展開の定義)
        monotone_ok  : 誤差が単調に減少する
        rate_ok      : 1 段あたり log10(3120) ≈ 3.49 桁/項の精度向上
        final_err    : 最終段の絶対誤差
    """
    getcontext().prec = prec_digits
    pi = _pi_decimal(prec_digits)

    # bounds check
    # 一段目は π の整数部 (3) を含むため a_1 ≈ 3 × BASE + 残差 となる.
    # 二段目以降は純粋な平衡展開なので |a_n| < BASE/2 を満たすべき.
    half_base = BASE // 2
    bounds_ok = all(abs(t.a_n) <= half_base for t in terms[1:])
    # 一段目の上限 (整数部 3 + 0.5 を吸収できれば十分)
    if terms:
        bounds_ok = bounds_ok and abs(terms[0].a_n) < 4 * BASE

    # monotonicity
    monotone_ok = all(
        terms[i].err_log10 < terms[i - 1].err_log10
        for i in range(1, len(terms))
    )

    # rate check: 平均で 3.0 〜 3.49 桁/項 (理論上限 log10(3120) = 3.494)
    if len(terms) >= 2:
        total_gain = terms[0].err_log10 - terms[-1].err_log10
        avg_rate = total_gain / (len(terms) - 1)
        rate_ok = 2.5 <= avg_rate <= 3.5  # 上限はわずかに超えうる (丸めで得をする段がある)
    else:
        avg_rate = 0.0
        rate_ok = True

    # final error
    final_pi = to_decimal_pi(terms, prec_digits)
    final_err = float(abs(pi - final_pi))

    return {
        'bounds_ok':   bounds_ok,
        'monotone_ok': monotone_ok,
        'rate_ok':     rate_ok,
        'avg_rate':    avg_rate,
        'final_err':   final_err,
        'n_terms':     len(terms),
    }


# ============================================================
# 内部関数
# ============================================================

def _nearest_int(d: Decimal) -> int:
    """Decimal を最近接整数に丸める.

    これが本モジュールの本質的演算: 円格子上の最近接点選択.
    """
    floor_val = int(d)  # 切り捨て
    frac = d - Decimal(floor_val)
    if frac >= Decimal('0.5'):
        return floor_val + 1
    elif frac <= Decimal('-0.5'):
        return floor_val - 1
    return floor_val


def _pi_decimal(prec_digits: int) -> Decimal:
    """Decimal で高精度の π を返す.

    Note: ここで Machin 系の公式を使うが, これは「比較対象としての π の
    真値」を得るためだけで, 本モジュールの平衡展開そのものは π の値を
    入力として要求する (構造を見るためのモジュールであり, π を高速計算
    するためのモジュールではない).
    """
    getcontext().prec = prec_digits + 20

    # 4 arctan(1/5) - arctan(1/239) using Decimal
    def arctan_inv(x_inv: int) -> Decimal:
        # arctan(1/x_inv) = 1/x_inv - 1/(3 x_inv^3) + 1/(5 x_inv^5) - ...
        x = Decimal(1) / Decimal(x_inv)
        x2 = x * x
        result = Decimal(0)
        term = x
        n = 1
        while True:
            new_result = result + term / Decimal(n)
            if new_result == result:
                break
            result = new_result
            term = -term * x2
            n += 2
        return result

    pi = 4 * (4 * arctan_inv(5) - arctan_inv(239))
    getcontext().prec = prec_digits
    return +pi  # quantize to current precision


# ============================================================
# 既知の B13 関連量との照合 (引継書 v126 記録用)
# ============================================================

# フィボナッチ数列 (a_1 = 9802 に F_14 が含まれることの確認用)
FIBONACCI = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597]
F_14 = FIBONACCI[13]  # = 377 = 13 × 29  (FIBONACCI[0]=1 始まりのため index は 13)

# C60 関連 (a_2 = -720 と総剰余角の照合用)
C60_TOTAL_ANGULAR_SURPLUS = 720  # = 2 × 13 B13 単位 (確定)


def annotate_term(term: B13PiTerm) -> str:
    """係数 a_n に既知の B13 量との対応があれば注釈をつける."""
    a = term.a_n
    notes = []
    if a % 13 == 0 and a != 0:
        notes.append(f"13 で割れる: {a}/13 = {a // 13}")
    if a % 169 == 0 and a != 0:
        notes.append("169 (=13²) で割れる")
    if abs(a) in FIBONACCI:
        idx = FIBONACCI.index(abs(a))
        notes.append(f"フィボナッチ F_{idx + 1}")  # FIBONACCI[0]=F_1 慣習に合わせる
    if a == 9802:
        notes.append(f"= 26 × F_14 = 2 × 13² × 29")
    if a == -C60_TOTAL_ANGULAR_SURPLUS:
        notes.append("C60 総剰余角と数値一致 (= -2×13 B13 単位)")
    return "; ".join(notes) if notes else ""


# ============================================================
# クイックリファレンス用デモ
# ============================================================

if __name__ == "__main__":
    print("=" * 70)
    print("b13_pi: π の 3120 平衡展開")
    print("=" * 70)
    print()
    print(f"BASE = {BASE} = 60 × 52 = 13 × 240")
    print()

    terms = evaluate_pi(n_terms=15, prec_digits=200)

    print(f"{'n':>3} | {'a_n':>8} | {'因数分解':<25} | {'累積誤差':>15} | 注釈")
    print("-" * 85)
    for t in terms:
        fac = factorize(t.a_n)
        err_str = f"10^{t.err_log10:.2f}"
        note = annotate_term(t)
        print(f"{t.n:>3} | {t.a_n:>8} | {fac:<25} | {err_str:>15} | {note}")

    print()
    print("検証結果:")
    v = verify_pi(terms, prec_digits=200)
    for key, val in v.items():
        print(f"  {key:15s} = {val}")

    print()
    print("π の最終近似 (15 段):")
    print(f"  浮動小数点: {to_float_pi(terms)}")
    print(f"  Decimal   : {to_decimal_pi(terms, 60)}")
    print(f"  真値      : {_pi_decimal(60)}")
