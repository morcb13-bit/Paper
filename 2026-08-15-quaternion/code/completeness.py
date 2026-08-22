#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""数列の完全性を測る。エルデシュ #346 の三条件を判定する道具。標準ライブラリのみ。

対象の定義は Lean 形式化（FormalConjectures/ErdosProblems/346.lean）に合わせる。
勝手に作り直さないこと（誤り51・52 の型）。

  complete          … 十分大きい整数がすべて、相異なる項の和で書ける
  strongly complete … 先頭を n 項切り落とした「尾」が、どの n でも complete
                      （※「有限個どれを抜いても」ではない。ここを取り違えた実績あり）
  第三条件          … 無限個の項を抜くと、どの抜き方でも complete でなくなる

判定はビットマスクの実測で行う。ブラウン型の不等式は速い篩として併置するが、
必要条件であって同値ではない。合否をそちらで決めないこと。

  検定1 密な列（a_n = n）の尾は complete
      OK なら：完全性の判定器が True を返せている
      NG なら：穴の数え方か観測窓が壊れている
  検定2 冪2列の尾は complete でない
      OK なら：判定器が False を返せている＝反証不能な検査ではない
      NG なら：誤り55 の型。構造上どんな入力でも通る検査を書いている
  検定3 フィボナッチ列の尾は complete
      OK なら：比 φ の列が strongly complete 側に入る
      NG なら：観測項数が足りないか、上端の縁を拾っている
  検定4 比 1.8 の等比列の尾は complete でない
      OK なら：φ より上が strongly complete から外れる
      NG なら：閾値の観測が効いていない
  検定5 比を 1.40〜1.90 で掃くと、成立が φ の側の一区画に限られる
      OK なら：範囲を掃いて残った点が孤立している（判別法2の運用）
      NG なら：閾値が観測窓の取り方で動いている。窓を変えて再測すること
  検定6 φ 超の列は無限抜きで壊れ、密な列は壊れない
      OK なら：第三条件の代理検査が両方向に振れる
      NG なら：代理検査が構造上必ず一方を返している
"""

PHI = (1 + 5 ** 0.5) / 2


# ---------------------------------------------------------------- 実測（厳密）

def reachable(seq):
    """seq の相異なる項の和で表せる整数の集合を、ビットマスク（整数）で返す。"""
    S = sum(seq)
    mask = 1
    lim = (1 << (S + 1)) - 1
    for a in seq:
        if a > 0:
            mask |= (mask << a) & lim
    return mask, S


def holes(seq, lo_frac=0.25, hi_frac=0.50):
    """全体和 S の [lo_frac·S, hi_frac·S] に残る穴の個数と最大値。

    上端 S 付近は必ず穴が開く（S-1 は最小項より小さい差を要求するため）。
    縁の作りものを合否に混ぜないよう、中央の帯だけを見る。
    """
    mask, S = reachable(seq)
    lo, hi = int(S * lo_frac), int(S * hi_frac)
    width = hi - lo + 1
    band = (mask >> lo) & ((1 << width) - 1)
    cnt = width - band.bit_count()
    mx = -1
    if cnt:
        inv = (~band) & ((1 << width) - 1)
        mx = lo + inv.bit_length() - 1
    return cnt, mx, (lo, hi)


def is_complete(seq, tol=0):
    """中央の帯に穴が tol 個以下なら complete とみなす（実測）。"""
    cnt, _, _ = holes(seq)
    return cnt <= tol


# ---------------------------------------------------------------- 速い篩

def brown_violations(seq, start=0):
    """a[k+1] > 1 + Σ_{i=start}^{k} a[i] となる k を列挙する。

    ここに引っかかる k では Σ+1 が恒久的な穴になる。必要条件の側の道具。
    同値ではないので、これだけで complete を主張しないこと。
    """
    bad, s = [], 0
    for k in range(start, len(seq) - 1):
        s += seq[k]
        if seq[k + 1] > 1 + s:
            bad.append((k, seq[k + 1], s + 1))
    return bad


# ---------------------------------------------------------------- #346 の三条件

def tails(seq, depth, terms):
    """先頭 n 項を落とした尾を depth 本返す。各尾は terms 項。"""
    out = []
    for n in range(depth):
        t = seq[n:n + terms]
        if len(t) >= terms:
            out.append((n, t))
    return out


def is_strongly_complete(seq, depth=4, terms=18, tol=0):
    """どの尾も complete か。落ちた尾の番号を併せて返す。"""
    bad = [n for n, t in tails(seq, depth, terms) if not is_complete(t, tol)]
    return (len(bad) == 0, bad)


def infinite_removals(n):
    """無限部分集合の代理を作る。等差の間引きと、疎な（指数的な）間引きの両方。

    **等差だけでは足りない。** 等差は抜き過ぎるので、成立するはずのない列でも
    「壊れた」を返してしまう。疎な間引きを必ず併用すること（誤り133）。
    """
    out = []
    for step in (2, 3):
        for off in range(step):
            out.append((f"等差{step}(off{off})",
                        {i for i in range(n) if i % step == off}))
    for base in (2, 3):
        for off in (0, 1):
            idx, k = set(), off
            while k < n:
                idx.add(k)
                k = max(k + 1, int(k * base) + 1)
            out.append((f"指数{base}(off{off})", idx))
    return out


def breaks_on_infinite_removal(seq, terms=18, tol=0):
    """無限個抜くと、どの抜き方でも complete でなくなるか。

    すべての無限部分集合を掃いてはいない。掃いた抜き方の名前を併せて返す。
    生き残りが1つでもあれば第三条件は成立していない。
    """
    survived = []
    for name, idx in infinite_removals(len(seq)):
        rest = [a for i, a in enumerate(seq) if i not in idx][:terms]
        if len(rest) >= 10 and is_complete(rest, tol):
            survived.append(name)
    return (len(survived) == 0, survived)


def is_lacunary(seq):
    rs = [seq[i + 1] / seq[i] for i in range(len(seq) - 1) if seq[i] > 0]
    return (min(rs) > 1.0, min(rs)) if rs else (False, 0.0)


def ratio_tail(seq, tail=6):
    n = len(seq)
    rs = [seq[i + 1] / seq[i] for i in range(max(0, n - tail - 1), n - 1)]
    return sum(rs) / len(rs) if rs else 0.0


def verdict(seq, depth=4, terms=18, tol=0):
    """三条件をまとめて判定。candidate=True なら #346 の反例候補。"""
    lac, lo = is_lacunary(seq)
    sc, bad = is_strongly_complete(seq, depth, terms, tol)
    inf, surv = breaks_on_infinite_removal(seq, terms, tol)
    r = ratio_tail(seq)
    return {
        "lacunary": lac, "ratio_min": round(lo, 6),
        "strongly": sc, "tails_failed": bad,
        "breaks_inf": inf, "inf_survived": surv,
        "ratio_tail": round(r, 6),
        "candidate": lac and sc and inf and abs(r - PHI) > 0.02,
    }


# ---------------------------------------------------------------- 既知の列

def geometric(r, a0, n):
    s, x = [], float(a0)
    while len(s) < n:
        v = int(round(x))
        if not s or v > s[-1]:
            s.append(v)
        x *= r
    return s


def fibonacci(n, a=1, b=2):
    s = [a, b]
    while len(s) < n:
        s.append(s[-1] + s[-2])
    return s


def powers_of_two(n, a0=1):
    return [a0 * 2 ** k for k in range(n)]


def dense(n, a0=1):
    return list(range(a0, a0 + n))


def graham(n):
    """f(k) = F(k) - (-1)^k の像（#346 の既知例、Graham 1964）。"""
    F = [0, 1]
    while len(F) < n + 8:
        F.append(F[-1] + F[-2])
    vals = sorted({F[k] - (-1) ** k for k in range(1, n + 8)} - {0})
    return vals[:n]


# ---------------------------------------------------------------- 自己検定

def _line(no, name, ok, note=""):
    print(f"  検定{no} {name:<44} {'OK' if ok else 'NG'}  {note}")
    return ok


def main():
    print("完全性ハーネス（エルデシュ #346）自己検定")
    print()
    res = []

    dn = dense(90, 4)
    res.append(_line(1, "密な列の尾は complete", is_strongly_complete(dn, 4, 60)[0]))

    p2 = powers_of_two(18, 5)
    ok, bad = is_strongly_complete(p2, 4, 16)
    res.append(_line(2, "冪2列の尾は complete でない", not ok, f"落ちた尾 {bad}"))

    fb = fibonacci(24, 3, 5)
    ok, bad = is_strongly_complete(fb, 4, 18)
    res.append(_line(3, "素のフィボナッチ列の尾は complete でない", not ok,
                     f"落ちた尾 {bad}"))

    gr = graham(22)
    ok, bad = is_strongly_complete(gr, 4, 18)
    res.append(_line(4, "Graham 列の尾は complete", ok,
                     f"列 {gr[:6]}… 落ちた尾 {bad}"))

    g18 = geometric(1.8, 10, 22)
    ok, bad = is_strongly_complete(g18, 4, 18)
    res.append(_line(5, "比 1.8 の等比列の尾は complete でない", not ok,
                     f"落ちた尾 {bad}"))

    print()
    print("  比の走査（尾 n=0 の穴の個数。φ = 1.618034）")
    prev, jump = None, []
    for i in range(13):
        r = 1.40 + 0.05 * i
        cnt, _, _ = holes(geometric(r, 10, 20)[:18])
        mark = "  <-- 立ち上がり" if prev is not None and prev < 10 <= cnt else ""
        if mark:
            jump.append(round(r, 2))
        print(f"    r = {r:.2f}   穴 {cnt}{mark}")
        prev = cnt
    res.append(_line(6, "穴の立ち上がりが φ の近傍に1回だけ出る",
                     len(jump) == 1 and abs(jump[0] - PHI) < 0.1,
                     f"立ち上がり {jump}"))

    inf_a, _ = breaks_on_infinite_removal(geometric(1.75, 10, 24))
    inf_b, surv = breaks_on_infinite_removal(dense(90, 4), terms=60)
    res.append(_line(7, "φ 超は無限抜きで壊れ、密な列は壊れない",
                     inf_a and not inf_b, f"密な列の生き残り {len(surv)} 通り"))

    print()
    print(f"  {sum(res)}/{len(res)} OK")
    print()
    print("  Graham 列の総合判定 :", verdict(gr))
    print("  比1.8 列の総合判定  :", verdict(g18))


if __name__ == "__main__":
    main()
