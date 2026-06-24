"""
readout_core.py  —  B13 v158 #4 directional readout (核ロジック)

v157 から逆算して確定した表現と #2 スコアを再現し、
その上に #4 directional readout を載せる。

run 表現:  [(residue, length), ...]  （min_run=3 で短ランを落とす）
#2 score:  coverage(最良の対蹠隣接ペア被覆率) × run_factor(2/max(2,#runs))   ← 無向き
#4 score:  同じ被覆だが、s-run が先・(-s)-run が後 の順序を要求           ← 方向忠実
"""

MOD = 13


def anti(s):
    """対蹠 residue: -s mod 13"""
    return (-s) % MOD


def total_len(runs):
    return sum(l for _, l in runs)


def num_runs(runs):
    return len(runs)


# ---- #2 無向き（v157 §3 をそのまま再現） ----
def coverage_undirected(runs, s, n_total=None):
    """最良の『対蹠で隣接』run ペアの被覆率。順序は問わない。"""
    if n_total is None:
        n_total = total_len(runs)
    if n_total == 0:
        return 0.0, None
    a, b = s % MOD, anti(s)
    best, best_pair = 0.0, None
    for i in range(len(runs) - 1):
        (r0, l0), (r1, l1) = runs[i], runs[i + 1]
        pair = {r0, r1}
        if pair == {a, b}:                      # 対蹠ペアが隣接（向き不問）
            cov = (l0 + l1) / n_total
            if cov > best:
                best, best_pair = cov, (runs[i], runs[i + 1])
    return best, best_pair


def run_factor(runs):
    return 2.0 / max(2, num_runs(runs))


def score_undirected(runs, s, n_total=None):
    cov, _ = coverage_undirected(runs, s, n_total)
    return cov * run_factor(runs)


# ---- #4 方向忠実（v158 の修復） ----
def coverage_directional(runs, s, gap_tol=1, n_total=None):
    """
    s-run が先、その直後または許容ギャップ gap_tol 内に (-s)-run が続く場合のみ credit。
    逆順 (-s)→s は credit しない。
    """
    if n_total is None:
        n_total = total_len(runs)
    if n_total == 0:
        return 0.0, None
    a, b = s % MOD, anti(s)
    best, best_pair = 0.0, None
    for i in range(len(runs)):
        if runs[i][0] != a:                     # 先頭は必ず s
            continue
        # i の直後〜gap_tol 内に -s run を探す
        for j in range(i + 1, min(i + 1 + gap_tol + 1, len(runs))):
            if runs[j][0] == b:                 # 後続が -s
                cov = (runs[i][1] + runs[j][1]) / n_total
                if cov > best:
                    best, best_pair = cov, (runs[i], runs[j])
                break                            # 最近接の -s で確定
    return best, best_pair


def score_directional(runs, s, gap_tol=1, n_total=None):
    cov, _ = coverage_directional(runs, s, gap_tol, n_total)
    return cov * run_factor(runs)


# =========================================================
#  Anchor unit-proof : v157 が記録した最悪 surrogate と real
# =========================================================
if __name__ == "__main__":
    print("=" * 64)
    print(" #4 directional readout — anchor unit-proof（fixture 非依存・確定）")
    print("=" * 64)

    # v157 §3: real s=1 の run は [(1,12),(12,12)]、cov1.0・runs2
    real_s1 = [(1, 12), (12, 12)]
    # v157 §3: s=1 の最悪 surrogate（薄合格の主因）
    worst_surr_s1 = [(12, 10), (1, 8)]          # n_total=18 → cov 0.75

    # n_total=24（全フレーム）基準。短ランは min_run=3 で落ち、被覆の分母は常に 24。
    cases = [
        ("real  s=1", real_s1, 1, 24),
        ("surr  s=1 (worst)", worst_surr_s1, 1, 24),       # 18/24 = 0.75
        ("surr  as s=12 trial", worst_surr_s1, 12, 24),    # -12 mod13 = 1 → 12→1 は正順
    ]

    print(f"\n{'case':22s} {'#2 undir cov':>13s} {'#2 score':>9s}"
          f" {'#4 dir cov':>11s} {'#4 score':>9s}")
    print("-" * 70)
    for name, runs, s, nt in cases:
        cu, pu = coverage_undirected(runs, s, nt)
        cd, pd = coverage_directional(runs, s, gap_tol=1, n_total=nt)
        su = score_undirected(runs, s, nt)
        sd = score_directional(runs, s, gap_tol=1, n_total=nt)
        print(f"{name:22s} {cu:13.3f} {su:9.3f} {cd:11.3f} {sd:9.3f}")

    print("\n判定:")
    print(" - real s=1:    無向き 1.0 / 方向 1.0   → 正順 1→12 を credit（維持）")
    print(" - 最悪surr s=1: 無向き 0.75 / 方向 0.0  → 逆順 12→1 を棄却（#2 の薄さの主因が消える）")
    print(" - 同 surr を s=12: 方向 0.75            → 12→1 は s=12 の正順、合格候補（定義整合）")
