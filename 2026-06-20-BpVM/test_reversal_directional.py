"""
test_reversal_directional.py
B13 / reversal_lane 切分テスト #4 （新規・engine 非依存・玩具）

#2 (test_reversal_localization.py) の localization_score は、対蹠ペアの
「隣接」だけを読み、向きを読まない（無向き antipodal adjacency）。
ロック定義は「stable run s に *続いて* stable run -s mod 13」＝向きを含む。
ここでは fixture（#1）と stable_runs（#2）を import で共有し、
readout の対蹠判定だけを方向忠実へ差し替える。動かすのは readout 一本だけ。

固定（触らない）: 信号生成・位相シャッフル・residue 抽出・stable_runs・null 天井(99pct)・n_surr。
差替え: 対蹠判定  vb==(-va)%P  →  (va==s) かつ (vb==(-s)%P)   ［s 先・-s 後のみ credit］
緑にするための調整はしない。surrogate・閾値は #2 と同一。出た数字をそのまま読む。
"""
import numpy as np
from test_reversal_isolation import P, make_signal, phase_shuffle, extract_residue
from test_reversal_localization import stable_runs, localization_score


def directional_localization_score(walk, s, min_run=3):
    """#2 と同一構造。違いは対蹠判定のみ: va==s%P かつ vb==(-s)%P（s 先・-s 後のみ credit）。"""
    segs = stable_runs(walk, min_run)
    n = len(segs)
    a, b = s % P, (-s) % P
    best_cov = 0.0
    for (va, la, _), (vb, lb, _) in zip(segs, segs[1:]):
        if va == a and vb == b:                 # 方向忠実
            best_cov = max(best_cov, (la + lb) / len(walk))
    run_factor = 2.0 / max(2, n)
    return best_cov * run_factor, n, best_cov


def describe(walk, min_run=3):
    segs = stable_runs(walk, min_run)
    return f"stable_runs={[(v,l) for v,l,_ in segs]}"


def calibrate(s, n_surr=400, seed=0):
    rng = np.random.default_rng(seed)
    x = make_signal(s, reverse=True)
    real_walk = extract_residue(x)
    real_u, _, _ = localization_score(real_walk)
    real_d, _, _ = directional_localization_score(real_walk, s)
    su, sd = [], []
    worst_u = (-1.0, None); worst_d = (-1.0, None)
    for _ in range(n_surr):
        w = extract_residue(phase_shuffle(x, rng))
        scu, _, _ = localization_score(w)
        scd, _, _ = directional_localization_score(w, s)
        su.append(scu); sd.append(scd)
        if scu > worst_u[0]: worst_u = (scu, w)
        if scd > worst_d[0]: worst_d = (scd, w)
    su, sd = np.array(su), np.array(sd)
    return dict(
        s=s, real_u=real_u, real_d=real_d,
        u_99=np.percentile(su, 99), u_max=su.max(), u_fracpos=(su > 0).mean(),
        d_99=np.percentile(sd, 99), d_max=sd.max(), d_fracpos=(sd > 0).mean(),
        worst_u=worst_u[1], worst_d=worst_d[1],
    )


if __name__ == "__main__":
    print("=" * 72)
    print("B13 reversal_lane 切分テスト #4 : directional readout（向き忠実）")
    print("=" * 72)
    print("固定: 信号・位相シャッフル・抽出・stable_runs・null 天井(99pct)・n_surr。")
    print("差替え: 対蹠判定のみ（s 先・-s 後）。surrogate は #2 と同一アンサンブル。\n")
    rows = []
    for s in (1, 2, 3):
        r = calibrate(s); rows.append(r)
        sep_u = r['real_u'] / (r['u_99'] + 1e-12)
        sep_d = r['real_d'] / (r['d_99'] + 1e-12)
        print(f"--- s={s} ---")
        print(f"  real        : undirected={r['real_u']:.4f}   directional={r['real_d']:.4f}")
        print(f"  surr 99pct  : undirected={r['u_99']:.4f}   directional={r['d_99']:.4f}")
        print(f"  surr max    : undirected={r['u_max']:.4f}   directional={r['d_max']:.4f}")
        print(f"  surr frac>0 : undirected={r['u_fracpos']:.3f}    directional={r['d_fracpos']:.3f}")
        print(f"  separation  : undirected={sep_u:.2f}x   directional={sep_d:.2f}x")
        print(f"  margin(99)  : undirected={r['real_u']-r['u_99']:+.4f}   directional={r['real_d']-r['d_99']:+.4f}")
        print(f"  margin(max) : undirected={r['real_u']-r['u_max']:+.4f}   directional={r['real_d']-r['d_max']:+.4f}")
        print(f"  worst surr(undir)={describe(r['worst_u'])}")
        print(f"  worst surr(dir) ={describe(r['worst_d'])}")
        print()
    print("=" * 72)
    print("総括（s=1＝最薄を読む）:")
    r1 = rows[0]
    print(f"  #2 無向き: surr max {r1['u_max']:.3f} / 99pct {r1['u_99']:.3f}  -> margin(99) {r1['real_u']-r1['u_99']:+.3f}")
    print(f"  #4 方向  : surr max {r1['d_max']:.3f} / 99pct {r1['d_99']:.3f}  -> margin(99) {r1['real_d']-r1['d_99']:+.3f}")
    print(f"  real は方向化で {r1['real_u']:.3f} -> {r1['real_d']:.3f}")
