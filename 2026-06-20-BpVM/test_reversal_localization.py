"""
test_reversal_localization.py
B13 / reversal_lane 切分テスト #2 （新規・engine 非依存・玩具）

#1 (test_reversal_isolation.py) の presence-readout は不合格として記録済み。
ここでは fixture（信号生成・位相シャッフル・residue 抽出）を #1 から import で
共有し、readout だけを order-localization 版に差し替える。

問い: presence でなく order-localization を読めば、surrogate は壊れるか。
固定した pass 条件（antenna 指定・出た数字をそのまま読む）:
  real      : ラン二本・全長被覆 1.0・境界一点  -> score 1.0
  surrogate : ラン三本以上・被覆低下・境界散在  -> score 沈む
  null 天井 : 同じ surrogate アンサンブルから取る
  steady    : 0/0 のまま
real が立ち surrogate 99pct/max が沈めば「機構候補」として前進。
max がまだ real に届けば reversal も §4 開け直し。通るとは declare しない。
"""
import numpy as np
from test_reversal_isolation import (
    P, make_signal, phase_shuffle, extract_residue, runs,
)

# ---------- order-localization readout（新規・固定） ----------
def stable_runs(walk, min_run=3):
    out, pos = [], 0
    for (v, l) in runs(walk):
        if l >= min_run:
            out.append((v, l, pos))
        pos += l
    return out

def localization_score(walk, min_run=3):
    """
    antenna 指定の三判別子をそのまま:
      coverage      = 最良の隣接対蹠安定ラン対が覆う割合
      run_factor    = 2 / max(2, 安定ラン総数)   （丁度二本=1.0、多いほど散在penalty）
    score = coverage * run_factor。対蹠でなければ 0、s=0 は除外。
    """
    segs = stable_runs(walk, min_run)
    n = len(segs)
    best_cov = 0.0
    for (va, la, _), (vb, lb, _) in zip(segs, segs[1:]):
        if va != 0 and vb == (-va) % P:
            best_cov = max(best_cov, (la + lb) / len(walk))
    run_factor = 2.0 / max(2, n)
    return best_cov * run_factor, n, best_cov

def describe(walk, min_run=3):
    segs = stable_runs(walk, min_run)
    return f"stable_runs={[(v,l) for v,l,_ in segs]}"

# ---------- 較正 + 判定 ----------
def calibrate_and_test(s, n_surr=400, seed=0):
    rng = np.random.default_rng(seed)
    x = make_signal(s, reverse=True)
    real_walk = extract_residue(x)
    real_sc, real_n, real_cov = localization_score(real_walk)

    surr_sc, worst = [], (-1, None)
    for _ in range(n_surr):
        w = extract_residue(phase_shuffle(x, rng))
        sc, _, _ = localization_score(w)
        surr_sc.append(sc)
        if sc > worst[0]:
            worst = (sc, w)
    surr_sc = np.array(surr_sc)
    ceiling = np.percentile(surr_sc, 99)
    smax = surr_sc.max()

    print(f"--- reversal localization (s={s}) ---")
    print(f"  real : score={real_sc:.4f}  runs={real_n}  coverage={real_cov:.3f}  ({describe(real_walk)})")
    print(f"  surr : mean={surr_sc.mean():.4f}  99pct={ceiling:.4f}  max={smax:.4f}  frac>0={(surr_sc>0).mean():.3f}")
    print(f"  worst surrogate: {describe(worst[1])}")
    verdict = ("前進(real>ceiling & max沈む)" if real_sc > ceiling and smax < real_sc
               else "max が real に届く -> §4 開け直し" if smax >= real_sc
               else "real>ceiling だが max は real 未満")
    print(f"  separation={real_sc/(ceiling+1e-12):.2f}x   verdict = {verdict}")
    print()
    return real_sc, ceiling, smax

def negative_control(s=3, n_surr=400, seed=1):
    rng = np.random.default_rng(seed)
    x = make_signal(s, reverse=False)
    rw = extract_residue(x)
    real_sc, *_ = localization_score(rw)
    surr = np.array([localization_score(extract_residue(phase_shuffle(x, rng)))[0] for _ in range(n_surr)])
    print(f"--- negative control: steady (s={s}, 反転なし) ---")
    print(f"  real score={real_sc:.4f}  surr mean={surr.mean():.4f}   期待 0/0")
    print()

if __name__ == "__main__":
    print("="*64)
    print("B13 reversal_lane 切分テスト #2 : order-localization readout")
    print("="*64)
    for s in (1, 2, 3):
        calibrate_and_test(s)
    negative_control()
