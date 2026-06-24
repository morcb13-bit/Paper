"""
run_engine_b13.py
engine_b13 を v159 harness に通し、自己検出版の real_d を測る。
緑化しない。検出が外れる所・real_d が落ちる所をそのまま出す。
"""
import numpy as np
from test_reversal_isolation import make_signal, phase_shuffle, extract_residue
from test_reversal_directional import directional_localization_score
from engine_reversal_eval import evaluate, baseline_gap99_at
from engine_b13 import engine_extract, detect_s

P = 13

print("=" * 78)
print("engine_b13 : 自己検出エンジンを v159 harness で測る")
print("=" * 78)

# ---- A. clean・全 s で検出正否と real_d ----
print("\n[A] clean 信号・全 s（s_true vs s_detected, real_d, gap99）")
print(f"  {'s_true':>6} {'s_det':>5} {'hit':>4} {'real_d':>7} {'surr99':>7} {'gap99':>7}  verdict")
n_hit = 0
for s_true in range(1, P):
    x = make_signal(s_true, reverse=True)
    r = evaluate(x, engine_extract, n_surr=300, seed=7)
    hit = (r['detected_s'] == s_true)
    n_hit += hit
    print(f"  {s_true:>6} {r['detected_s']:>5} {'✓' if hit else '✗':>4} "
          f"{r['real_d']:>7.3f} {r['surr_99']:>7.3f} {r['gap99']:>+7.3f}  {r['verdict'][:34]}")
print(f"  --> clean 検出正答 {n_hit}/{P-1}")

# ---- B. negative control: 反転無し steady ----
print("\n[B] negative control（反転なし steady, s=3）— 偽発火しないか")
xneg = make_signal(3, reverse=False)
walk_neg = extract_residue(xneg)
s_neg = detect_s(walk_neg)
rneg, _, cov = directional_localization_score(walk_neg, s_neg)
print(f"  walk      = {list(map(int, walk_neg))}")
print(f"  detect_s  = {s_neg}  real_d = {rneg:.3f}  cov = {cov:.3f}  （期待: real_d~0）")

# ---- C. noise 掃引・s=1（最薄）: 検出正答率 + real_d vs baseline ----
print("\n[C] noise 掃引 s=1（最薄）: 検出正答率と real_d の baseline 追従")
print(f"  {'sigma':>6} {'SNR(dB)':>8} {'det_hit%':>8} {'real_d_mean':>11} "
      f"{'real_d_min':>10} {'base_gap99':>10}  status")
S = 1
N_NOISE = 20
x0 = make_signal(S, reverse=True)
sig_rms = np.sqrt(np.mean(x0**2))
for sigma in [0.0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.00, 1.50]:
    rng = np.random.default_rng(2026)
    snr = 20*np.log10(sig_rms/sigma) if sigma > 0 else np.inf
    reals, hits = [], 0
    for _ in range(N_NOISE):
        xn = x0 + rng.normal(0, sigma, size=len(x0))
        walk = extract_residue(xn).astype(int)
        s_det = detect_s(walk)
        hits += (s_det == S)
        rd, _, _ = directional_localization_score(walk, s_det)  # 検出 s で採点
        reals.append(rd)
    reals = np.array(reals)
    base = baseline_gap99_at(snr)
    snrs = f"{snr:.1f}" if np.isfinite(snr) else "inf"
    status = ("real~1.0" if reals.mean() >= 0.9
              else "薄/崩れ" if reals.mean() >= 0.3 else "盲目")
    print(f"  {sigma:>6.2f} {snrs:>8} {100*hits/N_NOISE:>7.0f}% "
          f"{reals.mean():>11.3f} {reals.min():>10.3f} {base:>10.3f}  {status}")

print("\n" + "=" * 78)
print("読み: 検出が外れた所で real_d は #4 が間違った対蹠を探して落ちる（条件6＝盲目）。")
print("玩具(s ハードコード)はこの故障様式を隠していた。自己検出で初めて出る。")
print("=" * 78)
