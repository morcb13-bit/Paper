"""
run_compete_b13.py
B13 #B : detect_s 取り違えを *わざと* 踏ませ、harness が
  real-collapse（抽出劣化・det 正答のまま real_d 落ち）と
  detect-miss（誤 s ロック・det 失敗）
を本当に分離できるかを検める。緑化しない。出た数字をそのまま読む。

競合信号: 真の反転(s_true) に、別 s_decoy で反転する第二成分を amp で足す。
  sep=0  : 同一ビン → 位相が混ざり walk がどちらでもなくなる（混濁崩壊）
  sep>=2 : 別ビン   → argmax が切替わり walk が丸ごと decoy へ（切替崩壊）
"""
import numpy as np
from test_reversal_isolation import P, FRAME, HOP, N_EACH, DOM, extract_residue
from test_reversal_directional import directional_localization_score
from engine_b13 import detect_s, engine_extract
from engine_reversal_eval import evaluate
from engine_admissibility import admit

def _chirp(s, fbin, reverse=True):
    n_frames = 2*N_EACH
    N = n_frames*HOP + FRAME
    off = (s/P)/HOP
    half = N//2
    fc = np.empty(N)
    fc[:half] = fbin + off
    fc[half:] = fbin + (-off if reverse else off)
    return np.cos(2*np.pi*np.cumsum(fc))

def make_competing(s_true, s_decoy, amp, sep=0, reverse=True):
    """真 s_true(振幅1) + decoy s_decoy(振幅 amp, ビン +sep) の和。"""
    a = _chirp(s_true, DOM/FRAME, reverse)
    b = _chirp(s_decoy, (DOM+sep)/FRAME, reverse)
    return a + amp*b

def classify(hit, real_d):
    if hit and real_d >= 0.9:   return "clean"
    if hit and real_d < 0.9:    return "real-collapse"
    if (not hit) and real_d >= 0.9: return "DETECT-MISS(silent: real_d 健全のまま誤s)"
    return "detect-miss-collapse"

S_TRUE, S_DECOY = 1, 5     # 真=1(→12), decoy=5(→8)。両方とも s!=0 の対蹠対
print("="*86)
print(f"competing antipode: s_true={S_TRUE}(→{(-S_TRUE)%P})  s_decoy={S_DECOY}(→{(-S_DECOY)%P})")
print("harness 二列 (det_hit, real_d) で崩壊原因を名指せるか")
print("="*86)

for sep, tag in [(0, "sep=0  同一ビン(混濁)"), (3, "sep=3  別ビン(切替)")]:
    print(f"\n--- {tag} ---")
    print(f"  {'amp':>5} {'det_s':>5} {'hit':>4} {'real_d':>7} {'cov':>6}  mode")
    seen = set()
    for amp in [0.0, 0.3, 0.6, 0.9, 1.0, 1.1, 1.3, 1.6, 2.0]:
        x = make_competing(S_TRUE, S_DECOY, amp, sep=sep)
        walk = extract_residue(x).astype(int)
        s_det = detect_s(walk)
        hit = (s_det == S_TRUE)
        rd, _, cov = directional_localization_score(walk, s_det)
        mode = classify(hit, rd)
        seen.add(mode.split("(")[0])
        print(f"  {amp:>5.2f} {s_det:>5} {'✓' if hit else '✗':>4} {rd:>7.3f} {cov:>6.3f}  {mode}")
    print(f"  到達した mode: {sorted(seen)}")

# ---- 4象限がすべて埋まるか（分離能の証明）----
print("\n" + "="*86)
print("分離能チェック: (det_hit × real_d) の4象限がサーベイ内で全部出るか")
print("="*86)
quad = {}
for sep in (0, 1, 2, 3):
    for amp in np.linspace(0, 2.2, 23):
        x = make_competing(S_TRUE, S_DECOY, float(amp), sep=sep)
        walk = extract_residue(x).astype(int)
        s_det = detect_s(walk)
        hit = (s_det == S_TRUE)
        rd, _, _ = directional_localization_score(walk, s_det)
        key = ("det_hit" if hit else "det_miss", "real>=0.9" if rd >= 0.9 else "real<0.9")
        quad.setdefault(key, []).append((sep, round(float(amp),2), s_det, round(rd,3)))
for k in [("det_hit","real>=0.9"),("det_hit","real<0.9"),
          ("det_miss","real>=0.9"),("det_miss","real<0.9")]:
    ex = quad.get(k)
    label = {("det_hit","real>=0.9"):"clean",
             ("det_hit","real<0.9"):"real-collapse",
             ("det_miss","real>=0.9"):"DETECT-MISS(silent)",
             ("det_miss","real<0.9"):"detect-miss-collapse"}[k]
    if ex:
        sep,amp,sd,rd = ex[0]
        print(f"  [{'埋' if ex else '空'}] {k[0]:>8} × {k[1]:>9} = {label:<26} "
              f"例 sep={sep} amp={amp} det_s={sd} real_d={rd}  (n={len(ex)})")
    else:
        print(f"  [空] {k[0]:>8} × {k[1]:>9} = {label:<26} ← この象限に到達せず")

# ---- 全 harness 経由で admissibility が通ることも確認（silent ケースを一本）----
print("\n" + "="*86)
print("harness evaluate() 経由の確認（admissibility 通過 + verdict 文字列）")
print("="*86)
for sep, amp in [(3, 1.6), (0, 1.0)]:
    x = make_competing(S_TRUE, S_DECOY, amp, sep=sep)
    try:
        r = evaluate(x, engine_extract, n_surr=200, seed=11)
        hit = (r['detected_s'] == S_TRUE)
        print(f"  sep={sep} amp={amp}: detected_s={r['detected_s']} (hit={hit}) "
              f"real_d={r['real_d']:.3f} gap99={r['gap99']:+.3f}")
        print(f"     verdict='{r['verdict']}'")
        print(f"     → harness verdict は s 正否を見ない。det_hit を併読しないと "
              f"{'silent miss を見逃す' if (not hit and r['real_d']>=0.9) else '一致'}")
    except Exception as e:
        print(f"  sep={sep} amp={amp}: REJECT/err {e}")
