"""
bundle_b13.py
B-3 束定義 v0（director 固定・そのまま実装・広げない）。

prior（明示・仮定）:
  真の反転は独立観測をまたいで不変に残り、decoy は nuisance 軸を振ることで不変性を失う。
  B-3 はこの prior の *検査* であり、真を単一信号から生成する手続きではない。

束定義 v0:
  M = 3（最小の独立反復）
  nuisance 軸（束内で振る）: decoy_s, decoy_sep(bin), decoy_phase, decoy_amp, noise_seed
  固定          : s_true, true_bin, true_amp, true_phase law, extractor
  安定閾値      : unanimous（有効票すべて一致・1本でも割れたら BUNDLE_ABSTAIN）
  有効票        : 単体 verdict==FORWARD のものだけ。有効票<M なら BUNDLE_ABSTAIN。全本 BLIND なら BUNDLE_BLIND。
  合格意味      : BUNDLE_FORWARD = cross-signal stable dominant s（≠ 真・束内不変）
"""
import numpy as np
from test_reversal_isolation import P, FRAME, HOP, N_EACH, DOM
from engine_reversal_eval import evaluate
from engine_b13 import engine_extract
from guard_b13 import abstain_check, verdict3

M = 3

# ---- 固定 true_bin=DOM, true_amp=1, true_phase law=0 ----
def _chirp(s, fbin, reverse=True, phase=0.0):
    n_frames = 2 * N_EACH
    N = n_frames * HOP + FRAME
    off = (s / P) / HOP
    half = N // 2
    fc = np.empty(N)
    fc[:half] = fbin + off
    fc[half:] = fbin + (-off if reverse else off)
    return np.cos(2 * np.pi * np.cumsum(fc) + phase)


def make_member(s_true, decoy_s, sep, decoy_amp, decoy_phase,
                noise_sigma, noise_seed, true_reverse=True):
    """真(固定) + decoy(nuisance) + noise。extractor は触らない。"""
    a = _chirp(s_true, DOM / FRAME, reverse=true_reverse)                 # 固定
    b = _chirp(decoy_s, (DOM + sep) / FRAME, reverse=True, phase=decoy_phase)
    x = a + decoy_amp * b
    if noise_sigma > 0:
        x = x + np.random.default_rng(noise_seed).normal(0, noise_sigma, len(x))
    return x


def single_verdict(x, n_surr=120, seed=11):
    """単体 verdict3 + detected_s。"""
    r = evaluate(x, engine_extract, n_surr=n_surr, seed=seed)
    ab, s_list, _ = abstain_check(x)
    v = verdict3(r['real_d'], r['gap99'], ab)
    return v, int(r['detected_s']), r['real_d'], r['gap99']


def bundle_verdict(members, n_surr=120):
    """M 本の束に v0 規則を適用。返り値 (bundle_verdict, info)。"""
    rows = [single_verdict(x, n_surr=n_surr) for x in members]
    verdicts = [r[0] for r in rows]
    dets = [r[1] for r in rows]
    if all(v == "BLIND" for v in verdicts):
        bv = "BUNDLE_BLIND"
    else:
        votes = [dets[i] for i, v in enumerate(verdicts) if v == "FORWARD"]  # 有効票=FORWARD のみ
        if len(votes) < M:
            bv = "BUNDLE_ABSTAIN"           # 有効票 < M
        elif len(set(votes)) > 1:
            bv = "BUNDLE_ABSTAIN"           # 割れた（unanimous でない）
        else:
            bv = f"BUNDLE_FORWARD(s={votes[0]})"
    return bv, list(zip(verdicts, dets))
