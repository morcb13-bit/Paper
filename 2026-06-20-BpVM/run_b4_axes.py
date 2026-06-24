"""
run_b4_axes.py
B-4 : s-不変交絡を壊せる軸が decoy_s 以外にあるか。束定義は B-3 のまま・decoy_s 固定。
sep / phase / amp / seed / window / bin-neighborhood を個別に振り、
強支配・不変 decoy(s=5) が BUNDLE_FORWARD(誤) を維持するかを見る。

切り分け（重要）:
  ・S-SPLIT   : detected_s が束内で割れた → その軸が *不変 decoy の見かけ s を壊した*（新破壊軸）
  ・INCIDENTAL: detected_s は割れていないが単体 ABSTAIN/BLIND で有効票<M → 既存ガード由来（新軸でない）
  ・MAINTAINED: BUNDLE_FORWARD(誤) 維持 → その軸は無力
探索して無ければ「無い」と明記して終える。防御を作ったふりはしない。
"""
import numpy as np
from test_reversal_isolation import P, FRAME, HOP, DOM
from engine_reversal_eval import evaluate
from engine_b13 import detect_s
from guard_b13 import abstain_check, verdict3
from bundle_b13 import M, make_member

S_TRUE, DECOY_S = 1, 5
NS = 120

# ---- 可変 extractor（観測側 nuisance・凍結 instrument の定義は変えない別観測者）----
_WINS = {"hann": np.hanning, "hamming": np.hamming, "blackman": np.blackman, "box": np.ones}
def make_extractor(window="hann", bin_off=0):
    win = _WINS[window](FRAME)
    def extract(signal):
        starts = range(0, len(signal) - FRAME + 1, HOP)
        ph, mags = [], []
        for st in starts:
            X = np.fft.rfft(signal[st:st + FRAME] * win)
            ph.append(np.angle(X)); mags.append(np.abs(X))
        ph = np.array(ph); mags = np.array(mags)
        dom = int(np.argmax(mags[:, 1:].mean(0)) + 1 + bin_off)
        dom = max(1, min(dom, ph.shape[1] - 1))
        het = np.diff(np.unwrap(ph[:, dom])) - 2 * np.pi * dom / FRAME * HOP
        walk = (np.round(het / (2 * np.pi / P)).astype(int) % P)
        return walk, detect_s(walk)
    return extract

def single(signal, extractor, seed=11):
    r = evaluate(signal, extractor, n_surr=NS, seed=seed)
    ab, _, _ = abstain_check(signal)                 # ガードは凍結 hann（保守・誤 FWD を作れない）
    v = verdict3(r['real_d'], r['gap99'], ab)
    return v, int(r['detected_s'])

def bundle_and_classify(rows):
    verdicts = [v for v, _ in rows]; dets = [d for _, d in rows]
    fwd_s = [dets[i] for i, v in enumerate(verdicts) if v == "FORWARD"]
    if all(v == "BLIND" for v in verdicts):
        bv = "BUNDLE_BLIND"
    elif len(fwd_s) < M:
        bv = "BUNDLE_ABSTAIN"
    elif len(set(fwd_s)) > 1:
        bv = f"BUNDLE_ABSTAIN"
    else:
        bv = f"BUNDLE_FORWARD(s={fwd_s[0]})"
    # 機構ラベル
    s_split = len(set(dets)) > 1
    if bv.startswith("BUNDLE_FORWARD"):
        mech = "MAINTAINED（軸は無力）"
    elif s_split:
        mech = "S-SPLIT（★この軸が不変 decoy の見かけ s を壊した）"
    else:
        mech = "INCIDENTAL（s 不変・単体ガード/real_d 由来・新軸でない）"
    return bv, dets, verdicts, mech

frozen = make_extractor("hann", 0)
DA = 2.5  # 強支配

def members_signal_axis(axis):
    """信号側で 1 軸だけ振る。decoy_s=5 固定・他は固定。"""
    base = dict(sep=3, decoy_phase=0.0, decoy_amp=DA, noise_sigma=0.0, noise_seed=0)
    variants = {
        "sep":   [dict(base, sep=3), dict(base, sep=4), dict(base, sep=5)],
        "phase": [dict(base, decoy_phase=0.0), dict(base, decoy_phase=2.1), dict(base, decoy_phase=4.2)],
        "amp":   [dict(base, decoy_amp=2.0), dict(base, decoy_amp=2.5), dict(base, decoy_amp=3.0)],
        "seed":  [dict(base, noise_sigma=0.05, noise_seed=1),
                  dict(base, noise_sigma=0.05, noise_seed=2),
                  dict(base, noise_sigma=0.05, noise_seed=3)],
    }[axis]
    sigs = [make_member(S_TRUE, DECOY_S, p["sep"], p["decoy_amp"], p["decoy_phase"],
                        p["noise_sigma"], p["noise_seed"]) for p in variants]
    return [single(s, frozen) for s in sigs]

def members_obs_axis(axis):
    """観測側で 1 軸だけ振る。同一信号に別 extractor。"""
    sig = make_member(S_TRUE, DECOY_S, 3, DA, 0.0, 0.0, 0)
    exts = {
        "window":           [make_extractor("hann",0), make_extractor("hamming",0), make_extractor("blackman",0)],
        "bin-neighborhood": [make_extractor("hann",-1), make_extractor("hann",0), make_extractor("hann",+1)],
    }[axis]
    return [single(sig, e) for e in exts]

print("=" * 92)
print(f"B-4 : s-不変交絡(decoy_s={DECOY_S} 固定・amp={DA} 強支配) を壊せる軸はあるか  s_true={S_TRUE}")
print("=" * 92)
print(f"  {'axis':>16} {'detected_s':>16} {'verdicts':>22} {'bundle':>20}  機構")
results = {}
for axis in ["sep", "phase", "amp", "seed"]:
    rows = members_signal_axis(axis)
    bv, dets, vs, mech = bundle_and_classify(rows); results[axis] = mech
    print(f"  {axis:>16} {str(dets):>16} {str(vs):>22} {bv:>20}  {mech}")
for axis in ["window", "bin-neighborhood"]:
    rows = members_obs_axis(axis)
    bv, dets, vs, mech = bundle_and_classify(rows); results[axis] = mech
    print(f"  {axis:>16} {str(dets):>16} {str(vs):>22} {bv:>20}  {mech}")

print("\n" + "=" * 92)
splitters = [a for a, m in results.items() if m.startswith("S-SPLIT")]
if splitters:
    print(f"結論: decoy_s 以外に s-不変 decoy を壊す軸が見つかった: {splitters}")
else:
    print("結論: sep/phase/amp/seed/window/bin-neighborhood のいずれも S-SPLIT を起こさない。")
    print("      → B-3 の防御軸は decoy_s 変動のみ、と確定。これを A へ仕様として渡す。")
    print("      （棄権が出た軸があってもそれは既存ガード由来＝INCIDENTAL で、新破壊軸ではない。）")
print("=" * 92)
