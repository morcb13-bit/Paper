"""
test_reversal_margin_sweep_directional.py
B13 / reversal_lane #5 : margin sweep を locked instrument(#4 directional) で引き直す。

#3 (test_reversal_margin_sweep.py) と同じ口・同じ sigma 格子・同じ surrogate 手順。
差替えは readout だけ: localization_score(無向き) → directional_localization_score(向き忠実)。
判定天井は §5 規約どおり 99pct（max は監査ログ）。

これは engine の比較基準。engine の real_d が SNR 劣化下でこの曲線の近傍を保てば前進、
落ちれば「盲目」として落とす（surrogate は暴れない）。緑化はしない。曲線を読むだけ。
"""
import numpy as np
from test_reversal_isolation import make_signal, phase_shuffle, extract_residue
from test_reversal_directional import directional_localization_score

S = 1                      # 最薄の s=1 を測る（#3 と同じ）
N_NOISE = 20
N_SURR  = 200

def measure_at(sigma, seed=0):
    rng = np.random.default_rng(seed)
    x0 = make_signal(S, reverse=True)
    sig_rms = np.sqrt(np.mean(x0**2))
    real_scores, surr_maxes, surr_99s = [], [], []
    for _ in range(N_NOISE):
        n = rng.normal(0, sigma, size=len(x0))
        xn = x0 + n
        rs, _, _ = directional_localization_score(extract_residue(xn), S)
        real_scores.append(rs)
        sm = np.array([directional_localization_score(extract_residue(phase_shuffle(xn, rng)), S)[0]
                       for _ in range(N_SURR)])
        surr_maxes.append(sm.max()); surr_99s.append(np.percentile(sm, 99))
    real_scores = np.array(real_scores)
    snr_db = 20*np.log10(sig_rms/sigma) if sigma > 0 else np.inf
    return dict(
        sigma=sigma, snr_db=snr_db,
        real_mean=real_scores.mean(), real_min=real_scores.min(),
        surr_max=max(surr_maxes), surr_99=np.mean(surr_99s),
        gap99=real_scores.mean() - np.mean(surr_99s),
        gapmax=real_scores.mean() - max(surr_maxes),
    )

if __name__ == "__main__":
    print("="*86)
    print(f"reversal_lane margin sweep [DIRECTIONAL #4]  (s={S}, {N_NOISE} noise x {N_SURR} surr each)")
    print("="*86)
    print(f"{'sigma':>6} {'SNR(dB)':>8} {'real_mean':>10} {'real_min':>9} "
          f"{'surr_99':>8} {'surr_max':>9} {'GAP99':>7} {'GAPmax':>7}")
    print("-"*86)
    sigmas = [0.0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.00, 1.50]
    curve = []
    for i, sg in enumerate(sigmas):
        r = measure_at(sg, seed=100+i)
        curve.append((sg, r['snr_db'], r['gap99']))
        snr = f"{r['snr_db']:.1f}" if np.isfinite(r['snr_db']) else "inf"
        print(f"{sg:>6.2f} {snr:>8} {r['real_mean']:>10.3f} {r['real_min']:>9.3f} "
              f"{r['surr_99']:>8.3f} {r['surr_max']:>9.3f} {r['gap99']:>7.3f} {r['gapmax']:>7.3f}")
    print("-"*86)
    print("\ngap99 曲線 (sigma -> gap@99pct):")
    gmax = max(g for _,_,g in curve) or 1.0
    for sg, snr, g in curve:
        bar = "#" * int(round(max(g,0)/gmax*40))
        snrs = f"{snr:5.1f}" if np.isfinite(snr) else "  inf"
        print(f"  sigma={sg:>4.2f} SNR={snrs}dB | {bar} {g:.3f}")
    print("\n判定: gap99>0 を保つ SNR 範囲が engine の動作点を覆えば、#4 は SNR 下でも物差しとして立つ。")
