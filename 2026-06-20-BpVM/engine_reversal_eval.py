"""
engine_reversal_eval.py
B13 / reversal_lane engine 配線 harness（#4 凍結 instrument・admissibility 強制版）

locked 仕様（v159 入口条件）:
  raw time-series → admit() → (residue_walk, detected_s)
                 → directional_localization_score(walk, detected_s)
                 → phase_shuffle(raw) で surrogate → 99pct 天井
                 → #4 directional baseline の SNR→gap99 曲線に重ねる
プラグ点は一つ: engine_extract(signal) -> (residue_walk, detected_s)
失敗様式は real 盲目（条件6）。surrogate は暴れない＝偽陽性ではない。
"""
import numpy as np
from test_reversal_isolation import phase_shuffle
from test_reversal_directional import directional_localization_score
from engine_admissibility import admit, Inadmissible

CEIL_PCT = 99

# 条件5: #4 directional baseline（test_reversal_margin_sweep_directional.py の出力・固定）
# (SNR_dB, gap99) — engine の動作 SNR をここに重ねて読む
DIRECTIONAL_BASELINE = [
    (np.inf, 0.537), (23.0, 0.610), (17.0, 0.618), (11.0, 0.602),
    (6.1, 0.644), (3.0, 0.610), (-0.5, 0.469), (-3.0, 0.318), (-6.5, 0.107),
]


def evaluate(signal, engine_extract, n_surr=400, seed=0):
    """engine 信号一本に #4 を回す。s は engine 検出から取る（条件3,4）。"""
    walk, s = admit(signal, engine_extract)          # 条件1,3,4 を入口で強制
    real_d, n_runs, cov = directional_localization_score(walk, s)
    rng = np.random.default_rng(seed)
    surr = np.array([                                 # 条件2: null は raw 波形 phase-shuffle
        directional_localization_score(admit(phase_shuffle(signal, rng), engine_extract)[0], s)[0]
        for _ in range(n_surr)
    ])
    ceil99 = np.percentile(surr, CEIL_PCT)
    return dict(
        detected_s=s, real_d=real_d, runs=n_runs, cov=cov,
        surr_99=ceil99, surr_max=surr.max(), surr_fracpos=(surr > 0).mean(),
        gap99=real_d - ceil99, gapmax=real_d - surr.max(),
        verdict=_verdict(real_d, ceil99),
    )


def _verdict(real_d, ceil99):
    # 条件6: 失敗は real 盲目として扱う
    if real_d >= 0.9 and real_d > ceil99:
        return "前進: real_d~1.0 & gap99>0（reversal_lane §4 候補）"
    if real_d < 0.9:
        return "盲目: real_d 低下（surrogate 不変・偽陽性ではない）→ 落とす"
    return "薄い: real_d は天井超だが 1.0 未満（動作 SNR を baseline と照合）"


def baseline_gap99_at(snr_db):
    """条件5: 動作 SNR を baseline 曲線へ線形補間して重ねる。"""
    xs = [s for s, _ in DIRECTIONAL_BASELINE if np.isfinite(s)][::-1]
    ys = [g for s, g in DIRECTIONAL_BASELINE if np.isfinite(s)][::-1]
    if not np.isfinite(snr_db):
        return DIRECTIONAL_BASELINE[0][1]
    return float(np.interp(snr_db, xs, ys))


if __name__ == "__main__":
    from test_reversal_isolation import make_signal, extract_residue
    print("=" * 72)
    print("engine harness self-test（admissibility 強制・engine = 玩具で配線確認）")
    print("=" * 72)

    def toy_engine(sig):                  # ← 実機 engine の入口に差し替えるだけ
        return extract_residue(sig).astype(int), 1   # (residue_walk, detected_s)

    x = make_signal(1, reverse=True)
    r = evaluate(x, toy_engine)
    print(f"\n[clean] detected_s={r['detected_s']} real_d={r['real_d']:.3f} "
          f"runs={r['runs']} cov={r['cov']:.3f}")
    print(f"        surr_99={r['surr_99']:.3f} gap99={r['gap99']:+.3f}  "
          f"baseline_gap99(inf)={baseline_gap99_at(np.inf):.3f}")
    print(f"        verdict = {r['verdict']}")

    print("\n[reject 確認] 非準拠 engine は入口で弾かれる:")
    try:
        evaluate(x, lambda s: (extract_residue(s).astype(int), None))
    except Inadmissible as e:
        print(f"        REJECT: {e}")
