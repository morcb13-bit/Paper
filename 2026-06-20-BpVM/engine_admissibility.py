"""
engine_admissibility.py
B13 / reversal_lane engine-wiring admissibility gate（director 固定の6条件・実行可能版）

engine コードが来たとき、非準拠 engine を入口で弾く。誤った real_d を黙って出させない。

固定条件（v159 入口条件）:
  1. 入力は時系列波形であること（residue 列・特徴量は不可）。
  2. null は raw 波形への phase-shuffle で作ること。
  3. engine は (residue_walk, detected_s) を返すこと。
  4. detected_s が無い場合、全 s 探索は不可。別実験として天井再較正。
  5. 判定は #4 directional baseline の SNR→gap99 曲線に重ねる。
  6. 失敗様式は real 盲目として扱う。偽陽性ではない。
"""
import numpy as np

P = 13
FRAME = 256


class Inadmissible(Exception):
    pass


# ---- 条件1: 入力は時系列波形 ----
def check_timeseries(signal):
    if not isinstance(signal, np.ndarray):
        raise Inadmissible("条件1違反: 入力が np.ndarray でない（時系列波形が必要）")
    if signal.ndim != 1:
        raise Inadmissible(f"条件1違反: 入力が1次元でない（shape={signal.shape}）")
    if not np.issubdtype(signal.dtype, np.floating):
        raise Inadmissible(f"条件1違反: 入力が float でない（dtype={signal.dtype}）"
                           "／residue 列や整数特徴量は不可")
    if signal.size < 4 * FRAME:
        raise Inadmissible(f"条件1違反: 入力長 {signal.size} が短すぎる（波形でなく特徴量の疑い）")
    # residue 列が float 化されただけのものを弾く: 値域が [0,P) の整数値ばかりなら疑う
    vals = signal[np.isfinite(signal)]
    if vals.size and np.all(vals == np.round(vals)) and vals.min() >= 0 and vals.max() < P:
        raise Inadmissible("条件1違反: 入力が residue 列（0..12 の整数値）に見える。raw 波形が必要")
    return True


# ---- 条件3,4: engine は (residue_walk, detected_s) を返す ----
def check_engine_return(ret):
    if not (isinstance(ret, tuple) and len(ret) == 2):
        raise Inadmissible("条件3違反: engine_extract は (residue_walk, detected_s) の2要素を返すこと")
    walk, s = ret
    walk = np.asarray(walk)
    if not np.issubdtype(walk.dtype, np.integer):
        raise Inadmissible(f"条件3違反: residue_walk が整数列でない（dtype={walk.dtype}）")
    if walk.size == 0 or walk.min() < 0 or walk.max() >= P:
        raise Inadmissible(f"条件3違反: residue_walk の値域が [0,{P}) でない")
    # 条件4: detected_s 必須・全 s 探索不可
    if s is None:
        raise Inadmissible("条件4違反: detected_s が None。全 s 探索は不可。"
                           "別実験として多重比較込みで天井を再較正すること")
    s = int(s)
    if not (1 <= s <= P - 1):
        raise Inadmissible(f"条件4違反: detected_s={s} が 1..{P-1} の外（s=0 は readout 対象外）")
    return walk.astype(int), s


def admit(signal, engine_extract):
    """6条件のうち静的に検査できる 1,3,4 を入口で強制。通れば (walk, s) を返す。"""
    check_timeseries(signal)                 # 条件1
    ret = engine_extract(signal)             # 条件2は harness が phase_shuffle(signal) で担保
    walk, s = check_engine_return(ret)       # 条件3,4
    return walk, s                            # 条件5,6 は判定段（harness）で適用


# =========================================================
#  SELF-TEST: 準拠 engine を ACCEPT・非準拠を REJECT
# =========================================================
if __name__ == "__main__":
    from test_reversal_isolation import make_signal, extract_residue

    print("=" * 64)
    print("admissibility gate self-test")
    print("=" * 64)

    raw = make_signal(1, reverse=True)

    def good_engine(sig):
        return extract_residue(sig).astype(int), 1          # (walk, detected_s)

    def bad_no_s(sig):
        return extract_residue(sig).astype(int), None       # detected_s 無し

    def bad_walk_only(sig):
        return extract_residue(sig).astype(int)             # tuple でない

    cases = [
        ("準拠 engine（波形→walk,s）", raw, good_engine, True),
        ("非準拠: detected_s=None（全s探索）", raw, bad_no_s, False),
        ("非準拠: walk のみ返す", raw, bad_walk_only, False),
        ("非準拠: 入力が residue 列", extract_residue(raw).astype(float), good_engine, False),
    ]
    for name, sig, eng, expect_ok in cases:
        try:
            walk, s = admit(sig, eng)
            got = f"ADMIT (s={s}, walk_len={len(walk)})"
            ok = expect_ok
        except Inadmissible as e:
            got = f"REJECT: {e}"
            ok = not expect_ok
        mark = "✓" if ok else "✗ 期待外"
        print(f"  [{mark}] {name}\n        → {got}")
