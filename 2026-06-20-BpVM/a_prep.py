"""
a_prep.py
A-prep : 実観測 signal bundle の入口 harness（v164 仕様そのまま・新防御なし）。

役割（director 固定）:
  - 同一 target の実観測 bundle（まず K=3）を受け、各本に single_verdict/detected_s/real_d/guard_reason を出す。
  - 全本 FORWARD かつ detected_s 一致 → BUNDLE_FORWARD(s)／割れたら BUNDLE_ABSTAIN。
  - ログに必ず prior 宣言を刻む。
  - 実データなしに real_d の成功は書かない。玩具 self-test は *配線検査だけ*。

注意: BUNDLE_FORWARD は真の保証ではない（v162→v164）。意味は cross-signal stable dominant s。
"""
import numpy as np
from engine_admissibility import check_timeseries, Inadmissible
from engine_reversal_eval import evaluate
from engine_b13 import engine_extract
from guard_b13 import abstain_check, verdict3

K_EXPECT = 3

# director 固定・必ず刻む prior 宣言文
PRIOR_DECLARATION = (
    "BUNDLE_FORWARD は真の保証ではない。これは cross-signal stable dominant s の観測であり、"
    "偽 decoy が自然観測条件で s を保つ場合は通過し得る。"
    "B-4 により、この防御は decoy_s 自然変動 prior に依存する。"
)


def _verdict_with_reason(real_d, gap99, abstain, s_list):
    """ラベルは guard_b13.verdict3 に一本化（単一ソース・凍結 instrument）。
    ここは reason/diag 文字列を付すだけ。判定ロジックを再実装しない。
    （等価性は verify_verdict_equiv.py で網羅検査済み：160セル・境界含め不一致ゼロ。）"""
    v = verdict3(real_d, gap99, abstain)                 # ← 判定は委譲
    if v == "BLIND":
        reason = (f"real_d={real_d:.3f}<0.9 (real-collapse)" if real_d < 0.9
                  else f"gap99={gap99:+.3f}<=0 (null 未分離)")
    elif v == "ABSTAIN":
        reason = f"multibin s split {s_list}"
    else:
        reason = "clean (real_d>=0.9 & gap99>0 & single-bin s)"
    return v, reason


def single_member(signal, n_surr=200, seed=11):
    """実観測一本。条件1（raw 時系列）を入口で強制してから測る。"""
    check_timeseries(signal)                                  # 非波形は門前で弾く
    r = evaluate(signal, engine_extract, n_surr=n_surr, seed=seed)
    ab, s_list, _ = abstain_check(signal)
    v, reason = _verdict_with_reason(r['real_d'], r['gap99'], ab, s_list)
    return dict(verdict=v, detected_s=int(r['detected_s']),
                real_d=float(r['real_d']), gap99=float(r['gap99']), guard_reason=reason)


def run_bundle(signals, target="(unnamed target)", n_surr=200):
    """実観測 bundle 入口。v164 仕様で BUNDLE verdict を吐き、prior を必ず刻む。"""
    K = len(signals)
    rows = [single_member(s, n_surr=n_surr) for s in signals]

    verdicts = [r['verdict'] for r in rows]
    dets = [r['detected_s'] for r in rows]
    all_fwd = all(v == "FORWARD" for v in verdicts)
    s_unanimous = len(set(dets)) == 1
    if all_fwd and s_unanimous:
        bundle = f"BUNDLE_FORWARD(s={dets[0]})"
    else:
        bundle = "BUNDLE_ABSTAIN"
    all_blind = all(v == "BLIND" for v in verdicts)   # 診断のみ（headline は変えない）

    # ---- ログ（prior を必ず刻む）----
    print("-" * 84)
    print(f"[A-prep] target = {target}   K = {K}" + ("" if K == K_EXPECT else f"  (※ K_EXPECT={K_EXPECT})"))
    for i, r in enumerate(rows):
        print(f"  obs{i}: {r['verdict']:<8} s={r['detected_s']:<2} "
              f"real_d={r['real_d']:.3f} gap99={r['gap99']:+.3f}  reason: {r['guard_reason']}")
    print(f"  => {bundle}")
    if all_blind:
        print("  [diag] 全本 BLIND（反転なし）。s-曖昧でなく無反転——headline は ABSTAIN だが意味は別。")
    print(f"  [PRIOR] {PRIOR_DECLARATION}")
    print("-" * 84)
    return dict(target=target, K=K, bundle=bundle, rows=rows,
                all_blind=all_blind, prior=PRIOR_DECLARATION)


# =========================================================
#  WIRING SELF-TEST ONLY — 配線検査のみ。実 real_d の成功ではない。
#  玩具 bundle で「入口が正しい verdict へ配線されているか」だけを見る。
# =========================================================
if __name__ == "__main__":
    from bundle_b13 import make_member
    print("=" * 84)
    print("A-prep WIRING SELF-TEST（玩具・配線検査のみ／実データの結果ではない）")
    print("=" * 84)

    S_TRUE = 1
    def bundle_true_dominant():       # 真支配・decoy 弱遠・s 振る → BUNDLE_FORWARD(1) へ配線されるはず
        return [make_member(S_TRUE, ds, 6, 0.25, ph, 0.05, sd)
                for ds, ph, sd in [(9,0.0,1),(11,1.3,2),(4,2.6,3)]]
    def bundle_decoy_varied():        # decoy 支配・s 振る → BUNDLE_ABSTAIN へ配線されるはず
        return [make_member(S_TRUE, ds, sp, 2.3, ph, 0.05, sd)
                for ds, sp, ph, sd in [(12,5,0.0,1),(4,2,1.3,2),(7,3,2.6,3)]]
    def bundle_steady():              # 反転なし → 全 BLIND → ABSTAIN(+diag) へ配線されるはず
        return [make_member(S_TRUE, 5, 3, 0.0, 0.0, 0.05, sd, true_reverse=False) for sd in (1,2,3)]
    def bundle_invariant_decoy():     # 不変 decoy 支配 → BUNDLE_FORWARD(5 誤) へ配線されるはず（prior 境界の確認）
        return [make_member(S_TRUE, 5, 3, 2.5, ph, 0.05, sd) for ph, sd in [(0.0,1),(1.3,2),(2.6,3)]]

    for name, mk in [("真支配(期待 FORWARD s=1)", bundle_true_dominant),
                     ("decoy 振る(期待 ABSTAIN)", bundle_decoy_varied),
                     ("反転なし(期待 ABSTAIN+diag)", bundle_steady),
                     ("不変 decoy(期待 FORWARD s=5 誤・prior 境界)", bundle_invariant_decoy)]:
        print(f"\n### wiring: {name}")
        run_bundle(mk(), target=f"TOY::{name}", n_surr=150)

    print("\n注意: 上は配線検査のみ。BUNDLE_FORWARD が出ても玩具信号であり、実機 real_d の成功ではない。")
    print("      実データ投入は director の手番。run_bundle(実観測リスト, target) に差すだけ。")
