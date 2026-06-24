"""
run_b2_verdict3.py
B-2 検証: 三値 verdict が
  (1) clean 単成分 → FORWARD
  (2) silent miss → ABSTAIN（FORWARD で通さない＝欠陥が塞がる）
  (3) real-collapse / 混濁崩壊 → BLIND
  (4) benign 二成分（真が支配）での ABSTAIN 率＝保守ガードのコスト
を正しく吐くか。silent FORWARD が一件でも残れば B-2 失敗。出た数字をそのまま読む。
"""
import numpy as np
from test_reversal_isolation import make_signal
from engine_reversal_eval import evaluate, baseline_gap99_at
from engine_b13 import engine_extract
from guard_b13 import abstain_check, verdict3
from run_compete_b13 import make_competing, S_TRUE, S_DECOY, P

N_SURR = 150

def run_case(x, n_surr=N_SURR, seed=11):
    r = evaluate(x, engine_extract, n_surr=n_surr, seed=seed)   # real_d, gap99（dominant bin）
    ab, s_list, bins = abstain_check(x)                          # 多ビン s 非一意
    v = verdict3(r['real_d'], r['gap99'], ab)
    return r, ab, s_list, v

print("=" * 92)
print("B-2 三値 verdict 検証（BLIND > ABSTAIN > FORWARD・K=2・amp_ratio=0.5）")
print("=" * 92)

# (1) clean 単成分・全 s → FORWARD 期待
print("\n[1] clean 単成分・全 s → FORWARD 期待")
bad = 0
for s in range(1, P):
    x = make_signal(s, reverse=True)
    r, ab, s_list, v = run_case(x)
    ok = (v == "FORWARD")
    bad += (not ok)
    if not ok:
        print(f"  s={s}: v={v} real_d={r['real_d']:.3f} gap99={r['gap99']:+.3f} s_list={s_list}  ✗")
print(f"  → 全 s FORWARD: {'OK' if bad==0 else f'{bad} 件ずれ'}")

# (2)(4) competing sep=3（別ビン・silent miss が出る所）
print(f"\n[2/4] competing sep=3  s_true={S_TRUE} s_decoy={S_DECOY}（silent miss 域 + benign 域）")
print(f"  {'amp':>5} {'det_s':>5} {'real_d':>7} {'gap99':>7} {'abst':>5} {'s_list':>10}  {'verdict':>8}  真偽")
silent_forward = 0
benign_total = benign_abstain = 0
for amp in [0.0,0.3,0.6,0.8,0.9,1.0,1.1,1.3,1.6,2.0]:
    x = make_competing(S_TRUE, S_DECOY, amp, sep=3)
    r, ab, s_list, v = run_case(x)
    hit = (r['detected_s'] == S_TRUE)         # ground truth（玩具のみ既知）
    # 真偽ラベル: dominant が真か decoy か
    truth = "真支配" if hit else "decoy支配"
    # silent miss = dominant 誤 & real_d 健全 & verdict が FORWARD なら欠陥が残存
    if (not hit) and r['real_d'] >= 0.9 and v == "FORWARD":
        silent_forward += 1
    # benign = dominant が真（hit）かつ amp>0（競合存在）→ ここでの ABSTAIN は保守コスト
    if hit and amp > 0:
        benign_total += 1; benign_abstain += (v == "ABSTAIN")
    print(f"  {amp:>5.2f} {r['detected_s']:>5} {r['real_d']:>7.3f} {r['gap99']:>+7.3f} "
          f"{str(ab):>5} {str(s_list):>10}  {v:>8}  {truth}")
print(f"  → silent miss が FORWARD で残った件数: {silent_forward}（0 が必須）")
print(f"  → benign 競合域での ABSTAIN（保守コスト）: {benign_abstain}/{benign_total}")

# (3) sep=0 混濁崩壊 → BLIND 期待
print(f"\n[3] competing sep=0（混濁）→ 崩壊は BLIND 期待")
print(f"  {'amp':>5} {'real_d':>7} {'gap99':>7} {'abst':>5}  {'verdict':>8}")
for amp in [0.3,0.6,0.9,1.2]:
    x = make_competing(S_TRUE, S_DECOY, amp, sep=0)
    r, ab, s_list, v = run_case(x)
    print(f"  {amp:>5.2f} {r['real_d']:>7.3f} {r['gap99']:>+7.3f} {str(ab):>5}  {v:>8}")

# 負対照
print(f"\n[neg] steady（反転なし）→ BLIND 期待")
r, ab, s_list, v = run_case(make_signal(3, reverse=False))
print(f"  real_d={r['real_d']:.3f} gap99={r['gap99']:+.3f} abstain={ab} verdict={v}")

print("\n" + "=" * 92)
print("B-2 合否: silent miss が FORWARD で 0 件なら、二値 verdict の欠陥が塞がった。")
print("ABSTAIN は失敗でなく honesty 出力（構造が単一 s を強制していない）。")
print("=" * 92)
