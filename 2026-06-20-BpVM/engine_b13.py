"""
engine_b13.py
B13 / reversal_lane 実エンジン本体（v159 入口条件準拠・自己検出版）

玩具との差は一点だけだが本質的:
  玩具 good_engine : (extract_residue(sig), 1)        # detected_s をハードコード
  本エンジン       : (extract_residue(sig), detect_s) # s を walk から自己検出

これにより #4 は「検出した s」で採点する（条件3,4 のとおり）。
検出が外れれば directional readout は間違った対蹠を探し real_d が落ちる——
それを隠さない。出た数字をそのまま読む。

extract_residue / stable_runs / P は凍結 fixture を import で共有（触らない）。
numpy only.
"""
import numpy as np
from test_reversal_isolation import P, extract_residue
from test_reversal_localization import stable_runs


def detect_s(walk, min_run=3):
    """walk から s を自己検出。
    規則: 先行する安定ラン va に *続いて* 対蹠 vb==(-va)%P が来る対のうち、
          合計被覆が最大の対の va を s とする（向きを保つため va は先行側）。
          対蹠対が無ければ最頻の非零 residue へフォールバック（必ず 1..12 を返す）。
    surrogate でも必ず妥当な s を返す（admit が 1..12 を要求するため）。
    """
    segs = stable_runs(walk, min_run)           # [(v, l, pos), ...]
    best_len, best_s = 0, None
    for (va, la, _), (vb, lb, _) in zip(segs, segs[1:]):
        if va != 0 and vb == (-va) % P:
            if la + lb > best_len:
                best_len, best_s = la + lb, va
    if best_s is not None:
        return int(best_s)
    # フォールバック: 最頻の非零 residue（surrogate / 反転無し でも妥当値）
    nz = walk[walk != 0]
    if nz.size == 0:
        return 1
    s = int(np.bincount(nz, minlength=P).argmax())
    return s if 1 <= s <= P - 1 else 1


def engine_extract(signal):
    """v159 プラグ点。raw 波形 -> (residue_walk, detected_s)。
    detected_s は信号から検出（ハードコードしない）。"""
    walk = extract_residue(signal).astype(int)
    s = detect_s(walk)
    return walk, s
