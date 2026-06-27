"""
guard_b13.py
B-2 : 棄権ガード（多ビン s 非一意検出）＋三値 verdict。
director 固定の最小定義のみ実装。広げない。

ABSTAIN 定義（暫定・保守）:
  上位 K ビンを読む。bin0 は常に読む。2位以降は「mean|X| が 1位の amp_ratio 倍以上」の時だけ有効ビン。
  有効ビンが 2 本以上あり detected_s が一致しなければ → ABSTAIN（単一 s 非一意）。
  K=2, amp_ratio=0.5（2位が1位の50%以上なら読む）。

三値 verdict（precedence: BLIND > ABSTAIN > FORWARD）:
  BLIND  : real_d < 0.9  または gap99 <= 0           （崩壊は先に loud で落とす）
  ABSTAIN: 上記でなく、多ビン s 不一致               （健全域での silent miss を声出し）
  FORWARD: real_d>=0.9 かつ gap99>0 かつ non-abstain

ABSTAIN は失敗でなく honesty 出力＝「構造が単一答を強制していない」。
extract 系は凍結 fixture と同一の口（mean|X| ランキング・同じ het 式）。numpy only.
"""
import numpy as np
from test_reversal_isolation import P, FRAME, HOP
from engine_b13 import detect_s

K_DEFAULT = 2
AMP_RATIO_DEFAULT = 0.5


def _bin_profile(x):
    """凍結 extract_residue と同一の STFT。位相(frames×bins) と mean|X|(bins) を返す。"""
    win = np.hanning(FRAME)
    starts = range(0, len(x) - FRAME + 1, HOP)
    ph, mags = [], []
    for st in starts:
        Xf = np.fft.rfft(x[st:st + FRAME] * win)
        ph.append(np.angle(Xf)); mags.append(np.abs(Xf))
    return np.array(ph), np.array(mags).mean(0)        # mean|X| over frames


def _walk_at(ph, b):
    """ビン b の位相増分 -> residue walk（凍結 het 式と同一）。"""
    het = np.diff(np.unwrap(ph[:, b])) - 2 * np.pi * b / FRAME * HOP
    return np.round(het / (2 * np.pi / P)).astype(int) % P


def abstain_check(x, K=K_DEFAULT, amp_ratio=AMP_RATIO_DEFAULT):
    """多ビン s 非一意なら True。返り値 (abstain, s_list, valid_bins)。"""
    ph, mag = _bin_profile(x)
    mag = mag.copy(); mag[0] = 0.0                      # DC 除外（凍結口と同じ）
    order = np.argsort(mag)[::-1]                        # 強い順のビン
    top0 = mag[order[0]] + 1e-12
    s_list, used = [], []
    for rank, b in enumerate(order[:K]):
        if rank == 0 or mag[b] >= amp_ratio * top0:     # 1位は常に・2位以降は閾値超のみ
            s_list.append(detect_s(_walk_at(ph, b)))
            used.append(int(b))
    abstain = len(set(s_list)) >= 2
    return abstain, s_list, used


def verdict3(real_d, gap99, abstain):
    if real_d < 0.9 or gap99 <= 0:
        return "BLIND"
    if abstain:
        return "ABSTAIN"
    return "FORWARD"
