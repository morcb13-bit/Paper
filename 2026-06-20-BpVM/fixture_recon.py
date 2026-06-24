"""
fixture_recon.py — v157 §1 仕様からの fixture 再構成（読み・再構成グレード）

信号: dom bin=40・frame=256・hop=256(無重なり)・各半12フレーム
      residue 列 = [s]*12 → [(-s)%13]*12 のクリーン対蹠反転
抽出: 各フレーム rfft → bin40 の位相を 13 量子化 → residue
surrogate: 全信号を位相シャッフル（|.|保存・位相乱択）→ 再フレーム → 抽出

まず surrogate の run 構造を実見し、v157 の
  「位相シャッフルは長ラン・両 residue を保つ」「最悪 surr 被覆 ~0.75」
に質的に較正できるかを確かめる。
"""
import numpy as np

MOD = 13
FBIN = 40
FRAME = 256
NHALF = 12
NFRAMES = 2 * NHALF
MIN_RUN = 3


def make_signal(s, rng=None, noise=0.0):
    r = [s % MOD] * NHALF + [(-s) % MOD] * NHALF
    x = np.empty(NFRAMES * FRAME)
    n = np.arange(FRAME)
    for k, rk in enumerate(r):
        theta = 2 * np.pi * rk / MOD
        x[k * FRAME:(k + 1) * FRAME] = np.cos(2 * np.pi * FBIN * n / FRAME + theta)
    if noise > 0 and rng is not None:
        x = x + noise * rng.standard_normal(x.size)
    return x


def extract_residues(x):
    out = []
    for k in range(NFRAMES):
        frame = x[k * FRAME:(k + 1) * FRAME]
        X = np.fft.rfft(frame)
        ph = np.angle(X[FBIN])
        res = int(np.round(ph * MOD / (2 * np.pi))) % MOD
        out.append(res)
    return out


def phase_shuffle(x, rng):
    X = np.fft.rfft(x)
    mag = np.abs(X)
    rand_ph = rng.uniform(-np.pi, np.pi, size=X.shape)
    rand_ph[0] = 0.0
    if x.size % 2 == 0:
        rand_ph[-1] = 0.0  # Nyquist real
    Xs = mag * np.exp(1j * rand_ph)
    return np.fft.irfft(Xs, n=x.size)


def rle(seq, min_run=MIN_RUN):
    """run-length encode、min_run 未満のランは落とす。被覆分母は全フレーム数 NFRAMES。"""
    runs = []
    if not seq:
        return runs
    cur, ln = seq[0], 1
    for v in seq[1:]:
        if v == cur:
            ln += 1
        else:
            runs.append((cur, ln))
            cur, ln = v, 1
    runs.append((cur, ln))
    return [(r, l) for (r, l) in runs if l >= min_run]


if __name__ == "__main__":
    rng = np.random.default_rng(13)

    # --- clean 抽出が仕様どおりか（s=1..3 で狙い通り、s=4 は遷移許容） ---
    print("clean 抽出チェック:")
    for s in [1, 2, 3, 4]:
        res = extract_residues(make_signal(s))
        runs = rle(res)
        print(f"  s={s}: runs={runs}  expect [{s},{(-s)%MOD}]")

    # --- surrogate の run 構造を実見（s=1, 12本） ---
    print("\nsurrogate run 構造（s=1, phase-shuffle, 先頭12本）:")
    x = make_signal(1)
    longrun_hits = 0
    both_res_hits = 0
    for i in range(12):
        xs = phase_shuffle(x, rng)
        runs = rle(extract_residues(xs))
        has_both = ({r for r, _ in runs} >= {1, 12})
        max_run = max([l for _, l in runs], default=0)
        longrun_hits += (max_run >= 5)
        both_res_hits += has_both
        print(f"  surr#{i:2d}: runs={runs}  both(1&12)={has_both}  maxrun={max_run}")
    print(f"\n  長ラン(>=5)を持つ surrogate: {longrun_hits}/12")
    print(f"  対蹠両 residue を在庫: {both_res_hits}/12")
