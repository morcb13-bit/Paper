"""
test_reversal_isolation.py
B13 / reversal_lane 切分テスト（engine 非依存・玩具）

問い: 位相シャッフル surrogate は s -> -s mod 13 の対蹠反転を本当に壊すのか。
判定: real は反転を保つ / surrogate は反転を失う。この差が出れば reversal lane は前進可。

§2-B の轍を踏まないための固定規律:
  - 触ってよいのは信号生成と residue 抽出だけ（real が狙い通りの walk を持つように）。
  - surrogate（位相シャッフル）・readout（対蹠ペア検出）・null 天井（surrogate 99pct）は固定。
  - 出た数字をそのまま読む。pass するように調整しない。
numpy only.
"""
import numpy as np

P = 13
FRAME = 256
HOP   = 256
N_EACH = 12
DOM   = 40

# ---------- 信号生成（real のみ作り込む） ----------
def make_signal(s, reverse=True):
    n_frames = 2*N_EACH
    N = n_frames*HOP + FRAME
    fbin = DOM/FRAME
    off = (s/P)/HOP
    half = N//2
    fc = np.empty(N)
    fc[:half] = fbin + off
    fc[half:] = fbin + (-off if reverse else off)
    return np.cos(2*np.pi*np.cumsum(fc))

# ---------- 位相シャッフル surrogate（固定・|X| 保存） ----------
def phase_shuffle(x, rng):
    X = np.fft.rfft(x)
    mag = np.abs(X)
    phases = rng.uniform(-np.pi, np.pi, size=len(X))
    phases[0] = 0.0
    if len(x) % 2 == 0:
        phases[-1] = 0.0
    return np.fft.irfft(mag*np.exp(1j*phases), n=len(x))

# ---------- residue 抽出（real / surrogate で同一の口） ----------
def extract_residue(x):
    win = np.hanning(FRAME)
    starts = range(0, len(x)-FRAME+1, HOP)
    ph, mags = [], []
    for st in starts:
        X = np.fft.rfft(x[st:st+FRAME]*win)
        ph.append(np.angle(X)); mags.append(np.abs(X))
    ph = np.array(ph); mags = np.array(mags)
    dom = np.argmax(mags[:,1:].mean(0))+1
    het = np.diff(np.unwrap(ph[:,dom])) - 2*np.pi*dom/FRAME*HOP
    return np.round(het/(2*np.pi/P)).astype(int) % P

# ---------- reversal readout（固定・対蹠ペアのみを読む） ----------
def runs(walk):
    segs, i = [], 0
    while i < len(walk):
        j = i
        while j < len(walk) and walk[j] == walk[i]:
            j += 1
        segs.append((int(walk[i]), j-i)); i = j
    return segs

def reversal_score(walk, min_len=3):
    """s -> -s mod 13（s!=0）の隣接安定ラン対が覆う割合。対蹠でなければ 0。"""
    segs = [(v,l) for (v,l) in runs(walk) if l >= min_len]
    covered = 0
    for (va,la),(vb,lb) in zip(segs, segs[1:]):
        if va != 0 and vb == (-va) % P:
            covered += la + lb
    return covered / len(walk)

# ---------- 較正 + 判定 ----------
def calibrate_and_test(s, n_surr=400, seed=0, label=""):
    rng = np.random.default_rng(seed)
    x = make_signal(s, reverse=True)
    real_walk = extract_residue(x)
    real_score = reversal_score(real_walk)

    surr_scores = []
    for _ in range(n_surr):
        xs = phase_shuffle(x, rng)
        surr_scores.append(reversal_score(extract_residue(xs)))
    surr_scores = np.array(surr_scores)
    ceiling = np.percentile(surr_scores, 99)
    smax = surr_scores.max()
    sep = real_score / (ceiling + 1e-12)

    print(f"--- {label} (s={s}) ---")
    print(f"  real_walk      = {list(map(int, real_walk))}")
    print(f"  real_score     = {real_score:.4f}")
    print(f"  surr mean/99pct/max = {surr_scores.mean():.4f} / {ceiling:.4f} / {smax:.4f}")
    print(f"  surr frac >0   = {(surr_scores>0).mean():.3f}")
    print(f"  separation     = {sep:.2f}x   verdict = "
          f"{'REVERSAL (real>ceiling)' if real_score>ceiling and real_score>0 else 'null'}")
    return real_score, ceiling

# ---------- 負の対照: 反転なし steady walk ----------
def negative_control(s=3, n_surr=400, seed=1):
    rng = np.random.default_rng(seed)
    x = make_signal(s, reverse=False)       # steady [s]*24, 反転なし
    real_walk = extract_residue(x)
    real_score = reversal_score(real_walk)
    surr = np.array([reversal_score(extract_residue(phase_shuffle(x, rng))) for _ in range(n_surr)])
    print(f"--- negative control: steady walk (s={s}, 反転なし) ---")
    print(f"  real_walk  = {list(map(int, real_walk))}")
    print(f"  real_score = {real_score:.4f}   surr mean = {surr.mean():.4f}")
    print(f"  期待: 両方 ~0（steady は coupling でない・readout が誤発火しない）")

if __name__ == "__main__":
    print("="*60)
    print("B13 reversal_lane 切分テスト")
    print("="*60)
    for s in (1, 2, 3):
        calibrate_and_test(s, label=f"reversal")
        print()
    negative_control()
