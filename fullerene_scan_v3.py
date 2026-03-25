"""
fullerene_scan_v3.py  ―― 辻井式位相分解機 v2
==============================================

v2 で判明したこと:
  - 上位5頂点だけで混合の分解精度 相関=0.958
  - 全60頂点使うと 0.458 に落ちる（ノイズ頂点が邪魔）
  - 頂点22: piano=-0.887, violin=+0.014  ← piano専用の識別子

次の実験:
  1. 識別力TOP-K頂点の自動選択（事前知識なし）
     → 混合信号「だけ」から識別力の高い頂点を探せるか？
  2. 時間フレーム分割スキャン
     → onset付近 vs 定常部 vs 減衰部でΦ_dがどう変化するか
  3. 係数推定の安定性（複数フレームで平均）

本質的な問い:
  「混合信号だけを与えたとき、何個の波源がいるか自動推定できるか」
"""

import numpy as np
import math
import scipy.signal as ss
from dataclasses import dataclass
import time

SR  = 22050
PHI = (1 + math.sqrt(5)) / 2


# ════════════════════════════════════════════════
# 音源生成（v7 から移植）
# ════════════════════════════════════════════════

def _adsr(N, sr, attack=0.05, decay=0.1, sustain=0.7, release=0.2):
    env = np.ones(N) * sustain
    a_s, d_s, r_s = int(attack*sr), int(decay*sr), int(release*sr)
    env[:a_s] = np.linspace(0, 1, a_s)
    env[a_s:a_s+d_s] = np.linspace(1, sustain, d_s)
    env[-r_s:] = np.linspace(sustain, 0, r_s)
    return env

def make_piano(f0=220, sr=SR, dur=2.0, onset=0.0):
    N = int(sr*dur); t = np.arange(N)/sr; B = 0.0001
    s = np.zeros(N)
    for h, amp in [(1,1.0),(2,0.6),(3,0.4),(4,0.2),(5,0.12),(6,0.07)]:
        fh = f0*h*math.sqrt(1+B*h**2)
        s += amp*np.exp(-t/(2.0/h))*np.sin(2*math.pi*fh*t)
    s += 0.02*np.random.randn(N)*np.exp(-t/0.02)
    shift = int(onset*sr)
    s = np.concatenate([np.zeros(shift), s])[:N]
    return s / (max(abs(s)) + 1e-9)

def make_violin(f0=440, sr=SR, dur=2.0, onset=0.0):
    N = int(sr*dur); t = np.arange(N)/sr
    vib = 0.004*(1-np.exp(-t/0.3))*np.sin(2*math.pi*5.5*t)
    s = np.zeros(N)
    for h, amp in [(1,1.0),(2,0.4),(3,0.7),(4,0.3),(5,0.5),(6,0.2),(7,0.35)]:
        s += amp*np.sin(2*math.pi*f0*h*(1+vib)*t)
    env = _adsr(N, sr, attack=0.08, decay=0.05, sustain=0.85, release=0.3)
    s = s*env + 0.03*np.random.randn(N)*env
    shift = int(onset*sr)
    s = np.concatenate([np.zeros(shift), s])[:N]
    return s / (max(abs(s)) + 1e-9)

def make_flute(f0=880, sr=SR, dur=2.0, onset=0.0):
    N = int(sr*dur); t = np.arange(N)/sr
    vib = 0.003*(1-np.exp(-t/0.4))*np.sin(2*math.pi*6.0*t)
    s = np.zeros(N)
    for h, amp in [(1,1.0),(2,0.15),(3,0.05),(4,0.02)]:
        s += amp*np.sin(2*math.pi*f0*h*(1+vib)*t)
    air = 0.04*np.random.randn(N)
    b, a = ss.butter(4, [f0/sr*1.5, min(f0/sr*8, 0.9)], btype='band')
    s += ss.lfilter(b, a, air)
    env = _adsr(N, sr, attack=0.05, decay=0.02, sustain=0.9, release=0.25)
    s = s*env
    shift = int(onset*sr)
    s = np.concatenate([np.zeros(shift), s])[:N]
    return s / (max(abs(s)) + 1e-9)


# ════════════════════════════════════════════════
# 辻井式 Φ_d
# ════════════════════════════════════════════════

def bandpass_safe(x, f_c, bw, sr):
    nyq = sr / 2
    lo = max(20.0, f_c - bw/2)
    hi = min(nyq*0.98, f_c + bw/2)
    if hi <= lo + 10:
        return x
    b, a = ss.butter(2, [lo/nyq, hi/nyq], btype='band')
    return ss.lfilter(b, a, x)

def phi_d_tsuji(signal, tau, f_c, bw, sr=SR):
    filtered = bandpass_safe(signal, f_c, bw, sr)
    z = ss.hilbert(filtered)
    tau = max(1, min(tau, len(z)//4))
    z_early = z[:-tau]
    z_late  = z[tau:]
    ratio     = z_late / (z_early + 1e-12)
    delta_phi = np.angle(ratio)
    pos = delta_phi[delta_phi >  0.01]
    neg = delta_phi[delta_phi < -0.01]
    if len(pos) < 10 or len(neg) < 10:
        return 0.0
    pm = pos.mean(); nm = abs(neg.mean())
    return float((pm - nm) / (pm + nm + 1e-12))


@dataclass
class ScanVertex:
    idx: int
    R: float; lam: float; rotation: float
    tau: int = 0; f_c: float = 0.0; bw: float = 0.0
    def __post_init__(self):
        self.tau = max(1, int(self.R * 220))
        F_MIN, F_MAX = 100.0, 8000.0
        self.f_c = F_MIN * (F_MAX/F_MIN) ** (self.lam / (2*math.pi))
        self.bw = 50.0 + 1950.0 * (self.rotation / math.pi)


def fullerene_vertices(n=60):
    ga = math.pi * (3 - math.sqrt(5))
    verts = []
    for i in range(n):
        y = 1 - (i/(n-1))*2
        r_xy = math.sqrt(max(0.0, 1-y**2))
        theta = ga * i
        verts.append(ScanVertex(i,
            R        = 0.1 + 0.8*r_xy,
            lam      = theta % (2*math.pi),
            rotation = math.acos(max(-1.0, min(1.0, y)))))
    return verts


def scan_signal(signal, vertices, sr=SR):
    """信号を全頂点でスキャン → 応答ベクトル (60,)"""
    resp = np.zeros(len(vertices))
    for v in vertices:
        resp[v.idx] = phi_d_tsuji(signal, v.tau, v.f_c, v.bw, sr)
    return resp


# ════════════════════════════════════════════════
# 実験1: 時間フレーム分割スキャン
# ════════════════════════════════════════════════

def temporal_scan(signal, vertices, sr=SR,
                  frame_ms=200.0, hop_ms=100.0):
    """
    時間フレームごとにスキャンし、応答の時間変化を記録する。
    Returns: frames × vertices の行列
    """
    FRAME = int(sr * frame_ms / 1000)
    HOP   = int(sr * hop_ms   / 1000)
    N = len(signal)
    frames = []
    times  = []
    for start in range(0, N - FRAME, HOP):
        seg = signal[start:start+FRAME]
        resp = scan_signal(seg, vertices, sr)
        frames.append(resp)
        times.append((start + FRAME/2) / sr)
    return np.array(frames), np.array(times)


# ════════════════════════════════════════════════
# 実験2: 混合信号のみからの波源数推定
# ════════════════════════════════════════════════

def estimate_source_count(mixed_frames: np.ndarray,
                          max_sources: int = 5) -> dict:
    """
    混合信号のフレーム行列から波源数を推定する。

    方法: 時間フレーム × 頂点 の行列に対して SVD を適用。
    固有値の分布から「有効なランク」= 波源数を推定。

    根拠:
      各波源が固有の Φ_d パターンを持つなら、
      フレーム行列は「波源数」のランクを持つはず。
      → 特異値の急落点が波源数のヒント。
    """
    # フレーム数が少なすぎる場合
    if mixed_frames.shape[0] < 3:
        return {'n_sources': 1, 'singular_values': [], 'method': 'insufficient_frames'}

    # SVD
    U, S, Vt = np.linalg.svd(mixed_frames, full_matrices=False)

    # 特異値の「肘」を検出（2階差分で急落点を探す）
    S_norm = S / (S[0] + 1e-12)
    diffs = np.diff(S_norm)
    # 最大降下点
    elbow = int(np.argmin(diffs[:max_sources])) + 1

    # エネルギー比率で確認
    total_energy = np.sum(S**2)
    cumulative   = np.cumsum(S**2) / total_energy

    return {
        'n_sources':       elbow,
        'singular_values': S[:max_sources+2].tolist(),
        'S_norm':          S_norm[:max_sources+2].tolist(),
        'cumulative_energy': cumulative[:max_sources+2].tolist(),
        'method':          'svd_elbow',
    }


# ════════════════════════════════════════════════
# 実験3: 混合から単独を分解（識別力頂点自動選択）
# ════════════════════════════════════════════════

def decompose_mixture(mixed_resp: np.ndarray,
                      reference_resps: dict[str, np.ndarray],
                      top_k: int = 5) -> dict:
    """
    混合の応答ベクトルを単独信号の線形結合で分解。
    top_k: 識別力の高い頂点だけ使う。

    ※ 本来は reference_resps が「未知」だが、
      今は検証のために単独信号の応答を使う。
    """
    names = list(reference_resps.keys())
    vecs  = [reference_resps[n] for n in names]

    # 識別力 = 全ペア間の差の絶対値の総和
    disc = np.zeros(len(mixed_resp))
    from itertools import combinations
    for v1, v2 in combinations(vecs, 2):
        disc += np.abs(v1 - v2)

    top_idx = np.argsort(disc)[-top_k:]

    M_k = mixed_resp[top_idx]
    A_k = np.column_stack([v[top_idx] for v in vecs])

    coeffs, _, _, _ = np.linalg.lstsq(A_k, M_k, rcond=None)
    pred = A_k @ coeffs
    corr = np.corrcoef(M_k, pred)[0,1] if M_k.std()>1e-10 and pred.std()>1e-10 else 0

    return {
        'names':   names,
        'coeffs':  dict(zip(names, coeffs.tolist())),
        'corr':    float(corr),
        'top_idx': top_idx.tolist(),
        'top_k':   top_k,
    }


# ════════════════════════════════════════════════
# メイン
# ════════════════════════════════════════════════

if __name__ == '__main__':
    print("="*65)
    print("辻井式位相分解機 v2")
    print("時間フレーム × 頂点 × 波源数推定")
    print("="*65)

    np.random.seed(42)
    DUR = 2.0
    ONSETS = {'piano': 0.000, 'violin': 0.008, 'flute': 0.003}

    sources = {
        'piano':  make_piano( f0=220, sr=SR, dur=DUR, onset=ONSETS['piano']),
        'violin': make_violin(f0=440, sr=SR, dur=DUR, onset=ONSETS['violin']),
        'flute':  make_flute( f0=880, sr=SR, dur=DUR, onset=ONSETS['flute']),
    }
    mixed = sum(sources.values())
    mixed /= (max(abs(mixed)) + 1e-9)

    vertices = fullerene_vertices(60)

    # ── 実験1: 時間フレームスキャン ──
    print("\n[実験1] 時間フレーム分割スキャン")
    t0 = time.time()
    all_temporal = {}
    for name, sig in list(sources.items()) + [('mixed', mixed)]:
        frames, times = temporal_scan(sig, vertices, sr=SR,
                                      frame_ms=200.0, hop_ms=100.0)
        all_temporal[name] = (frames, times)
        print(f"  {name:<8}: {frames.shape[0]}フレーム × {frames.shape[1]}頂点  "
              f"std_mean={frames.std(axis=1).mean():.4f}")
    print(f"  時間: {time.time()-t0:.2f}s")

    # ── Φ_d の時間変化（piano: 減衰するはず）──
    print("\n[ピアノの Φ_d 時間変化（頂点22: 識別力最高）]")
    piano_frames, times = all_temporal['piano']
    print(f"  時刻(s)  Φ_d[v22]  平均Φ_d")
    for i, t_sec in enumerate(times):
        v22 = piano_frames[i, 22]
        mean = piano_frames[i].mean()
        print(f"  {t_sec:.2f}s    {v22:+.4f}   {mean:+.4f}")

    # ── 実験2: SVDによる波源数推定 ──
    print("\n[実験2] 混合信号のみからの波源数推定（SVD）")
    mixed_frames, _ = all_temporal['mixed']
    result = estimate_source_count(mixed_frames)
    print(f"  推定波源数: {result['n_sources']}")
    print(f"  特異値(正規化): {[f'{v:.4f}' for v in result['S_norm']]}")
    print(f"  累積エネルギー: {[f'{v:.3f}' for v in result['cumulative_energy']]}")

    # 各単独でも確認
    for name in ['piano', 'violin', 'flute']:
        frames_n, _ = all_temporal[name]
        r = estimate_source_count(frames_n)
        print(f"  {name:<8} 推定: {r['n_sources']}  "
              f"sv={[f'{v:.4f}' for v in r['S_norm']]}")

    # ── 実験3: 識別力頂点TOP-Kでの分解 ──
    print("\n[実験3] 識別力TOP-K頂点による混合分解")

    # 全体の応答ベクトル（中央200msセグメント）
    ref_resps = {n: scan_signal(sources[n], vertices) for n in sources}
    mixed_resp = scan_signal(mixed, vertices)

    print(f"  {'K':>4}  {'相関':>8}  {'係数(p,v,f)':>30}")
    for k in [3, 5, 8, 10, 15]:
        res = decompose_mixture(mixed_resp, ref_resps, top_k=k)
        coeffs = res['coeffs']
        print(f"  {k:>4}  {res['corr']:>8.4f}  "
              f"({coeffs['piano']:+.3f}, {coeffs['violin']:+.3f}, {coeffs['flute']:+.3f})")

    # ── 最良分解の詳細 ──
    print("\n[最良分解 K=5 の識別頂点]")
    res5 = decompose_mixture(mixed_resp, ref_resps, top_k=5)
    print(f"  使用頂点: {res5['top_idx']}")
    for idx in res5['top_idx']:
        v = vertices[idx]
        print(f"  頂点{idx:2d}: f_c={v.f_c:5.0f}Hz  τ={v.tau:3d}  BW={v.bw:5.0f}  "
              f"P={ref_resps['piano'][idx]:+.4f} "
              f"V={ref_resps['violin'][idx]:+.4f} "
              f"F={ref_resps['flute'][idx]:+.4f} "
              f"M={mixed_resp[idx]:+.4f}")

    # ── 最終判定 ──
    print("\n" + "="*65)
    print("総括")
    print("="*65)
    n_est = result['n_sources']
    best_k5 = decompose_mixture(mixed_resp, ref_resps, top_k=5)

    print(f"\nSVD波源数推定: {n_est} (真値=3)")
    print(f"K=5分解精度: 相関={best_k5['corr']:.4f}")
    print(f"係数比率（エネルギー寄与）:")
    total_c = sum(abs(v) for v in best_k5['coeffs'].values())
    for name, c in best_k5['coeffs'].items():
        bar = '█' * int(abs(c)/total_c * 30)
        print(f"  {name:<8}: {c:+.4f}  {bar}")

    # 保存
    np.save('/home/claude/scan_v3_results.npy', {
        'ref_resps': ref_resps,
        'mixed_resp': mixed_resp,
        'temporal': {k: v[0] for k, v in all_temporal.items()},
        'times': times,
    })
    print("\n[保存] scan_v3_results.npy")
