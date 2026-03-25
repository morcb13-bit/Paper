"""
pipeline_v2.py  ―― 引けるものを引いて残差を帰属する
=======================================================

設計思想:
  全部分けようとしない。
  強いものから順に引いて、引けなくなったら止める。
  残差は引き済み波源のエンベロープ相関で帰属させる。
  帰属できないものは残差として返す。

3フェーズ:
  Phase 1: pull_one() を繰り返す（改善量が閾値を下回ったら停止）
  Phase 2: 残差をフレームエンベロープ相関で既存波源に帰属
  Phase 3: 最終評価（残差は独立した成分として返す）

結果（強いビブラートあり合成音）:
  piano: 相関=0.877  SDR=3.4dB  ← 事前知識なし
  flute: 相関=0.846  SDR=3.4dB  ← 残差帰属で復元
  violin: 残差として返す（次セッションの課題）

  v7手動プロファイル比較:
    piano  0.846/3.3 → 0.877/3.4  超えた
    flute  0.459/0.6 → 0.846/3.4  大幅改善
"""

import numpy as np
import scipy.signal as ss
import math
from scipy.signal import find_peaks

SR  = 22050
PHI = (1 + math.sqrt(5)) / 2


# ════════════════════════════════════════════════
# プリミティブ
# ════════════════════════════════════════════════

def get_coeff(sig: np.ndarray, freq: float, sr: int = SR) -> complex:
    N = len(sig); t = np.arange(N) / sr
    return (2*np.dot(sig, np.cos(2*math.pi*freq*t))/N +
            1j*2*np.dot(sig, np.sin(2*math.pi*freq*t))/N)

def extract_sinusoid(N: int, freq: float, coeff: complex, sr: int = SR) -> np.ndarray:
    t = np.arange(N) / sr
    return coeff.real*np.cos(2*math.pi*freq*t) + coeff.imag*np.sin(2*math.pi*freq*t)

def amp_blackman(sig: np.ndarray, f: float, sr: int = SR) -> float:
    N = len(sig); t = np.arange(N)/sr; w = np.blackman(N); wg = np.mean(w)
    return math.sqrt(
        (2*np.dot(sig, np.cos(2*math.pi*f*t)*w)/(N*wg))**2 +
        (2*np.dot(sig, np.sin(2*math.pi*f*t)*w)/(N*wg))**2)

def fib_converge(sig: np.ndarray, f0: float, sr: int = SR,
                 n_iter: int = 40, window: float = 10.0) -> float:
    """黄金比探索でf0を精密化"""
    al = amp_blackman(sig, max(0.5, f0-0.1), sr)
    ar = amp_blackman(sig, f0+0.1, sr)
    lo, hi = (f0, f0+window) if ar > al else (max(0.5, f0-window), f0)
    cf = f0
    for _ in range(n_iter):
        if hi-lo < 0.05: break
        m1 = hi-(hi-lo)/PHI; m2 = lo+(hi-lo)/PHI
        if amp_blackman(sig, m1, sr) < amp_blackman(sig, m2, sr):
            lo = m1; cf = m2
        else:
            hi = m2; cf = m1
    return cf

def inst_freq_median(sig: np.ndarray, f_approx: float, sr: int = SR) -> float:
    """
    瞬時周波数の中央値でf0推定（ビブラート対応）
    ビブラートは対称に揺れるので中央値は真のf0に収束する
    """
    nyq = sr/2
    lo = max(20.0, f_approx*0.5); hi = min(nyq*0.98, f_approx*2.0)
    if hi > lo+10:
        b, a = ss.butter(2, [lo/nyq, hi/nyq], btype='band')
        filtered = ss.lfilter(b, a, sig)
    else:
        filtered = sig
    z = ss.hilbert(filtered)
    phase = np.unwrap(np.angle(z))
    inst_f = np.diff(phase) * sr / (2*math.pi)
    mask = (inst_f > f_approx*0.5) & (inst_f < f_approx*2.0)
    if mask.sum() < 10:
        return f_approx
    return float(np.median(inst_f[mask]))

def sdr_db(ref: np.ndarray, est: np.ndarray) -> float:
    noise = ref - est
    return 10*math.log10((np.mean(ref**2)+1e-12) / (np.mean(noise**2)+1e-12))


# ════════════════════════════════════════════════
# Phase 1: 1波源を検出して引く
# ════════════════════════════════════════════════

def pull_one(residual: np.ndarray, sr: int = SR):
    """
    残差から最も強い波源を1つ検出・抽出する。

    処理:
      1. FFT最大振幅ピーク → f_cand
      2. fib_converge と 瞬時周波数中央値 で f0 を両推定
         → 振幅が大きい方を採用（ビブラートに強い）
      3. テンプレート構築（先頭50msで倍音係数）
      4. フレームベース serial subtraction

    Returns:
      (f0, extracted_signal), new_residual
      検出できない場合は None, residual
    """
    N = len(residual)
    WIN   = int(sr * 0.05)   # テンプレート構築用 50ms
    FRAME = int(sr * 0.05)
    HOP   = int(sr * 0.025)

    # FFT最大振幅ピーク
    spec = np.abs(np.fft.rfft(residual * np.hanning(N)))
    freqs = np.fft.rfftfreq(N, 1/sr)
    peaks_idx, _ = find_peaks(spec, height=spec.max()*0.03,
                               distance=int(20/(sr/N)))
    if len(peaks_idx) == 0:
        return None, residual

    f_cand = float(freqs[peaks_idx[np.argmax(spec[peaks_idx])]])
    if not (50 < f_cand < 4000):
        return None, residual

    # f0 推定: fib_converge vs 瞬時周波数中央値
    f0_fib  = fib_converge(residual[:WIN], f_cand, sr, window=f_cand*0.15)
    f0_inst = inst_freq_median(residual, f_cand, sr)
    f0 = f0_fib if (amp_blackman(residual[:WIN], f0_fib, sr) >=
                    amp_blackman(residual[:WIN], f0_inst, sr)) else f0_inst

    # テンプレート構築
    c_f0 = get_coeff(residual[:WIN], f0, sr)
    if abs(c_f0) < 1e-6:
        return None, residual

    tmpl = {}
    for h in range(1, 8):
        freq = f0 * h
        if freq > sr/2: break
        tmpl[h] = get_coeff(residual[:WIN], freq, sr) / (c_f0 + 1e-12)

    # フレームベース抽出
    out = np.zeros(N); ov = np.zeros(N); hann = np.hanning(FRAME)
    for start in range(0, N-FRAME, HOP):
        seg = residual[start:start+FRAME].copy()
        c_f0_ = get_coeff(seg, f0, sr)
        if abs(c_f0_) < 1e-8: continue
        for h, ratio in tmpl.items():
            freq = f0 * h
            if freq > sr/2: break
            c_use = 0.6*c_f0_*ratio + 0.4*get_coeff(seg, freq, sr)
            if abs(c_use) < 1e-7: continue
            comp = extract_sinusoid(FRAME, freq, c_use, sr)
            out[start:start+FRAME] += comp * hann
            seg -= comp
        ov[start:start+FRAME] += hann

    ov = np.maximum(ov, 1e-9)
    out /= ov
    return (f0, out), residual - out


# ════════════════════════════════════════════════
# Phase 2: 残差をエンベロープ相関で帰属
# ════════════════════════════════════════════════

def rms_envelope(sig: np.ndarray, starts: list, frame_len: int) -> np.ndarray:
    return np.array([math.sqrt(np.mean(sig[s:s+frame_len]**2)+1e-12)
                     for s in starts])

def assign_residual(residual: np.ndarray,
                    extracted: list,
                    sr: int = SR) -> list:
    """
    残差のフレームエンベロープと各抽出波源のエンベロープを比較し、
    最も相関が高い波源に残差を統合する。

    相関が閾値（0.1）未満なら統合せずそのまま返す。
    """
    if not extracted or residual.std() < 1e-6:
        return extracted

    N = len(residual)
    FRAME = int(sr * 0.05); HOP = int(sr * 0.025)
    starts = list(range(0, N-FRAME, HOP))

    res_env = rms_envelope(residual, starts, FRAME)
    best_idx, best_corr = -1, 0.0

    for i, (f0, sig) in enumerate(extracted):
        ext_env = rms_envelope(sig, starts, FRAME)
        if res_env.std() < 1e-10 or ext_env.std() < 1e-10:
            continue
        c = float(np.corrcoef(res_env, ext_env)[0, 1])
        if abs(c) > abs(best_corr):
            best_corr = c; best_idx = i

    if best_idx >= 0 and abs(best_corr) > 0.1:
        f0, sig = extracted[best_idx]
        extracted[best_idx] = (f0, sig + residual)
        print(f"  残差 → f0={f0:.0f}Hz に帰属（エンベロープ相関={best_corr:+.4f}）")
    else:
        print(f"  残差は未帰属（最大相関={best_corr:.4f} < 0.1）")

    return extracted


# ════════════════════════════════════════════════
# メイン: 3フェーズ分離
# ════════════════════════════════════════════════

def separate(mixed: np.ndarray,
             sr: int = SR,
             max_sources: int = 8,
             rms_threshold: float = 0.003,
             verbose: bool = True) -> dict:
    """
    混合信号から波源を分離する。

    Parameters
    ----------
    mixed         : 混合信号
    sr            : サンプルレート
    max_sources   : 最大波源数（上限、実際は残差改善で自動停止）
    rms_threshold : この改善量を下回ったら停止
    verbose       : 進捗表示

    Returns
    -------
    {
      'sources': [(f0, signal), ...],  # 抽出できた波源
      'residual': np.ndarray,           # 分けられなかった残差
    }
    """
    residual = mixed.copy()
    extracted = []

    # Phase 1: 引けるだけ引く
    if verbose: print("[Phase 1] 波源を順番に引く")
    for i in range(max_sources):
        rms_before = math.sqrt(np.mean(residual**2))
        result, residual_new = pull_one(residual, sr)

        if result is None:
            if verbose: print(f"  [{i+1}] 検出なし → 停止")
            break

        f0, out = result
        rms_after = math.sqrt(np.mean(residual_new**2))
        improvement = rms_before - rms_after

        if improvement < rms_threshold:
            if verbose:
                print(f"  [{i+1}] f0={f0:.0f}Hz  改善={improvement:.4f} < {rms_threshold} → 停止")
            break

        residual = residual_new
        extracted.append((f0, out))

        if verbose:
            print(f"  [{i+1}] f0={f0:.0f}Hz  "
                  f"残差RMS: {rms_before:.4f}→{rms_after:.4f}  "
                  f"改善={improvement:.4f}")

    if verbose:
        print(f"  → {len(extracted)}波源抽出  残差RMS={math.sqrt(np.mean(residual**2)):.4f}")

    # Phase 2: 残差を帰属
    if verbose: print("\n[Phase 2] 残差を既存波源に帰属")
    extracted = assign_residual(residual, extracted, sr)

    # 帰属後の残差
    final_residual = residual.copy()
    if extracted:
        # 帰属済みなら残差はゼロ化（assign_residual内で統合済み）
        # 未帰属の場合はそのまま
        all_assigned = any(abs(np.corrcoef(
            rms_envelope(residual, list(range(0,len(residual)-int(sr*0.05),int(sr*0.025))), int(sr*0.05)),
            rms_envelope(sig, list(range(0,len(sig)-int(sr*0.05),int(sr*0.025))), int(sr*0.05))
        )[0,1]) > 0.1 for _, sig in extracted)
        if all_assigned:
            final_residual = np.zeros(len(mixed))

    return {
        'sources': extracted,
        'residual': final_residual,
    }


# ════════════════════════════════════════════════
# テスト音源生成
# ════════════════════════════════════════════════

def _adsr(N, sr, attack=0.05, decay=0.1, sustain=0.7, release=0.2):
    env = np.ones(N)*sustain
    a_s, d_s, r_s = int(attack*sr), int(decay*sr), int(release*sr)
    env[:a_s] = np.linspace(0, 1, a_s)
    env[a_s:a_s+d_s] = np.linspace(1, sustain, d_s)
    env[-r_s:] = np.linspace(sustain, 0, r_s)
    return env

def make_piano(f0=220, sr=SR, dur=2.0, onset=0.0):
    N = int(sr*dur); t = np.arange(N)/sr; B = 0.0001; s = np.zeros(N)
    for h, amp in [(1,1.0),(2,0.6),(3,0.4),(4,0.2),(5,0.12),(6,0.07)]:
        fh = f0*h*math.sqrt(1+B*h**2)
        s += amp*np.exp(-t/(2.0/h))*np.sin(2*math.pi*fh*t)
    s += 0.02*np.random.randn(N)*np.exp(-t/0.02)
    shift = int(onset*sr)
    s = np.concatenate([np.zeros(shift), s])[:N]
    return s / (max(abs(s))+1e-9)

def make_violin(f0=440, sr=SR, dur=2.0, onset=0.0, vib_depth=0.010):
    N = int(sr*dur); t = np.arange(N)/sr
    vib = vib_depth*(1-np.exp(-t/0.1))*np.sin(2*math.pi*5.5*t)
    s = np.zeros(N)
    for h, amp in [(1,1.0),(2,0.4),(3,0.7),(4,0.3),(5,0.5),(6,0.2),(7,0.35)]:
        s += amp*np.sin(2*math.pi*f0*h*(1+vib)*t)
    env = _adsr(N, sr, 0.08, 0.05, 0.85, 0.3)
    s = s*env + 0.03*np.random.randn(N)*env
    shift = int(onset*sr)
    s = np.concatenate([np.zeros(shift), s])[:N]
    return s / (max(abs(s))+1e-9)

def make_flute(f0=880, sr=SR, dur=2.0, onset=0.0, vib_depth=0.008):
    N = int(sr*dur); t = np.arange(N)/sr
    vib = vib_depth*(1-np.exp(-t/0.2))*np.sin(2*math.pi*6.0*t)
    s = np.zeros(N)
    for h, amp in [(1,1.0),(2,0.15),(3,0.05),(4,0.02)]:
        s += amp*np.sin(2*math.pi*f0*h*(1+vib)*t)
    air = 0.04*np.random.randn(N)
    b, a = ss.butter(4, [f0/sr*1.5, min(f0/sr*8, 0.9)], btype='band')
    s += ss.lfilter(b, a, air)
    env = _adsr(N, sr, 0.05, 0.02, 0.9, 0.25)
    s = s*env
    shift = int(onset*sr)
    s = np.concatenate([np.zeros(shift), s])[:N]
    return s / (max(abs(s))+1e-9)


# ════════════════════════════════════════════════
# メイン実行
# ════════════════════════════════════════════════

if __name__ == '__main__':
    print("="*60)
    print("pipeline_v2: 引けるものを引いて残差を帰属する")
    print("="*60)

    np.random.seed(42)
    DUR = 2.0
    ONSETS = {'piano': 0.000, 'violin': 0.008, 'flute': 0.003}

    # 強いビブラートあり（実演奏に近い）
    sources = {
        'piano':  make_piano( f0=220, sr=SR, dur=DUR, onset=ONSETS['piano']),
        'violin': make_violin(f0=440, sr=SR, dur=DUR, onset=ONSETS['violin'], vib_depth=0.010),
        'flute':  make_flute( f0=880, sr=SR, dur=DUR, onset=ONSETS['flute'],  vib_depth=0.008),
    }
    mixed = sum(sources.values())
    mixed /= (max(abs(mixed)) + 1e-9)

    print(f"\nビブラート深度: violin=±1.0%  flute=±0.8%")
    print(f"onset: piano=0ms  violin=8ms  flute=3ms\n")

    # 分離
    result = separate(mixed, SR)

    # 評価
    print("\n[評価]")
    print(f"  {'推定f0':>8}  {'最良一致':>8}  {'相関':>8}  {'SDR':>8}")
    print("  " + "-"*40)

    for f0, sig in result['sources']:
        best_n, best_c, best_s = '', 0.0, -99.0
        for name, ref in sources.items():
            c = float(np.corrcoef(ref, sig)[0,1]) if sig.std() > 1e-10 else 0.0
            if c > best_c:
                best_c = c; best_s = sdr_db(ref, sig); best_n = name
        print(f"  {f0:>8.0f}Hz  {best_n:>8}  {best_c:>8.4f}  {best_s:>8.1f}dB")

    res = result['residual']
    if res.std() > 1e-4:
        rms_r = math.sqrt(np.mean(res**2))
        rms_i = math.sqrt(np.mean(mixed**2))
        print(f"\n  残差RMS: {rms_r:.4f} (入力比{rms_r/rms_i*100:.1f}%)")
        print("  残差の内容:")
        for name, ref in sources.items():
            c = float(np.corrcoef(ref, res)[0,1]) if res.std() > 1e-10 else 0.0
            print(f"    vs {name}: 相関={c:+.4f}")

    print()
    print("[参考] v7 手動プロファイル: piano=0.846/3.3  violin=0.220/0.2  flute=0.459/0.6")
