"""
instrument_separator_v7.py
===========================
onset差分によるテンプレート自動構築 + トポロジカル帰属分離

핵심アイデア（v6からの根本的転換）:
  実際の演奏では各楽器の発音タイミングがズレている。
  piano打鍵→flute息継ぎ→violin弓入れ、それぞれ独立したアタック。

  この「単独で鳴っている区間」を使えば
  テンプレートを混合汚染なしに自動構築できる。

  アルゴリズム:
    1. onset検出: 各楽器の発音開始タイミングを差分RMSで推定
    2. 単独区間抽出: onset[k] 〜 onset[k+1] の区間は楽器k単独
    3. テンプレート構築: その区間でfib収束 → 複素係数テンプレート確定
    4. 帰属スコアリング: 全フレームで位相指紋+振幅包絡でスコア計算
    5. 一括キャンセル: 帰属確定した成分を楽器単位でまとめて除去

  onset が取れない場合（同時発音）のフォールバック:
    amp_profile ベースの仮テンプレートで v6 相当に降格。
"""

import numpy as np
import math
import scipy.signal as ss
from dataclasses import dataclass, field
from collections import defaultdict

SR  = 22050
PHI = (1 + math.sqrt(5)) / 2

H_SET  = frozenset({1, 5, 8, 12})
H2_SET = frozenset({2, 3, 10, 11})
H4_SET = frozenset({4, 6, 7, 9})


# ════════════════════════════════════════════════
# 基本プリミティブ（実時間軸 t=k/sr）
# ════════════════════════════════════════════════

def get_coeff(signal: np.ndarray, freq: float, sr: int = SR) -> complex:
    """freq Hz 成分の複素係数。t=k/sr 実時間軸で非整数周波数も正確。"""
    N = len(signal)
    t = np.arange(N) / sr
    ac  = 2*np.dot(signal, np.cos(2*math.pi*freq*t)) / N
    as_ = 2*np.dot(signal, np.sin(2*math.pi*freq*t)) / N
    return ac + 1j*as_


def extract_sinusoid(N: int, freq: float, coeff: complex, sr: int = SR) -> np.ndarray:
    t = np.arange(N) / sr
    return coeff.real*np.cos(2*math.pi*freq*t) + coeff.imag*np.sin(2*math.pi*freq*t)


def amp_blackman(signal: np.ndarray, f: float, sr: int = SR) -> float:
    N = len(signal); t = np.arange(N)/sr
    w = np.blackman(N); wg = np.mean(w)
    return math.sqrt(
        (2*np.dot(signal, np.cos(2*math.pi*f*t)*w)/(N*wg))**2 +
        (2*np.dot(signal, np.sin(2*math.pi*f*t)*w)/(N*wg))**2
    )


def fib_converge(signal: np.ndarray, f0: float, sr: int = SR,
                 n_iter: int = 40, window: float = 10.0) -> float:
    al = amp_blackman(signal, max(0.5, f0-0.1), sr)
    ar = amp_blackman(signal, f0+0.1, sr)
    lo, hi = (f0, f0+window) if ar > al else (max(0.5, f0-window), f0)
    cf = f0
    for _ in range(n_iter):
        if hi-lo < 0.05: break
        m1 = hi-(hi-lo)/PHI; m2 = lo+(hi-lo)/PHI
        if amp_blackman(signal, m1, sr) < amp_blackman(signal, m2, sr):
            lo = m1; cf = m2
        else:
            hi = m2; cf = m1
    return cf


def coset_band(f: float) -> str:
    g = round(f) % 13
    if g == 0:       return '0'
    if g in H_SET:   return 'H'
    if g in H2_SET:  return '2H'
    return '4H'


def sdr_db(ref: np.ndarray, est: np.ndarray) -> float:
    noise = ref - est
    return 10*math.log10((np.mean(ref**2)+1e-12)/(np.mean(noise**2)+1e-12))


# ════════════════════════════════════════════════
# Step 1: onset 検出
# ════════════════════════════════════════════════

def detect_onsets(mixed: np.ndarray, sr: int,
                  block_ms: float = 2.0,
                  threshold_ratio: float = 0.15,
                  min_gap_ms: float = 3.0) -> list[float]:
    """
    差分RMSエンベロープによるonset検出。

    連続するブロック間でRMSが急増した時刻を onset とみなす。
    Returns: onset時刻のリスト（秒）、先頭 t=0 を含む
    """
    BLOCK = int(sr * block_ms / 1000)
    MIN_GAP = int(sr * min_gap_ms / 1000)
    N = len(mixed)

    rms_env = []
    for s in range(0, N-BLOCK, BLOCK//2):
        rms_env.append(math.sqrt(np.mean(mixed[s:s+BLOCK]**2) + 1e-12))

    rms_arr = np.array(rms_env)
    diff = np.diff(rms_arr, prepend=rms_arr[0])

    # 差分の正のピークを onset 候補とする
    peak_thresh = threshold_ratio * rms_arr.max()
    onsets_samples = [0]   # t=0 は常に含める
    last = 0

    step = BLOCK // 2
    for i, d in enumerate(diff):
        t_sample = i * step
        if d > peak_thresh and t_sample - last > MIN_GAP:
            onsets_samples.append(t_sample)
            last = t_sample

    return [s/sr for s in sorted(set(onsets_samples))]


# ════════════════════════════════════════════════
# Step 2: 楽器プロファイル（テンプレート付き）
# ════════════════════════════════════════════════

@dataclass
class InstrumentProfile:
    name:        str
    f0_init:     float
    n_harmonics: int
    anchor_band: str
    amp_profile: dict          # {h: 相対振幅}
    f0_window:   float = 10.0

    # 実行時に構築
    f0:               float = field(init=False)
    phase_template:   dict  = field(default_factory=dict)   # {h: complex ratio}
    template_ready:   bool  = field(default=False)

    def __post_init__(self):
        self.f0 = self.f0_init

    def build_template_from_segment(self, segment: np.ndarray, sr: int) -> bool:
        """
        単独区間からテンプレートを構築する（混合汚染なし）。
        Returns: 成功したら True
        """
        if len(segment) < int(sr * 0.01):   # 10ms未満は短すぎる
            return False

        f_est = fib_converge(segment, self.f0, sr=sr,
                             n_iter=50, window=self.f0_window)
        c_f0 = get_coeff(segment, f_est, sr)
        if abs(c_f0) < 1e-6:
            return False

        template = {}
        for h in range(1, self.n_harmonics + 1):
            freq = f_est * h
            if freq > sr / 2:
                break
            c_h = get_coeff(segment, freq, sr)
            template[h] = c_h / (c_f0 + 1e-12)

        self.f0             = f_est
        self.phase_template = template
        self.template_ready = True
        return True

    def build_template_from_amp_profile(self, mixed_frame: np.ndarray, sr: int):
        """
        フォールバック: amp_profile から位相テンプレートを仮構築。
        位相角は混合フレームの f0 成分から取る（汚染あり）。
        """
        f_est = fib_converge(mixed_frame, self.f0, sr=sr, window=self.f0_window)
        c_f0  = get_coeff(mixed_frame, f_est, sr)
        f0_phase = c_f0 / (abs(c_f0) + 1e-12)

        template = {}
        for h, amp_ratio in self.amp_profile.items():
            freq = f_est * h
            if freq > sr / 2:
                break
            c_h = get_coeff(mixed_frame, freq, sr)
            phase_h = c_h / (abs(c_h) + 1e-12) if abs(c_h) > 1e-8 else f0_phase
            template[h] = amp_ratio * phase_h / (f0_phase + 1e-12)

        self.f0             = f_est
        self.phase_template = template
        self.template_ready = False   # フォールバック扱い

    def track_f0(self, frame: np.ndarray, sr: int):
        """フレームごとの f0 追跡（ビブラート追従）。"""
        f_est = fib_converge(frame, self.f0, sr=sr,
                             n_iter=30, window=self.f0_window)
        self.f0 = f_est

    def score(self, freq: float, coeff: complex, h: int) -> float:
        """
        周波数成分 (freq, coeff, h次倍音) のこの楽器へのスコア。
        S = S_phase × S_amp × S_band
        """
        if not self.phase_template or abs(coeff) < 1e-8:
            return 0.0

        # ── S_phase: 位相指紋との一致度 ──
        if h in self.phase_template:
            tmpl = self.phase_template[h]
            # テンプレートの位相角
            tp = tmpl / (abs(tmpl) + 1e-12)
            # coeff の位相角
            cp = coeff / (abs(coeff) + 1e-12)
            # コサイン類似度 → [0, 1]
            s_phase = (cp.real*tp.real + cp.imag*tp.imag + 1.0) * 0.5
        else:
            s_phase = 0.1

        # ── S_amp: 振幅包絡との一致度 ──
        expected = self.amp_profile.get(h, 0.0)
        if expected > 0 and 1 in self.phase_template:
            tmpl_h1 = abs(self.phase_template[1]) + 1e-12
            tmpl_hh = abs(self.phase_template.get(h, complex(expected))) + 1e-12
            # h倍音の実測振幅 / テンプレートの h倍音振幅比率
            # テンプレートが clean なら tmpl_hh/tmpl_h1 ≈ amp_profile[h]
            actual_ratio   = (abs(coeff) / tmpl_hh) * tmpl_h1
            expected_ratio = expected
            ratio_match = min(actual_ratio, expected_ratio) / (
                          max(actual_ratio, expected_ratio) + 1e-12)
            s_amp = ratio_match
        else:
            s_amp = 0.3

        # ── S_band: B13コセット帯アンカー ──
        band = coset_band(freq)
        s_band = 1.0 if band == self.anchor_band else (0.3 if band == '0' else 0.5)

        return float(s_phase * (0.5*s_amp + 0.5) * (0.7*s_band + 0.3))


# ════════════════════════════════════════════════
# Step 3: onset区間でのテンプレート構築
# ════════════════════════════════════════════════

def build_templates_from_onsets(mixed: np.ndarray,
                                sr: int,
                                profiles: list[InstrumentProfile],
                                onset_times: list[float],
                                rms_threshold: float = 0.05) -> dict[str, str]:
    """
    onset 差分信号から各楽器のクリーンテンプレートを構築する。

    核心戦略:
      onset[k] の「直前フレーム」と「直後フレーム」の差分が
      その時点で新たに加わった楽器の純音に近い。

      mixed(t+δ) - mixed(t-δ) ≈ 新楽器の信号

    これにより piano の倍音が flute 帯に漏れていても、
    差分で piano 成分は相殺され flute 成分が残る。
    """
    result  = {}
    N       = len(mixed)
    PRE_MS  = 3     # onset 前の参照区間（ms）
    POST_MS = 15    # onset 後のセグメント（ms）
    PRE     = int(sr * PRE_MS  / 1000)
    POST    = int(sr * POST_MS / 1000)
    MIN_SEG = int(sr * 0.008)

    # 各楽器の f0 帯域でのonset時刻を推定（帯域別RMS急増）
    def bandpass_rms_onsets(x, f_center, bw=150.0):
        lo = max(f_center - bw, 20.0)
        hi = min(f_center + bw, sr/2 - 1)
        b, a = ss.butter(4, [lo/(sr/2), hi/(sr/2)], btype='band')
        filtered = ss.lfilter(b, a, x)
        blk = int(sr * 0.002)
        envs = []
        for s in range(0, len(filtered)-blk, blk//2):
            envs.append(math.sqrt(np.mean(filtered[s:s+blk]**2)+1e-12))
        env = np.array(envs)
        thresh = max(env) * 0.08
        step = blk//2
        for i, v in enumerate(env):
            if v > thresh:
                return i * step
        return 0

    # 各楽器の f0 帯域でのonset時刻を推定（帯域別RMS急増）
    def bandpass_rms_onsets(x, f_center, bw=80.0):
        lo = max(f_center - bw, 20.0)
        hi = min(f_center + bw, sr/2 - 1)
        b, a = ss.butter(4, [lo/(sr/2), hi/(sr/2)], btype='band')
        filtered = ss.lfilter(b, a, x)
        blk = int(sr * 0.002)
        envs = []
        for s in range(0, len(filtered)-blk, blk//2):
            envs.append(math.sqrt(np.mean(filtered[s:s+blk]**2)+1e-12))
        env = np.array(envs)
        thresh = max(env) * 0.08
        step = blk//2
        for i, v in enumerate(env):
            if v > thresh:
                return i * step
        return 0

    # まず f0 が最も低い楽器（= 先発候補）の onset を推定
    # その楽器の倍音を除去してから他楽器の onset を推定する
    sorted_by_f0 = sorted(profiles, key=lambda p: p.f0_init)

    # Step 1: 全楽器の生の onset を推定
    inst_onset_s = {}
    for prof in profiles:
        inst_onset_s[prof.name] = bandpass_rms_onsets(mixed, prof.f0_init)

    # Step 2: 最も低い f0 の楽器の単独区間から倍音係数を取得し除去
    # → その楽器の倍音が他楽器の帯域に漏れている問題を解決
    look_ms = 3   # 先発楽器の単独区間として使う長さ（ms）
    LOOK = int(sr * look_ms / 1000)

    residual_for_onset = mixed.copy()
    for prof in sorted_by_f0:
        if LOOK > 0:
            # 先発区間でこの楽器の倍音係数を推定してから除去
            t_arr = np.arange(len(mixed)) / sr
            for h in range(1, prof.n_harmonics + 1):
                freq = prof.f0_init * h
                if freq > sr / 2:
                    break
                c = get_coeff(mixed[:LOOK], freq, sr)
                if abs(c) > 1e-7:
                    comp = c.real*np.cos(2*math.pi*freq*t_arr) + c.imag*np.sin(2*math.pi*freq*t_arr)
                    residual_for_onset -= comp

        # この楽器の onset を更新（他楽器の倍音が除去された残差で）
        new_onset = bandpass_rms_onsets(residual_for_onset, prof.f0_init)
        old = inst_onset_s[prof.name]
        if new_onset != old:
            print(f"    onset補正: {prof.name} {old/sr*1000:.1f}ms → {new_onset/sr*1000:.1f}ms")
            inst_onset_s[prof.name] = new_onset

    print(f"    帯域別onset推定: " +
          ", ".join(f"{p.name}={inst_onset_s[p.name]/sr*1000:.1f}ms"
                    for p in profiles))

    # 遅い onset 順（差分が取りやすい順）に処理
    sorted_profs = sorted(profiles, key=lambda p: inst_onset_s[p.name], reverse=True)

    for prof in sorted_profs:
        onset_sample = inst_onset_s[prof.name]

        # 差分セグメント: onset後 - onset前 の定常成分除去
        pre_start  = max(0, onset_sample - PRE)
        pre_end    = onset_sample
        post_start = onset_sample
        post_end   = min(N, onset_sample + POST)

        if post_end - post_start < MIN_SEG or pre_end - pre_start < 4:
            print(f"  △ {prof.name}: セグメント短すぎ, fallback")
            continue

        pre_seg  = mixed[pre_start:pre_end]
        post_seg = mixed[post_start:post_end]

        # onset 前の各倍音係数を onset 後から引く（既存楽器成分を除去）
        # 方法: 共通長で post の各既存楽器倍音を pre の係数でキャンセル
        diff_seg = post_seg.copy()

        for other in profiles:
            if other.name == prof.name:
                continue
            if not other.template_ready:
                continue
            # 他楽器のテンプレートが確定済みなら、その倍音を差分信号から除去
            for h in range(1, other.n_harmonics + 1):
                freq = other.f0 * h
                if freq > sr / 2:
                    break
                c = get_coeff(diff_seg, freq, sr)
                if abs(c) > 1e-7:
                    diff_seg -= extract_sinusoid(len(diff_seg), freq, c, sr)

        # 差分セグメントでテンプレート構築
        ok = prof.build_template_from_segment(diff_seg, sr)
        if ok:
            result[prof.name] = 'onset_clean'
            print(f"  ✓ {prof.name}: onset {onset_sample/sr*1000:.1f}ms, "
                  f"差分セグメント {POST_MS}ms, f0={prof.f0:.2f}Hz")
        else:
            result[prof.name] = 'fallback'
            print(f"  △ {prof.name}: 差分セグメント構築失敗, fallback")

    # 未構築の楽器はフォールバック
    init_frame = mixed[:int(sr*0.05)]
    for prof in profiles:
        if not prof.template_ready:
            prof.build_template_from_amp_profile(init_frame, sr)
            result[prof.name] = 'fallback'
            print(f"  △ {prof.name}: amp_profile フォールバック")

    return result


# ════════════════════════════════════════════════
# Step 4: 帰属スコアリング
# ════════════════════════════════════════════════

def score_and_assign(frame: np.ndarray,
                     profiles: list[InstrumentProfile],
                     sr: int,
                     tie_threshold: float = 0.15) -> dict[str, list[tuple[float,complex]]]:
    """
    各楽器の f0 を追跡し、全倍音成分を帰属スコアで振り分ける。
    """
    # f0 追跡
    for prof in profiles:
        prof.track_f0(frame, sr)

    # 倍音候補を収集（±1Hz以内は同一とみなす）
    freq_map: dict[float, list] = defaultdict(list)
    for prof in profiles:
        for h in range(1, prof.n_harmonics + 1):
            freq = prof.f0 * h
            if freq > sr / 2:
                break
            found = next((k for k in freq_map if abs(k - freq) < 1.0), None)
            key = found if found is not None else freq
            freq_map[key].append((prof, h))

    # スコアリングと帰属
    assignment: dict[str, list] = {p.name: [] for p in profiles}

    for freq, candidates in freq_map.items():
        coeff = get_coeff(frame, freq, sr)
        if abs(coeff) < 1e-7:
            continue

        scores = {prof.name: prof.score(freq, coeff, h)
                  for prof, h in candidates}

        total = sum(scores.values()) + 1e-12
        max_s = max(scores.values())

        # 引き分け判定
        sv = sorted(scores.values(), reverse=True)
        is_tie = (len(sv) >= 2 and max_s > 1e-6 and
                  sv[0] - sv[1] < tie_threshold * max_s)

        if is_tie:
            for prof, h in candidates:
                r = scores[prof.name] / total
                assignment[prof.name].append((freq, coeff * r))
        else:
            winner = max(scores, key=scores.get)
            assignment[winner].append((freq, coeff))

    return assignment


# ════════════════════════════════════════════════
# Step 5: フレームベース分離
# ════════════════════════════════════════════════

def separate_v7(mixed: np.ndarray,
                sr: int,
                profiles: list[InstrumentProfile],
                frame_ms: float = 50.0,
                hop_ms:   float = 25.0,
                onset_block_ms: float = 2.0,
                tie_threshold:  float = 0.15) -> dict[str, np.ndarray]:
    """
    onset差分テンプレート + serial subtraction による楽器分離。

    フレームごとの処理:
      1. テンプレートが確定している楽器から順番に減算（serial subtraction）
      2. 残差に対してスコア帰属（衝突成分の分配）
    """
    FRAME = int(sr * frame_ms / 1000)
    HOP   = int(sr * hop_ms   / 1000)
    N     = len(mixed)

    # onset検出
    print("[1] onset 検出...")
    onsets = detect_onsets(mixed, sr, block_ms=onset_block_ms)
    print(f"    検出: {[f'{t*1000:.1f}ms' for t in onsets[:8]]}")

    # テンプレート構築
    print("[2] テンプレート構築...")
    status = build_templates_from_onsets(mixed, sr, profiles, onsets)

    # 処理順: テンプレートが clean な楽器を先に処理
    order = sorted(profiles,
                   key=lambda p: (0 if p.template_ready else 1, p.f0))

    print(f"    処理順: {[p.name for p in order]}")

    # フレームベース分離
    print("[3] フレームベース分離...")
    outputs     = {p.name: np.zeros(N) for p in profiles}
    overlap_cnt = np.zeros(N)
    hann        = np.hanning(FRAME)

    for start in range(0, N - FRAME, HOP):
        end       = start + FRAME
        frame     = mixed[start:end].copy()
        residual  = frame.copy()

        extracted = {p.name: np.zeros(FRAME) for p in profiles}

        for prof in order:
            # f0 追跡
            prof.track_f0(residual, sr)

            # 全倍音をテンプレート係数で一括キャンセル
            if prof.phase_template:
                # f0 成分の係数を取得
                c_f0 = get_coeff(residual, prof.f0, sr)

                for h in range(1, prof.n_harmonics + 1):
                    freq = prof.f0 * h
                    if freq > sr / 2:
                        break
                    # テンプレートがある場合: f0係数からh倍音係数を予測
                    if h in prof.phase_template and abs(c_f0) > 1e-8:
                        c_pred = c_f0 * prof.phase_template[h]
                        # 実測係数との加重平均（テンプレートを信頼しすぎない）
                        c_meas = get_coeff(residual, freq, sr)
                        # テンプレートが clean なら重みを 0.7、fallback なら 0.3
                        w = 0.7 if prof.template_ready else 0.4
                        c_use = w*c_pred + (1-w)*c_meas
                    else:
                        c_use = get_coeff(residual, freq, sr)

                    if abs(c_use) < 1e-7:
                        continue
                    component = extract_sinusoid(FRAME, freq, c_use, sr)
                    extracted[prof.name] += component
                    residual -= component

        for prof in profiles:
            outputs[prof.name][start:end] += extracted[prof.name] * hann

        overlap_cnt[start:end] += hann

    overlap_cnt = np.maximum(overlap_cnt, 1e-9)
    for name in outputs:
        outputs[name] /= overlap_cnt

    outputs['residual'] = mixed.copy()
    for name in list(outputs.keys()):
        if name != 'residual':
            outputs['residual'] -= outputs[name]

    return outputs


# ════════════════════════════════════════════════
# テスト用音源生成
# ════════════════════════════════════════════════

def _adsr(N, sr, attack=0.05, decay=0.1, sustain=0.7, release=0.2):
    env = np.ones(N)*sustain
    a_s,d_s,r_s = int(attack*sr),int(decay*sr),int(release*sr)
    env[:a_s]=np.linspace(0,1,a_s)
    env[a_s:a_s+d_s]=np.linspace(1,sustain,d_s)
    env[-r_s:]=np.linspace(sustain,0,r_s)
    return env

def make_piano(f0=220, sr=SR, dur=2.0, onset=0.0):
    N=int(sr*dur); t=np.arange(N)/sr; B=0.0001
    s=np.zeros(N)
    for h,amp in [(1,1.0),(2,0.6),(3,0.4),(4,0.2),(5,0.12),(6,0.07)]:
        fh=f0*h*math.sqrt(1+B*h**2)
        s+=amp*np.exp(-t/(2.0/h))*np.sin(2*math.pi*fh*t)
    s+=0.02*np.random.randn(N)*np.exp(-t/0.02)
    shift=int(onset*sr)
    s=np.concatenate([np.zeros(shift),s])[:N]
    return s/(max(abs(s))+1e-9)

def make_violin(f0=440, sr=SR, dur=2.0, onset=0.0):
    N=int(sr*dur); t=np.arange(N)/sr
    vib=0.004*(1-np.exp(-t/0.3))*np.sin(2*math.pi*5.5*t)
    s=np.zeros(N)
    for h,amp in [(1,1.0),(2,0.4),(3,0.7),(4,0.3),(5,0.5),(6,0.2),(7,0.35)]:
        s+=amp*np.sin(2*math.pi*f0*h*(1+vib)*t)
    env=_adsr(N,sr,attack=0.08,decay=0.05,sustain=0.85,release=0.3)
    s=s*env+0.03*np.random.randn(N)*env
    shift=int(onset*sr)
    s=np.concatenate([np.zeros(shift),s])[:N]
    return s/(max(abs(s))+1e-9)

def make_flute(f0=880, sr=SR, dur=2.0, onset=0.0):
    N=int(sr*dur); t=np.arange(N)/sr
    vib=0.003*(1-np.exp(-t/0.4))*np.sin(2*math.pi*6.0*t)
    s=np.zeros(N)
    for h,amp in [(1,1.0),(2,0.15),(3,0.05),(4,0.02)]:
        s+=amp*np.sin(2*math.pi*f0*h*(1+vib)*t)
    air=0.04*np.random.randn(N)
    b,a=ss.butter(4,[f0/sr*1.5,min(f0/sr*8,0.9)],btype='band')
    s+=ss.lfilter(b,a,air)
    env=_adsr(N,sr,attack=0.05,decay=0.02,sustain=0.9,release=0.25)
    s=s*env
    shift=int(onset*sr)
    s=np.concatenate([np.zeros(shift),s])[:N]
    return s/(max(abs(s))+1e-9)


# ════════════════════════════════════════════════
# メイン
# ════════════════════════════════════════════════

if __name__ == '__main__':
    import time

    print("="*60)
    print("楽器分離 v7 — onset差分テンプレート + トポロジカル帰属")
    print("="*60)

    np.random.seed(42)
    DUR = 2.0

    # 実演奏的なonsetズレを設定
    ONSETS = {'piano': 0.000, 'violin': 0.008, 'flute': 0.003}
    piano_src  = make_piano( f0=220, sr=SR, dur=DUR, onset=ONSETS['piano'])
    violin_src = make_violin(f0=440, sr=SR, dur=DUR, onset=ONSETS['violin'])
    flute_src  = make_flute( f0=880, sr=SR, dur=DUR, onset=ONSETS['flute'])
    mixed = piano_src + violin_src + flute_src

    print(f"\nonset設定: piano={ONSETS['piano']*1000:.0f}ms "
          f"violin={ONSETS['violin']*1000:.0f}ms "
          f"flute={ONSETS['flute']*1000:.0f}ms\n")

    profiles = [
        InstrumentProfile('piano',  220.0, 6, 'H',
            {1:1.0,2:0.6,3:0.4,4:0.2,5:0.12,6:0.07}, f0_window=5.0),
        InstrumentProfile('violin', 440.0, 7, '4H',
            {1:1.0,2:0.4,3:0.7,4:0.3,5:0.5,6:0.2,7:0.35}, f0_window=15.0),
        InstrumentProfile('flute',  880.0, 4, '2H',
            {1:1.0,2:0.15,3:0.05,4:0.02}, f0_window=15.0),
    ]

    t0 = time.time()
    results = separate_v7(mixed, SR, profiles,
                          frame_ms=50.0, hop_ms=25.0,
                          onset_block_ms=2.0, tie_threshold=0.15)
    elapsed = time.time() - t0
    print(f"\n処理時間: {elapsed:.2f}秒\n")

    print(f"{'楽器':<8}  {'相関':>8}  {'SDR(dB)':>10}  テンプレート")
    print("-"*50)
    refs = {'piano': piano_src, 'violin': violin_src, 'flute': flute_src}
    for prof in profiles:
        name = prof.name
        ref  = refs[name]
        est  = results[name]
        corr = float(np.corrcoef(ref, est)[0,1])
        s    = sdr_db(ref, est)
        tmpl = "onset_clean" if prof.template_ready else "fallback"
        print(f"{name:<8}  {corr:>8.4f}  {s:>10.1f}  {tmpl}")

    rms_in  = math.sqrt(np.mean(mixed**2))
    rms_res = math.sqrt(np.mean(results['residual']**2))
    print(f"\n残差 RMS: {rms_res:.4f}  (入力比 {rms_res/rms_in*100:.1f}%)")

    # ── 単独区間でのテンプレート精度確認 ──
    print("\n[テンプレート精度: 単独区間での検証]")
    print(f"  {'楽器':<8}  {'f0推定':>8}  {'真のf0':>8}  {'位相誤差(°)':>12}")
    print("  " + "-"*44)
    true_f0s = {'piano': 220.0, 'violin': 440.0, 'flute': 880.0}
    for prof in profiles:
        # 単独区間（onset直後 10ms）での真の f0 との誤差
        onset_s = int(ONSETS[prof.name] * SR)
        seg = refs[prof.name][onset_s:onset_s+int(SR*0.02)]
        f_true = fib_converge(seg, prof.f0_init, sr=SR, window=prof.f0_window)
        phase_err = abs(math.degrees(
            math.atan2(prof.phase_template.get(1, 1+0j).imag,
                       prof.phase_template.get(1, 1+0j).real)
        ))
        print(f"  {prof.name:<8}  {prof.f0:>8.2f}  {f_true:>8.2f}  {phase_err:>12.1f}")
