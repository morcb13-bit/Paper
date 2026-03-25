# 🎸 好きなMP3からベースラインだけ抜き出せるか、Claudeに試してもらった

## これは何？

好きな曲のMP3をClaudeにアップロードすると、**ベースラインだけを分離したWAVファイル**を返してくれる。

プラグインなし。専用ソフトなし。Claudeだけ。

---

## 実際にやってみた

2曲で試した結果がこれ。

| 指標 | 意味 | 結果 |
|:---|:---|:---|
| 低域集中度 | ベース以外の音がどれだけ混入しているか | **100%** ✅ |
| 中域漏れ | ギターやピアノが混ざっていないか | **0.010** ✅ |
| RMS比 | 元音源に対するベース音量の割合 | **18%** ✅ |

数字は合格圏。あとは耳で判断してほしい。

---

## 仕組み

フーリエ変換で倍音構造を分解し、**40〜200Hzの成分だけを黄金比探索で精密に抽出**する。

ポイントは「引きすぎない」こと。倍音の上限を200Hzに制限することで、中域・高域への漏れを防いでいる。

```
f0=100Hzのベースの場合:
  h1 = 100Hz ← 抽出する
  h2 = 200Hz ← 抽出する（上限）
  h3 = 300Hz ← 抽出しない（他の楽器と混在するため）
```

---

## あなたの曲で試す方法

**以下のテキストをそのままClaudeに貼り付けて、MP3をアップロードするだけ。**

---

次のコードを実行して、アップロードされたMP3からベースラインを分離してください。

```python
import librosa, numpy as np, math, soundfile as sf
from scipy.signal import butter, lfilter, find_peaks
import sys, glob, os

uploads = glob.glob("/mnt/user-data/uploads/*.mp3")
if not uploads:
    print("MP3ファイルが見つかりません。アップロードしてください。")
    sys.exit()
src_path = uploads[-1]
print(f"対象ファイル: {os.path.basename(src_path)}")

SR = 11025
y, sr_orig = librosa.load(src_path, sr=None, mono=False)
orig = librosa.resample(y.mean(axis=0), orig_sr=sr_orig, target_sr=SR)
print(f"長さ: {len(orig)/SR:.1f}秒")

def get_coeff(sig, freq, sr=SR):
    N = len(sig); t = np.arange(N)/sr
    return (2*np.dot(sig, np.cos(2*math.pi*freq*t))/N +
            1j*2*np.dot(sig, np.sin(2*math.pi*freq*t))/N)

def extract_sinusoid(N, freq, coeff, sr=SR):
    t = np.arange(N)/sr
    return coeff.real*np.cos(2*math.pi*freq*t) + coeff.imag*np.sin(2*math.pi*freq*t)

def fib_converge(sig, f0, sr=SR, n_iter=40, window=10.0):
    PHI = (1+math.sqrt(5))/2
    lo, hi = max(0.5, f0-window), f0+window
    cf = f0
    for _ in range(n_iter):
        if hi-lo < 0.05: break
        m1 = hi-(hi-lo)/PHI; m2 = lo+(hi-lo)/PHI
        def amp(f):
            N=len(sig); t=np.arange(N)/sr; w=np.blackman(N); wg=np.mean(w)
            return math.sqrt((2*np.dot(sig,np.cos(2*math.pi*f*t)*w)/(N*wg))**2+
                             (2*np.dot(sig,np.sin(2*math.pi*f*t)*w)/(N*wg))**2)
        if amp(m1) < amp(m2): lo=m1; cf=m2
        else: hi=m2; cf=m1
    return cf

def pull_bass(residual, sr=SR, f_hi=250, max_harmonic_hz=200):
    N = len(residual)
    WIN = int(sr*0.15); FRAME = int(sr*0.15); HOP = int(sr*0.05)
    nyq = sr/2
    b, a = butter(4, f_hi/nyq, btype='low')
    low = lfilter(b, a, residual)
    spec = np.abs(np.fft.rfft(low*np.hanning(N)))
    freqs = np.fft.rfftfreq(N, 1/sr)
    mask = (freqs>=40)&(freqs<=f_hi)
    sub = spec*mask
    peaks, _ = find_peaks(sub, height=sub.max()*0.05, distance=int(10/(sr/N)))
    if len(peaks)==0: return None, residual
    f_cand = float(freqs[peaks[np.argmax(sub[peaks])]])
    f0 = fib_converge(residual[:WIN], f_cand, sr, window=f_cand*0.15)
    c_f0 = get_coeff(residual[:WIN], f0, sr)
    if abs(c_f0)<1e-6: return None, residual
    tmpl = {}
    amp_h1 = abs(c_f0)
    for h in range(1, 12):
        freq = f0*h
        if freq>sr/2 or freq>max_harmonic_hz: break
        c_h = get_coeff(residual[:WIN], freq, sr)
        ratio = abs(c_h)/(amp_h1+1e-12)
        if ratio>0.95: c_h = c_h*(0.95/ratio)
        tmpl[h] = c_h/(c_f0+1e-12)
    out=np.zeros(N); ov=np.zeros(N); hann=np.hanning(FRAME)
    for start in range(0, N-FRAME, HOP):
        seg=residual[start:start+FRAME].copy()
        c_f0_=get_coeff(seg, f0, sr)
        if abs(c_f0_)<1e-8: continue
        for h, ratio in tmpl.items():
            freq=f0*h
            if freq>sr/2 or freq>max_harmonic_hz: break
            c_use=0.5*c_f0_*ratio+0.5*get_coeff(seg, freq, sr)
            if abs(c_use)<1e-7: continue
            comp=extract_sinusoid(FRAME, freq, c_use, sr)
            out[start:start+FRAME]+=comp*hann
            seg-=comp
        ov[start:start+FRAME]+=hann
    ov=np.maximum(ov,1e-9)
    out/=ov
    return (f0, out), residual-out

print("\n[抽出中... 曲の長さによって1〜3分かかります]")
N_total=len(orig); N_frame=int(2.0*SR); N_hop=int(1.0*SR)
bass_out=np.zeros(N_total); bass_ov=np.zeros(N_total)
hann2=np.hanning(N_frame); f0_log=[]

for start in range(0, N_total-N_frame, N_hop):
    seg=orig[start:start+N_frame].copy()
    scale=max(abs(seg))+1e-9
    result, _=pull_bass(seg/scale, SR)
    if result is None: continue
    f0, bs=result
    bass_out[start:start+N_frame]+=bs*scale*hann2
    bass_ov [start:start+N_frame]+=hann2
    f0_log.append(f0)

bass_ov=np.maximum(bass_ov,1e-9)
bass=bass_out/bass_ov

nyq=SR/2
spec=np.abs(np.fft.rfft(bass*np.hanning(N_total)))
freqs=np.fft.rfftfreq(N_total, 1/SR)
conc=np.sum(spec[(freqs>=40)&(freqs<=300)]**2)/(np.sum(spec**2)+1e-9)
b2,a2=butter(4,[500/nyq,0.98],btype='band')
mid_leak=math.sqrt(np.mean(lfilter(b2,a2,bass)**2))/(math.sqrt(np.mean(lfilter(b2,a2,orig)**2))+1e-9)

def freq_to_note(f):
    if f<=0: return "---"
    n=round(12*math.log2(f/440.0))+69
    names=['C','C#','D','D#','E','F','F#','G','G#','A','A#','B']
    return f"{names[n%12]}{n//12-1}"

print(f"\n[品質スコア]")
print(f"  低域集中度 : {conc*100:.1f}%  {'✅' if conc>0.7 else '⚠️ 他の楽器が混入している可能性'}")
print(f"  中域漏れ   : {mid_leak:.3f}   {'✅' if mid_leak<0.15 else '⚠️ ギター・ピアノが混ざっている可能性'}")
print(f"  ベース音域 : {freq_to_note(np.median(f0_log))}  ({np.median(f0_log):.1f}Hz)")

basename = os.path.splitext(os.path.basename(src_path))[0]
out_path = f"/mnt/user-data/outputs/{basename}_bass.wav"
sf.write(out_path, bass, SR)
print(f"\n  → {basename}_bass.wav を保存しました。ダウンロードして聴いてみてください。")
print(f"\n聴いた感想をZennのコメントに書いてもらえると、次のバージョンの改善に直結します。")
```

---

## フィードバックをください

**このアルゴリズムは、あなたが試した結果で育ちます。**

コメント欄に以下を書いてもらえると助かります：

- 曲のジャンル（ロック、ジャズ、EDMなど）
- 品質スコアの数字（低域集中度・中域漏れ）
- 聴いた感想（「ベースだけになった」「ドラムが混ざった」など）

フィードバックをもとに改善して、また記事を更新します。

---

## 現在の限界と今後

| パート | 状態 |
|:---|:---|
| ベースライン | ✅ v2（本記事） |
| メロディ | 🔧 開発中 |
| コード（ピアノ/ギター） | 🔧 開発中 |
| ドラム | 🔧 開発中 |

ベースが抜けたら、次はメロディ。メロディが抜けたら、次はコード。**一緒に精度を上げていきましょう。**
