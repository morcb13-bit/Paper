## B13で音楽信号を分離する ー Claude music separator　ー
## Claude に読み込ませて試せます

B13phaseを使った音楽信号分離器です。
この記事の Python コード（面倒なら記事全体でもよい）を Claude に読み込ませ、手持ちの楽曲をアップロードすると、楽器ごとに分離したWAVデータが生成できます。

### 楽器とボーカルの分離ではなく、音楽から、楽器を引き算していくイメージです。
### 最後に残ったrest.midには、主に、ボーカルや打楽器などの信号が残ります。

※ 現状は、最小限の機能しか入れてないので、楽器の分離も不充分です。
　 あなた自身が、Claudeに適切な指示を出して学習させることで、もっと色々楽しめるようになるかもしれません。

```
手順:
1. この記事の下記のコードをコピーする（記事を丸ごとコピーでも大丈夫です）
2. Claude の会話に貼り付けるか、ファイルとして添付
3. WAV または MP3/WMA 形式の楽曲をアップロード
4. 「この曲で楽器分離をやってください」と指示
```

Claude は自動的に以下を実行します。

```
1. ffmpeg で WAV に変換（sr=11025Hz, モノラル）
2. f0 安定スコアをスキャン（数秒）
3. 上位 f0 をグループとして分離 WAV を出力
4. rest（残り成分）も出力
```

:::message alert
著作権のある楽曲は個人で楽しむ範囲でのみ使用してください。分離した音声の公開・配布は著作権法上の問題になる可能性があります。
:::

---

## 安定 f0 スコアリング

楽曲全体を通じて「どの f0 が何回出現し、どれだけのパワーを持っていたか」を集計します。

出現頻度が高く、平均パワーが大きい f0 ＝ 曲全体に渡って鳴り続ける楽器、です。

```python
import numpy as np
import wave
from collections import defaultdict

BASE = 3120
sr   = 11025

def scan_stable_f0s(wav_path, n_groups=8, f0_min=20, f0_max=200):
    with wave.open(wav_path) as wf:
        raw = wf.readframes(wf.getnframes())
    sig = np.frombuffer(raw, dtype=np.int16).astype(float)
    n_blocks = len(sig) // BASE

    f0_count       = defaultdict(int)
    f0_total_power = defaultdict(float)

    for b in range(n_blocks):
        seg     = sig[b*BASE:(b+1)*BASE]
        F       = np.fft.rfft(seg, n=BASE)
        power   = np.abs(F)**2
        total_p = power.sum()
        if total_p == 0:
            continue
        remaining = power.copy()

        for _ in range(n_groups):
            best_f0, best_score = None, 0.0
            for f0 in range(f0_min, f0_max + 1):
                score, hit = 0.0, 0
                for k in range(1, 9):
                    w = f0 * k
                    if w >= len(remaining):
                        break
                    if remaining[w] > 0:
                        hit += 1
                    score += remaining[w]
                if hit >= 3 and score > best_score:
                    best_score, best_f0 = score, f0
            if best_f0 is None:
                break
            f0_count[best_f0]       += 1
            f0_total_power[best_f0] += best_score / total_p
            for k in range(1, 9):
                w = best_f0 * k
                if w < len(remaining):
                    remaining[w] = 0

    # 安定スコア = 出現回数 × 平均パワー
    stable = sorted(
        f0_count.keys(),
        key=lambda f: -(f0_count[f] * f0_total_power[f] / max(f0_count[f], 1))
    )
    return stable, f0_count, f0_total_power
```

手持ちのインストゥルメンタル曲で試すと、こういう結果になります（曲名は非公開）。

```
f0安定ランキング（11分の管弦楽曲）:
  w32 (113.1Hz)  出現= 963  avg=14.24%  score=13718  ← 低域弦コード
  w21 ( 74.2Hz)  出現= 964  avg=10.62%  score=10242  ← ベース
  w47 (166.1Hz)  出現= 936  avg= 8.62%  score= 8069  ← 中域コード
  w95 (335.7Hz)  出現= 553  avg= 7.32%  score= 4047  ← メロディ
  w63 (222.6Hz)  出現= 374  avg= 9.67%  score= 3617  ← 中高域
```

---

## 分離：FFT ビンの割り当て

安定 f0 が分かったら、各 f0 の倍音列に対応する FFT ビンをグループごとに取り出すだけです。

```python
def separate_instruments(wav_path, groups, output_prefix):
    """
    groups: {'楽器名': f0_winding番号} の辞書
    例: {'bass': 21, 'low_chord': 32, 'melody': 95}
    """
    with wave.open(wav_path) as wf:
        raw = wf.readframes(wf.getnframes())
        sr  = wf.getframerate()
    sig      = np.frombuffer(raw, dtype=np.int16).astype(float)
    n_blocks = len(sig) // BASE

    # 各グループの倍音列 winding 番号セット
    def wset(f0, w_max=500):
        return {f0*k for k in range(1, 10) if f0*k <= w_max}

    group_wsets = {name: wset(f0) for name, f0 in groups.items()}

    out      = {name: np.zeros(n_blocks * BASE) for name in groups}
    out['rest'] = np.zeros(n_blocks * BASE)

    for b in range(n_blocks):
        seg    = sig[b*BASE:(b+1)*BASE]
        F      = np.fft.rfft(seg, n=BASE)
        F_rest = F.copy()

        for name, ws in group_wsets.items():
            F_g = np.zeros_like(F)
            for w in ws:
                if w < len(F):
                    F_g[w]    = F[w]
                    F_rest[w] = 0
            out[name][b*BASE:(b+1)*BASE] = np.fft.irfft(F_g, n=BASE)

        out['rest'][b*BASE:(b+1)*BASE] = np.fft.irfft(F_rest, n=BASE)

    def save(arr, path):
        mx = np.abs(arr).max()
        if mx > 0:
            arr = arr / mx * 28000
        arr = np.clip(arr, -32767, 32767).astype(np.int16)
        with wave.open(path, 'w') as wf:
            wf.setnchannels(1)
            wf.setsampwidth(2)
            wf.setframerate(sr)
            wf.writeframes(arr.tobytes())

    for name, arr in out.items():
        save(arr, f'{output_prefix}_{name}.wav')
        print(f'saved: {output_prefix}_{name}.wav')
```

---

## 倍音パターン指紋

面白いことに、楽器種別によって「どの倍音 k にパワーが集中するか」が異なります。

合成音（C4）で確認した倍音パワー比（k=1〜8、正規化）：

```python
FINGERPRINTS = {
    'flute':   [0.976, 0.021, 0.002, 0.000, 0.000, 0.000, 0.000, 0.000],
    'guitar':  [0.613, 0.216, 0.094, 0.051, 0.021, 0.005, 0.000, 0.000],
    'violin':  [0.406, 0.256, 0.191, 0.093, 0.032, 0.013, 0.007, 0.003],
    'piano':   [0.071, 0.772, 0.120, 0.029, 0.007, 0.002, 0.000, 0.000],
    'trumpet': [0.019, 0.473, 0.166, 0.217, 0.067, 0.035, 0.014, 0.007],
}
```

