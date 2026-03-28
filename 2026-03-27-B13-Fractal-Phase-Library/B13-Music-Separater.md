# coneFFT で楽器を分離する — 倍音列の整数構造を使う

:::message
この記事は「[40×40×40格子からフラクタル・フラーレンへ](https://zenn.dev/)」の続きです。フラーレン変換（`winding_int`、`signal_to_fullerene`、`resonance_class`）の基礎は前の記事を参照してください。
:::

---

## アイデア：倍音列は整数比で「揃っている」

楽器の音はほぼ必ず倍音列を持ちます。

```
ピアノの C4 (261.6Hz) の倍音列:
  261.6Hz, 523.2Hz, 784.8Hz, 1046.4Hz, ...
  → f0, 2f0, 3f0, 4f0, ...
```

winding 変換の文脈で言えば、基音が winding 番号 `w_f0` に対応するなら、倍音は `2*w_f0`、`3*w_f0`、`4*w_f0` … に対応します。

**整数比で並んでいる。**

これは非常に強い構造的制約です。同じ `w` を共有しない音源は、winding 空間（= FFT空間）で自然に分離できます。

```
ピアノ C4: w=74, 148, 222, 296, 370, ...  (74の倍数)
ボーカル D4: w=83, 166, 249, 332, 415, ... (83の倍数)
ベース  G2: w=47,  94, 141, 188, 235, ...  (47の倍数)
```

これらは重なりません。

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

- **フルート**: k=1 が 97.6%、ほぼ純音
- **ピアノ**: k=2 が 77.2%（k=1 より強い）
- **トランペット**: k=2 と k=4 が突出
- **バイオリン・ギター**: k=1〜4 が均等（自然減衰型）

コサイン類似度で見ると：

```
piano  ↔ trumpet: 0.907   ← 「偶数倍音型」クラスタ
violin ↔ guitar:  0.945   ← 「自然減衰型」クラスタ
flute  ↔ guitar:  0.937   ← 同じく自然減衰型
piano  ↔ flute:   0.112   ← 正反対の構造
```

**バイオリンとギターが高類似**なのは、どちらも筐体（胴体）共鳴が k=1 を際立たせるからです。弦の振動 + 共鳴胴という構図が共通しています。

この指紋ベクトルを使うと、「検出した f0 候補が何の楽器か」をパターンマッチングで識別できる可能性があります（未実装・次の課題）。

---

## 打楽器は rest に自然に残る

実際の楽曲（和太鼓を含む管弦楽）で試したとき、打撃音が `rest` トラックによく残ることを確認しました。

これは偶然ではありません。

```
和太鼓の打撃音:
  広帯域ノイズ + 胴鳴り
  → 倍音列の整数比構造を持たない
  → どのグループにも引っかからない
  → rest に自然に残る
```

逆に言えば、`rest` トラックは「調性を持たない音」のコレクションです。打楽器・打撃音・環境ノイズが集まります。

```
倍音列グループ = 調性楽器トラック
rest          = 非調性音トラック（打楽器・ノイズ）
```

これは設計したものではなく、**倍音構造の有無という物理的性質から自然に出てきた分離**です。

---

## Claude に読み込ませて試す

この記事で紹介したコードを Claude に読み込ませ、手持ちの楽曲をアップロードすると、楽器ごとの分離 WAV を生成できます。

```
手順:
1. 上記のコードを Python ファイルとして保存
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

## 現状の限界と次の課題

**現状の限界：**

```
- 固定 f0 を取るため、コードチェンジで同楽器が別 f0 に移動すると分離が崩れる
- 同じ f0 を複数楽器が共有するケース（ユニゾン）では分離不可
- greedy 方式のため先に取ったグループに後のグループの成分が混入する
```

**次の課題：**

```
1. フレームごとの f0 動的追跡（コードチェンジ対応）
2. 倍音パターンマッチング（FINGERPRINTS との照合）
3. coneFFT の resonance_label と楽器種別の対応検証
   → Layer1 の z5/z13 比率差が楽器の「音色」と対応するか？
```

---

## まとめ

| 項目 | 内容 |
|:--|:--|
| 分離原理 | 倍音列の整数比構造（f0, 2f0, 3f0 ...）を FFT ビン空間で識別 |
| f0 検出 | 安定スコア（出現頻度 × 平均パワー）で曲全体の主要 f0 を抽出 |
| 分離方式 | 各 f0 の倍音列に対応するビンを取り出して逆 FFT |
| 打楽器 | 整数比構造なし → rest に自然分離 |
| 楽器指紋 | 倍音パワー比ベクトルとコサイン類似度でクラスタ化 |
| 処理速度 | 2500 ブロック（約 12 分）が 1 秒以下 |

倍音列という「整数の構造」を使うことで、外部の学習データも複雑なモデルも使わず、代数的に楽器を分離できます。完全な分離はまだ先ですが、「万能巻き取り機」の方向性として面白い結果が出始めています。

---

*続きは次のノートで。*
