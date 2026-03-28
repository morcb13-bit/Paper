# 40×40×40格子からフラクタル・フラーレンへ — 波を「構造で巻き取る」

:::message
**2026-03-28 更新:** 実データ（フルート・バイオリン・トランペット・ピアノの C6）での検証結果を追記しました。テストデータの理論がほぼそのまま成立することを確認しています。
:::

---

## はじめに

フーリエ変換は「波を巻き取る」操作です。

3Blue1Brown の動画で見たことがある方も多いと思います。信号を単位円に巻き付けて重心を取ると、巻き付け周波数が信号の周波数と一致したとき重心が大きくずれる。それが周波数の検出です。

この考え方は美しい。ただ、一つ根本的な制約があります。

**「窓」が必要なこと。**

ビブラートのかかったバイオリン、チャープ信号、周波数が時間とともに変化する電磁波。こういう非定常な信号に対して、フーリエ変換は「窓を短くすれば時間分解能が上がるが周波数精度が落ちる」というトレードオフを抱えています。これはハイゼンベルクの不確定性原理と同じ構造です。

この記事では、そのトレードオフを回避する別のアプローチを探る過程を書きます。

出発点は「40×40×40 の整数格子」です。

---

## 1. 40倍則と整数の世界

```python
# step_rule_40: 40倍則
def step_rule_40(T, U):
    return [6*t + 2*u for t, u in zip(T, U)], \
           [2*t - 6*u for t, u in zip(T, U)]

# 逆変換も整数で閉じる
def step_rule_40_inv(T, U):
    return [(6*t + 2*u) // 40 for t, u in zip(T, U)], \
           [(2*t - 6*u) // 40 for t, u in zip(T, U)]
```

変換行列 $M = \begin{pmatrix} 6 & 2 \\ 2 & -6 \end{pmatrix}$ について $M^2 = 40I$ が成り立ちます。

なぜ40か。

```
40 = det(M) の絶対値
40 = lcm(8, 5)         ← 2³ と 5 の最小公倍数
40 = gcd(6²+2², 2²+6²) ← ノルムの共有因子
```

これは偶然ではありません。3次元に展開すると $40 \times 40 \times 40 = 64000$ という格子が現れ、これが後で切断20面体（フラーレン）と結びつきます。

---

## 2. なぜフラーレンか

40×40×40 の3次元整数格子をそのまま使うと計算コストが高い。そこで問いが立ちます。

**「この格子の対称性を保ちながら、より少ない自由度で表現できないか？」**

切断20面体（サッカーボールの形、C₆₀ フラーレン）の頂点数は **60**。

```
BASE = 3120 = 60 × 52
           = 頂点数 × 活性状態数
           = 2⁴ × 3 × 5 × 13
```

この `BASE = 3120` という数が鍵です。

```
3120 / 5  = 624  → 72°  (Z₅: ペンタゴンの対称性)
3120 / 13 = 240  → Z₁₃ の 1 ステップ
3120 / 60 = 52   → 6°   (最小単位)
```

40×40×40 の3次元格子を、60頂点×52状態のフラーレン格子に「引き直す」。これが設計の核心です。

```
従来の3次元格子:
  各点が独立 → 64000 自由度

フラーレン格子:
  60 頂点が切断20面体上に配置
  → 対称性を内蔵した 3120 自由度
  → 計算が「構造を知っている」
```

---

## 3. 整数で波を巻き取る

フラーレン格子上での正弦波テーブル：

```python
import math

BASE = 3120
COS_TABLE = [round(BASE * math.cos(2 * math.pi * i / BASE)) for i in range(BASE)]
SIN_TABLE = [round(BASE * math.sin(2 * math.pi * i / BASE)) for i in range(BASE)]
```

浮動小数点を使うのはこのテーブル生成の一度だけ。以降はすべて整数演算です。

```python
def winding_int(signal, w):
    """周波数 w で信号を巻き取る。浮動小数点ゼロ。"""
    cx, cy = 0, 0
    for k, s in enumerate(signal):
        p = (w * k) % BASE
        cx += s * COS_TABLE[p]
        cy += s * SIN_TABLE[p]
    return cx, cy  # Python 任意精度整数

def mag_sq(cx, cy):
    return cx * cx + cy * cy
```

確認してみます。

```python
N = BASE  # 3120点

# 混合信号: freq=3 と freq=7
signal = [SIN_TABLE[3*k % BASE] + SIN_TABLE[7*k % BASE] for k in range(N)]

for w in range(10):
    cx, cy = winding_int(signal, w)
    m2 = mag_sq(cx, cy)
    if m2 > 10**15:
        print(f"w={w}: HIT  |重心|² = {m2:.2e}")
    elif m2 > 0:
        print(f"w={w}: ノイズ |重心|² = {m2:.2e}")
```

```
w=3: HIT   |重心|² = 2.31e+20
w=7: HIT   |重心|² = 2.31e+20
その他: ノイズ ~10⁸
```

**S/N 比: 10¹²。** 浮動小数点なしで、周波数分離が完璧に動きます。

---

## 4. 窓が必要なのはなぜか — そして回避策

ここで問題に戻ります。

```python
# ビブラート信号: 周波数が 5±1 で揺らぐ
phase = 0.0
vibrato = []
for k in range(BASE):
    f_k = 5 + math.sin(2 * math.pi * k / 200)
    vibrato.append(round(BASE * math.sin(2 * math.pi * phase / BASE)))
    phase += f_k
```

この信号を窓なしで巻き取ると、ピークが「ぼける」。

```
w=4: 中程度
w=5: ピーク（でも純音より低い）
w=6: 中程度
```

窓を短くすれば追随できますが、周波数精度が落ちる。これがトレードオフです。

**回避策の核心:**

「外から巻き取る」のをやめて、**計算空間そのものを信号の位相空間にする**。

フラーレンの 60 頂点は、アトラクタ上で次の状態を保ちます（核61）：

```
60 セルの LSD 値 = 等間隔格子（間隔 = 52 = BASE/60）
位相速度         = 1LSD 周期(139 ステップ)で Sturmian 列(79 or 80) 格子進む
```

つまり60個のセルが、フラーレン上に 0〜59 の格子点として「等間隔に」並んでいる。この格子の動きが、信号の周波数と「共鳴」するかどうかを見れば、窓なしで周波数の性質が分かります。

---

## 5. 共鳴条件の代数的導出

フラーレンの対称性群：

```
Z₅  (ペンタゴン): BASE/5  = 624
Z₁₃ (B13 の核):  BASE/13 = 240
Z₅ × Z₁₃ ≅ Z₆₅  (gcd(5,13) = 1 より)
```

信号の基音 $f_0$ の倍音列 $f_0, 2f_0, 3f_0, \ldots$ は、$(f_0 \bmod 5,\ f_0 \bmod 13)$ の2次元格子上を周回します。

**定理: フラーレン共鳴の分類**

$$\text{軌道サイズ} = \frac{65}{\gcd(f_0 \bmod 65,\ 65)}$$

| $\gcd(f_0, 65)$ | 共鳴クラス | 軌道 | 意味 |
|:-:|:-:|:-:|:--|
| 1  | 完全共鳴 (full) | 65点 | 全対称性モードを励起 |
| 5  | Pentagon 共鳴 | 13点 | Z₅ 節点固定、Z₁₃ が自由 |
| 13 | Z₁₃ 共鳴 | 5点 | Z₁₃ 原点固定、Z₅ が自由 |
| 65 | 直流 (dc) | 1点 | 原点固定 |

```python
import math

def resonance_label(f0):
    g = math.gcd(f0 % 65, 65)
    if g == 1:  return 'full'
    if g == 5:  return 'pentagon'
    if g == 13: return 'z13'
    return 'dc'

for f0 in [1, 5, 13, 65]:
    print(f"f₀={f0:3d}: {resonance_label(f0)}")
```

```
f₀=  1: full
f₀=  5: pentagon
f₀= 13: z13
f₀= 65: dc
```

---

## 6. signal_to_fullerene: フラーレン座標変換

信号を Z₅×Z₁₃ の65点格子に投影します。

```python
def signal_to_fullerene(signal, max_harmonics=13):
    """
    信号をフラーレン座標（Z₅×Z₁₃ スペクトル）に変換。

    Returns:
        spectrum[m5][m13] = (cx, cy)  完全整数
        m5  ∈ {0,1,2,3,4}    (Z₅ 座標)
        m13 ∈ {0,1,...,12}   (Z₁₃ 座標)
    """
    spectrum = [[(0, 0)] * 13 for _ in range(5)]
    for w in range(1, max_harmonics + 1):
        cx, cy = winding_int(signal, w)
        m5, m13 = w % 5, w % 13
        ocx, ocy = spectrum[m5][m13]
        spectrum[m5][m13] = (ocx + cx, ocy + cy)
    return spectrum
```

**共鳴クラスの判定 — z5_excess 方式:**

単純に「Z₅=0 行のパワー比が閾値以下なら full」と判定すると、`max_harmonics=13` のとき構造的な問題が生じます。`max_harmonics=13` では 65 点格子のうち 13 点しか使われず、そのうち Z₅=0 行には w=5 と w=10 の 2 点が必ず入ります。均等分布なら 2/13 ≈ 15.4% になるので、**純粋な full 信号でも旧閾値 0.1 を超えてしまいます。**

解決策は「背景レベル（均等分布時の期待値）を差し引いた超過量」で判定することです。

```python
def resonance_class(spectrum, max_harmonics=13):
    """
    z5_excess 方式でフラーレン共鳴クラスを判定。

    z5_background = (max_harmonics // 5)  / max_harmonics  # 均等分布の期待値
    z5_excess     = z5_node_ratio - z5_background          # 背景差し引き後

    z5_excess > +0.70 → pentagon
    z5_excess < -0.05 かつ z13_excess < -0.05 → full
    それ以外 → mixed

    ⚠️ max_harmonics は signal_to_fullerene() に渡した値と一致させること。
    """
    z5_0  = sum(mag_sq(*spectrum[0][m13]) for m13 in range(13))
    z13_0 = sum(mag_sq(*spectrum[m5][0])  for m5  in range(5))
    total = sum(mag_sq(*spectrum[m5][m13])
                for m5 in range(5) for m13 in range(13))

    if total == 0:
        return {'label': 'zero', 'z5_node_ratio': 0.0, 'z13_node_ratio': 0.0,
                'z5_excess': 0.0, 'z13_excess': 0.0}

    z5_r  = z5_0  / total
    z13_r = z13_0 / total

    mh     = max_harmonics
    z5_bg  = (mh // 5)  / mh   # ≈ 0.20（mh が大きいとき）、2/13 ≈ 0.154（mh=13）
    z13_bg = (mh // 13) / mh

    z5_excess  = z5_r  - z5_bg
    z13_excess = z13_r - z13_bg

    if z5_excess  > 0.70:                              label = 'pentagon'
    elif z13_excess > 0.70:                            label = 'z13'
    elif z5_excess < -0.05 and z13_excess < -0.05:    label = 'full'
    else:                                              label = 'mixed'

    return {
        'label':          label,
        'z5_node_ratio':  z5_r,   'z13_node_ratio': z13_r,
        'z5_excess':      z5_excess, 'z13_excess':  z13_excess,
        'z5_bg':          z5_bg,   'z13_bg':        z13_bg,
        'total':          total,
    }
```

**検証:**

```python
N = BASE

sig_full  = [SIN_TABLE[k    % BASE] for k in range(N)]  # f₀=1
sig_penta = [SIN_TABLE[5*k  % BASE] for k in range(N)]  # f₀=5
sig_z13   = [SIN_TABLE[13*k % BASE] for k in range(N)]  # f₀=13

for label, sig in [('f₀=1', sig_full), ('f₀=5', sig_penta), ('f₀=13', sig_z13)]:
    sp = signal_to_fullerene(sig, max_harmonics=13)
    rc = resonance_class(sp, max_harmonics=13)
    print(f"{label}: z5_excess={rc['z5_excess']:+.4f}  → {rc['label']}")
```

```
f₀=1:  z5_excess=-0.1538  → full
f₀=5:  z5_excess=+0.8462  → pentagon
f₀=13: z5_excess=-0.1538  → z13
```

3クラスが完全に分離されました。

---

## 7. ビブラートの代数的検出

```python
# ビブラート深さ = Z₅節点比の低下量
phase = 0.0
sig_penta = [SIN_TABLE[5*k % BASE] for k in range(N)]
vibrato   = []
for k in range(N):
    f_k = 5 + math.sin(2 * math.pi * k / 200)
    vibrato.append(round(BASE * math.sin(2 * math.pi * phase / BASE)))
    phase += f_k

sp_pure = signal_to_fullerene(sig_penta, max_harmonics=13)
sp_vib  = signal_to_fullerene(vibrato,   max_harmonics=13)

rc_pure = resonance_class(sp_pure, max_harmonics=13)
rc_vib  = resonance_class(sp_vib,  max_harmonics=13)

depth = rc_pure['z5_node_ratio'] - rc_vib['z5_node_ratio']
print(f"ビブラート深さ指標: {depth:.6f}")
```

```
ビブラート深さ指標: 0.000692
```

**窓関数なし、時間分解能のトレードオフなし。** ビブラートの深さが代数的な一つの数として現れます。

---

## 8. フラーレン格子の全体像

```
        Z₁₃座標 (0〜12)
        ────────────────────────────►
   Z₅  [0]  ・  ・  ・  ・  ・ ← Pentagon 共鳴（行）
   座  [1]  ・  ・  ・  ・  ・
   標  [2]  ・  ・  ・  ・  ・
   (  [3]  ・  ・  ・  ・  ・
   0  [4]  ・  ・  ・  ・  ・
   〜       ↑
   4)       Z₁₃ 共鳴（列）

  f₀=5  の倍音: Z₅=0 行に集中 → Pentagon 共鳴
  f₀=13 の倍音: Z₁₃=0 列に集中 → Z₁₃ 共鳴
  f₀=1  の倍音: 65 点全体に分散 → 完全共鳴
```

この格子は「周波数空間」ではなく**「対称性空間」**です。信号を外から巻き取るのではなく、フラーレンの構造が信号の代数的性質を直接読み出します。

---

## 9. 4H回避という自然フィルタ

B13 フラクタルエンジンには「4H 回避」という性質があります（核57）。

Fibonacci 数列を mod 13 で見ると、{4,6,7,9}（4H 帯）に一度も入らない。

```python
def fib_mod13_sequence(n):
    a, b = 0, 1
    for _ in range(n):
        yield a % 13
        a, b = b, a + b

four_h = {4, 6, 7, 9}
seq    = list(fib_mod13_sequence(28))
print([v for v in seq if v in four_h])  # []
```

```
[]  ← 4H 帯への到達はゼロ
```

これは「フラーレンアトラクタが 4H 帯の倍音に自然に応答しない」ことを意味します。楽器音・声・電磁波の倍音列のうち、4H 帯に落ちる成分は構造的にフィルタされる。

この性質は設計ではなく、Fibonacci 数列の代数的性質から来る **必然** です。

---

## 10. 実データへの適用 — N≠BASE の問題と解決

ここまではすべて N=BASE=3120 のテストデータの話でした。実データは異なるサンプリングレートと長さを持ちます。

### 発生する問題: 非整数周期とスペクトル漏れ

サンプリングレート sr=11025Hz の WAV で C6（1046.5Hz）を扱う場合、1周期のサンプル数は：

```
11025 / 1046.5 = 10.535... サンプル（非整数）
```

これは Sturmian 数列の構造を持ちます（連分数展開: [10; 1, 1, 6, 1, 1]）。非整数周期をそのまま `winding_int` にかけると、**スペクトル漏れ（spectral leakage）** が生じます。具体的には、C6 に対応する w=296 の隣の w=295（Z₅=0 の pentagon 座席）にエネルギーが流れ込み、本来 full のはずの信号が pentagon と誤判定されます。

### 解決策: 2層設計

```
Layer 1 (常用):  max_harmonics = 13
  → 低域（〜46Hz）のエネルギー分布を見る
  → C6 の基音は範囲外だが、楽器種別の低域構造差として現れる
  → 高速・固定コスト

Layer 2 (f0 既知時): max_harmonics = w_f0 * 2
  → w_f0 = round(BASE * f0 / sr)
  → C6: w_f0 = round(3120 × 1046.5 / 11025) = 296、mh = 592
  → 基音の両隣をカバーすることで漏れが相殺され、正確な判定が得られる
  → *1 では漏れが残る。*2 が必須
```

```python
import wave, struct

def read_wav_scaled(path, base=3120):
    """WAV ファイルを BASE スケールの整数列に変換。16bit モノラル専用。"""
    with wave.open(path) as w:
        sr  = w.getframerate()
        raw = w.readframes(w.getnframes())
        sw  = w.getsampwidth()
    samples = struct.unpack(f'<{len(raw)//sw}h', raw)
    return [round(v * base / 32768) for v in samples], sr

def wav_to_fullerene(path, f0_hz=None, layer=1, base=3120):
    """
    WAV ファイルをフラーレン時系列フレームに変換。
    layer=1: max_harmonics=13（低域、高速）
    layer=2: max_harmonics=w_f0*2（基音帯、f0_hz 必須）
    """
    samples, sr = read_wav_scaled(path, base)
    if layer == 1:
        mh = 13
    else:
        w_f0 = round(base * f0_hz / sr)
        mh   = w_f0 * 2
    n_blocks = len(samples) // base
    frames   = []
    for i in range(n_blocks):
        seg = samples[i * base: (i+1) * base]
        sp  = signal_to_fullerene(seg, max_harmonics=mh)
        rc  = resonance_class(sp, max_harmonics=mh)
        frames.append({'frame': i, 'spectrum': sp, 'resonance': rc})
    return frames, mh, sr
```

:::message
**frame_size = BASE = 3120 に固定する理由:**
`winding_int` の直交性は 1 周期以上の積分が成立して初めて現れます。`frame_size < 1周期` だと判定が崩壊し、全フレームが mixed になります。実データで N=3120 ではない場合でも、BASE 単位でブロックに切って処理します。
:::

---

## 11. 実データ検証: 4楽器 C6

フルート・バイオリン・トランペット・ピアノの C6 音（sr=11025Hz）で検証しました（各3〜7秒のモノラル WAV）。

### Layer 2（mh=592）での成分比率

```
楽器       pentagon%  z13%   full%
flute        5.6%     2.5%   91.8%
violin       9.2%     1.5%   90.7%
trumpet      8.7%     0.3%   91.0%
piano        3.0%     1.0%   96.0%
```

4楽器ともに full が主体（91〜96%）でした。pentagon 比率は violin が最大で、piano が最小です。これは倍音構造の楽器差がフラーレン格子上に現れたものです。

### 成分分離の完全性

pentagon 成分を切り出して引き算すると、残渣の Z₅=0 行パワーは **完全にゼロ** になります。

```python
def spectrum_mask_z5(spectrum):
    """Z₅=0 行のみ残す（pentagon 成分を取り出す）。"""
    return [
        [(spectrum[m5][m13] if m5 == 0 else (0, 0)) for m13 in range(13)]
        for m5 in range(5)
    ]

def spectrum_subtract(sp_a, sp_b):
    """スペクトルの差分。完全整数演算。"""
    return [
        [(sp_a[m5][m13][0] - sp_b[m5][m13][0],
          sp_a[m5][m13][1] - sp_b[m5][m13][1])
         for m13 in range(13)]
        for m5 in range(5)
    ]

# pentagon 成分を除去したあとの Z₅=0 行残留
sp_pent     = spectrum_mask_z5(sp)
sp_residual = spectrum_subtract(sp, sp_pent)

z5_residual = sum(mag_sq(*sp_residual[0][m13]) for m13 in range(13))
print(z5_residual)  # → 0  （全楽器で確認済み）
```

### 4楽器合成・逐次分離テスト

4楽器の信号を加算して合成し、既知のスペクトルを1つずつ引いていったとき、最後の残渣がゼロになるかを確認しました。

```
操作                     残りパワー    全体比
合成(4楽器)              3.530e+18    100.0%
  - trumpet              3.529e+18     99.98%
  - violin               2.064e+18     58.47%
  - flute                1.692e+15      0.05%
  - piano                0.000e+00      0.00%  ← 完全消滅 ✓
非ゼロセル数                0 / 65             ← 全格子点ゼロ ✓
```

整数演算のため、位相も含めて誤差なく分離できています。

これはスペクトル空間での「完全分離」が実データでも成立することの証明です。

---

## まとめ

**確立したこと（テストデータ＋実データで検証済み）：**

| 項目 | 結果 |
|:--|:--|
| 整数演算巻き取り | S/N 比 10¹²、浮動小数点ゼロ |
| 共鳴クラス分類 | Z₅×Z₁₃≅Z₆₅ の代数的分類 |
| 3クラス完全分離 | z5_excess 方式で mh 依存なく判定 |
| ビブラート検出 | 窓不要、代数的深さ指標 |
| N≠BASE 対応 | Layer1/Layer2 の2層設計 |
| 実データ分離 | 4楽器合成を逐次引き算 → 残渣ゼロ ✓ |

**これから：**

- 楽器の pentagon/z13/full 比率の **時系列変化**（アタック〜サスティン〜リリース）
- 同一楽器で音高を変えたときの比率変化（音高依存か楽器固有か）
- MRC / クライオ EM 画像への適用

---

## コード（自己完結版）

以下は標準ライブラリのみで動作する最小実装です。

```python
import math

BASE = 3120
COS_TABLE = [round(BASE * math.cos(2*math.pi*i/BASE)) for i in range(BASE)]
SIN_TABLE = [round(BASE * math.sin(2*math.pi*i/BASE)) for i in range(BASE)]

def winding_int(signal, w):
    cx, cy = 0, 0
    for k, s in enumerate(signal):
        p = (w * k) % BASE
        cx += s * COS_TABLE[p]
        cy += s * SIN_TABLE[p]
    return cx, cy

def mag_sq(cx, cy): return cx*cx + cy*cy

def signal_to_fullerene(signal, max_harmonics=13):
    spectrum = [[(0, 0)] * 13 for _ in range(5)]
    for w in range(1, max_harmonics + 1):
        cx, cy = winding_int(signal, w)
        m5, m13 = w % 5, w % 13
        ocx, ocy = spectrum[m5][m13]
        spectrum[m5][m13] = (ocx + cx, ocy + cy)
    return spectrum

def resonance_class(spectrum, max_harmonics=13):
    z5_0  = sum(mag_sq(*spectrum[0][j]) for j in range(13))
    z13_0 = sum(mag_sq(*spectrum[i][0]) for i in range(5))
    total = sum(mag_sq(*spectrum[i][j]) for i in range(5) for j in range(13))
    if not total:
        return {'label': 'zero', 'z5_node_ratio': 0.0, 'z13_node_ratio': 0.0,
                'z5_excess': 0.0, 'z13_excess': 0.0}
    z5_r, z13_r = z5_0 / total, z13_0 / total
    mh = max_harmonics
    z5_exc  = z5_r  - (mh // 5)  / mh
    z13_exc = z13_r - (mh // 13) / mh
    if   z5_exc  >  0.70:                    label = 'pentagon'
    elif z13_exc >  0.70:                    label = 'z13'
    elif z5_exc  < -0.05 and z13_exc < -0.05: label = 'full'
    else:                                    label = 'mixed'
    return {'label': label, 'z5_node_ratio': z5_r, 'z13_node_ratio': z13_r,
            'z5_excess': z5_exc, 'z13_excess': z13_exc}

def spectrum_subtract(sp_a, sp_b):
    return [[(sp_a[i][j][0]-sp_b[i][j][0], sp_a[i][j][1]-sp_b[i][j][1])
             for j in range(13)] for i in range(5)]

def spectrum_mask_z5(sp):
    return [[(sp[i][j] if i == 0 else (0, 0)) for j in range(13)] for i in range(5)]
```

---

*続きは次のノートで。*
