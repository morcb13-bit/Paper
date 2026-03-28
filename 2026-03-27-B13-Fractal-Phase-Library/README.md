# b13phase — B13 Fractal Phase Library

**BASE = 3120** を基数とするフラーレン座標変換ライブラリ。  
信号を Z₅×Z₁₃ の65点格子（切断20面体の対称性）に投影し、共鳴クラスを整数演算のみで判定します。

```
S/N比: 10¹²   浮動小数点演算: ゼロ（テーブル生成時の1回のみ）
```

---

## 概要

フーリエ変換の「窓のトレードオフ（時間分解能 vs 周波数精度）」を回避するアプローチです。

信号を外から巻き取るのではなく、**フラーレン（C₆₀）の対称性群 Z₅×Z₁₃≅Z₆₅ を計算空間として使う**ことで、信号の代数的性質（共鳴クラス）を窓なしで判定します。

```
BASE = 3120 = 2⁴ × 3 × 5 × 13 = 60頂点 × 52活性状態

共鳴クラス分類:
  gcd(f₀, 65) = 1  → full     （全対称性モードを励起）
  gcd(f₀, 65) = 5  → pentagon （Z₅節点固定）
  gcd(f₀, 65) = 13 → z13      （Z₁₃原点固定）
  gcd(f₀, 65) = 65 → dc       （直流）
```

詳細は [Zenn 記事](https://zenn.dev/) を参照してください。

---

## 動作確認環境

- Python 3.10 以上
- 標準ライブラリのみ（`wave`, `struct`, `math`）
- WAV 解析時は `numpy` を推奨（f0 推定に使用）

---

## インストール

```bash
git clone https://github.com/yourname/b13phase.git
cd b13phase
python -m pytest  # テスト（省略可）
```

pip パッケージは未公開です。`b13phase_v075/` ディレクトリをそのまま使ってください。

```python
import sys
sys.path.insert(0, '/path/to/b13phase')
import b13phase_v075 as b13
```

---

## クイックスタート

### テストデータで共鳴クラスを確認する

```python
import b13phase_v075 as b13

BASE      = b13.BASE        # 3120
COS_TABLE = b13.COS_TABLE   # 整数コサインテーブル
N         = BASE

# 3種類の純音を生成
sig_full  = [COS_TABLE[(1  * k) % BASE] for k in range(N)]  # f₀=1  → full
sig_pent  = [COS_TABLE[(5  * k) % BASE] for k in range(N)]  # f₀=5  → pentagon
sig_z13   = [COS_TABLE[(13 * k) % BASE] for k in range(N)]  # f₀=13 → z13

for name, sig in [('f₀=1', sig_full), ('f₀=5', sig_pent), ('f₀=13', sig_z13)]:
    sp = b13.signal_to_fullerene(sig, max_harmonics=13)
    rc = b13.resonance_class(sp, max_harmonics=13)
    print(f"{name}: {rc['label']:8s}  z5_excess={rc['z5_excess']:+.4f}")
```

```
f₀=1:  full      z5_excess=-0.1538
f₀=5:  pentagon  z5_excess=+0.8462
f₀=13: z13       z5_excess=-0.1538
```

### WAV ファイルを解析する

```python
import b13phase_v075 as b13

# Layer 1: 低域（mh=13）— 高速、楽器種別の粗い分類
frames, mh, sr = b13.wav_to_frames('flute.wav', layer=1)

# Layer 2: 基音帯（mh=w_f0×2）— f0 既知時の正確な判定
frames, mh, sr = b13.wav_to_frames('flute.wav', f0_hz=1046.5, layer=2)

summary = b13.timeseries_summary(frames)
print(summary['label_counts'])   # {'full': 2, 'mixed': 9} など
print(summary['label_sequence']) # フレームごとのラベル列
```

### 成分分解

```python
sp  = b13.signal_to_fullerene(signal, max_harmonics=592)
dec = b13.spectrum_decompose(sp)
print(f"pentagon={dec['pentagon_ratio']*100:.1f}%  "
      f"z13={dec['z13_ratio']*100:.1f}%  "
      f"full={dec['full_ratio']*100:.1f}%")
# pentagon=5.6%  z13=2.5%  full=91.8%  (フルート C6)
```

### 成分分離（引き算）

```python
# pentagon 成分を除去する
sp_pent     = b13.spectrum_mask_z5(sp)
sp_residual = b13.spectrum_subtract(sp, sp_pent)

# Z₅=0 行の残留パワー = 0（完全除去）
residual = sum(b13.mag_sq(*sp_residual[0][m13]) for m13 in range(13))
print(residual)  # 0
```

---

## 主要 API

### `signal_to_fullerene(signal, max_harmonics=13)`

整数信号を Z₅×Z₁₃ の 65 点フラーレン座標に変換します。

| 引数 | 型 | 説明 |
|:--|:--|:--|
| `signal` | `List[int]` | 整数値の信号列（BASE スケール推奨） |
| `max_harmonics` | `int` | 展開する倍音数。実データでは `w_f0*2` を使う |

戻り値: `spectrum[m5][m13] = (cx, cy)` — 完全整数、m5∈{0..4}、m13∈{0..12}

---

### `resonance_class(spectrum, max_harmonics=13)`

フラーレンスペクトルから共鳴クラスを判定します。

**z5_excess 方式**（背景レベルを差し引いた超過量で判定）:

```
z5_background = (max_harmonics // 5) / max_harmonics
z5_excess     = z5_node_ratio - z5_background

z5_excess > +0.70 → pentagon
z5_excess < -0.05 かつ z13_excess < -0.05 → full
それ以外 → mixed
```

⚠️ `max_harmonics` は `signal_to_fullerene()` に渡した値と一致させること。mh=13 のとき Z₅=0 行の均等分布期待値は 2/13≈15.4% となり、旧来の「z5 < 0.1 → full」では full が必ず mixed と誤判定されます。

戻り値の主要フィールド: `label`, `z5_node_ratio`, `z13_node_ratio`, `z5_excess`, `z13_excess`, `total`

---

### `wav_to_frames(path, f0_hz=None, layer=1)`

WAV ファイルをフラーレン時系列フレームに変換します。

| `layer` | `max_harmonics` | 用途 |
|:-:|:--|:--|
| 1 | 13 | 低域（〜46Hz）の楽器種別分類。高速・固定コスト |
| 2 | `w_f0 × 2` | 基音帯の正確な判定。`f0_hz` 必須 |

`w_f0 = round(BASE × f0_hz / sr)` — Layer 2 では `w_f0` だけでなく `w_f0×2` まで必要（スペクトル漏れの相殺のため）。

---

### 成分分離ユーティリティ

| 関数 | 説明 |
|:--|:--|
| `spectrum_mask_z5(sp)` | Z₅=0 行のみ残す（pentagon 成分） |
| `spectrum_mask_z13(sp)` | Z₁₃=0 列のみ残す（z13 成分） |
| `spectrum_decompose(sp)` | pentagon / z13 / dc / full 比率を返す |
| `spectrum_subtract(sp_a, sp_b)` | スペクトル差分（完全整数） |
| `spectrum_add(sp_a, sp_b)` | スペクトル和（完全整数） |

---

### その他

| 関数 | 説明 |
|:--|:--|
| `winding_int(signal, w)` | 周波数 w での整数巻き取り → `(cx, cy)` |
| `mag_sq(cx, cy)` | `cx² + cy²`（整数） |
| `resonance_label(f0)` | `gcd(f0%65, 65)` から共鳴クラス名を返す |
| `resonance_orbit_size(f0)` | Z₆₅ 上の軌道サイズ |
| `analyze_timeseries(signal, ...)` | 時系列フレーム分析 |
| `timeseries_summary(frames)` | フレーム集計・遷移抽出 |
| `read_wav_scaled(path)` | WAV → BASE スケール整数列 |
| `estimate_f0_fft(samples, sr)` | numpy.fft で基音推定 |

---

## 実データ検証結果（C6, sr=11025Hz）

フルート・バイオリン・トランペット・ピアノの C6（1046.5Hz）で検証しました。

### Layer 2（mh=592）での成分比率

| 楽器 | pentagon% | z13% | full% |
|:--|:-:|:-:|:-:|
| flute | 5.6% | 2.5% | **91.8%** |
| violin | 9.2% | 1.5% | **90.7%** |
| trumpet | 8.7% | 0.3% | **91.0%** |
| piano | 3.0% | 1.0% | **96.0%** |

### 4楽器合成・逐次分離

4楽器を加算合成したスペクトルから既知スペクトルを1つずつ引いていくと、最後の残渣パワーが完全にゼロになることを確認しました（非ゼロセル 0/65）。

```
合成(4楽器)  3.530e+18  100.0%
- trumpet    3.529e+18   99.98%
- violin     2.064e+18   58.47%
- flute      1.692e+15    0.05%
- piano      0.000e+00    0.00%  ← 完全消滅 ✓
```

整数演算のため、位相も含めて誤差ゼロで分離できます。

---

## N≠BASE の実データへの対応

実データのサンプリングレートは BASE=3120 と一致しません。

```
C6(1046.5Hz), sr=11025Hz の場合:
  1周期 = 11025 / 1046.5 = 10.535... サンプル（非整数・Sturmian 構造）
  → w=296(C6) の隣 w=295(pentagon 座席) にスペクトル漏れが発生
  → max_harmonics=w_f0 のみだと pentagon と誤判定
  → max_harmonics=w_f0×2 で漏れが相殺 → full ✓
```

`frame_size = BASE = 3120` は固定。`max_harmonics` の選択で解析レイヤーを使い分けます。

---

## モジュール構成

```
b13phase_v075/
├── __init__.py              エクスポート定義（v0.8.0）
├── constants.py             BASE, STEP_Z5/Z13, コセット定義
├── level0_table.py          COS_TABLE, SIN_TABLE（3120エントリ）
├── fullerene_transform.py   signal_to_fullerene, resonance_class 等
├── audio_io.py              read_wav_scaled, wav_to_frames 等
├── fib_drive.py             FibDriveEngine（Fibonacci/Pell駆動）
├── fractal_engine.py        閾値交差・発火応答
├── z13_structure.py         Z₁₃コセット代数
├── truncated_icosahedron.py 頂点座標・中心角
├── great_circle.py          大円解析・φ算術
├── cell.py                  B13Cell, B13Fractal
├── connection.py            隣接グラフ・伝播
├── phase_digits.py          任意精度BASE進演算
├── phase_packed_u64.py      packed u64演算
└── evaluator_proto.py       多桁位相評価器
```

---

## ライセンス

MIT License

---

## 関連

- Zenn 記事: 「40×40×40格子からフラクタル・フラーレンへ — 波を「構造で巻き取る」」
- 引継書: `Claude-2026-03-28-v80.md`（開発経緯・設計判断の記録）
