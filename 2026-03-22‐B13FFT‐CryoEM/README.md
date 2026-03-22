---
title: "B13-FFT × クライオEM 実データ解析：アポフェリチン殻構造の選択的抽出"
emoji: "🔬"
type: "tech"
topics: ["cryo-EM", "FFT", "Python", "数理生物学", "画像解析"]
published: false
---

# B13-FFT × クライオEM 実データ解析

**B13格子フィルタによるアポフェリチン殻構造の選択的抽出**

バージョン: v63 (2026-03-22)

---

## 概要

B13-FFT は整数論的構造 `BASE = 3120 = 2⁴ × 3 × 5 × 13` に基づくフーリエビン選択フィルタである。通常の低域/高域フィルタとは異なり、B13格子に対応する特定周波数ビンのみを元画像から差し引いた**残渣**を出力する。

アポフェリチン（外径 ~120 Å、殻厚 ~20 Å）に対してスケール 119 Å を設定すると：

- 殻構造シグナルが残渣に選択的に保持される
- 背景ノイズおよびCTFリング成分が除去される
- 動径プロファイルのゼロ交差が殻壁位置（40 Å / 60 Å）に収束する

本リポジトリは、実クライオEMデータ（EMPIAR-10146相当、1.34 Å/px）に対してこの解析を実行する完全な実装を提供する。

---

## 理論背景

### B13ビンの定義

B13格子点の方位角インデックスは以下で与えられる：

```
az_idx = ((5k + 12j) mod 60) × 52  mod 3120
    k ∈ {0, 1, ..., 11},  j ∈ {0, 1, ..., 4}
```

物理スケール $s$ [Å] とピクセルサイズ $d$ [Å/px] から `scale_factor = round(s / d)` を計算し、ビンをピクセル空間に写像する。アポフェリチン解析では $s = 119$ Å、$d = 1.34$ Å/px として `scale_factor = 89`。1240 px 行に対して **100ビン（全体の 8.1%）** が選択される。

### 残渣の計算

行・列を独立に処理し、その平均を残渣とする：

```
residual_row = image − IFFT(FFT_row  ×  B13_mask_row)
residual_col = image − IFFT(FFT_col  ×  B13_mask_col)
residual_2d  = (residual_row + residual_col) / 2
edge_product = |residual_row| × |residual_col|
```

`edge_product` は円形構造の境界を強調し、粒子検出の根拠となる。

### Fib/Pell の 4H 回避定理（v63 完全証明）

$\mathbb{Z}_{13}^*$ の部分群を $H = \{1, 5, 8, 12\}$、そのコセットを $2H = \{2, 3, 10, 11\}$、$4H = \{4, 6, 7, 9\}$ と定義する。

**定理:** すべての $n \geq 0$ に対して $F(n) \bmod 13 \notin 4H$、および $P(n) \bmod 13 \notin 4H$

**証明（3ステップ）:**

**Step 1** — 生成行列のスカラー冪（直接計算）

$$M^7 \equiv 8I \pmod{13}, \quad 8 \in H$$
$$N^7 \equiv 5I \pmod{13}, \quad 5 \in H$$

**Step 2** — $H$ の乗法閉包（数値確認済み）

$$H \times (H \cup 2H) = H \cup 2H \quad (4H \text{ を含まない})$$

**Step 3** — 帰納的完結

$n = 7q + r$ と書くと：

$$[F(n+1), F(n)]^T = M^n \cdot [1,0]^T = 8^q \cdot (M^r \cdot [1,0]^T)$$

$8^q \in H$（Step 1）、$M^r \cdot [1,0]^T \in \{0\} \cup H \cup 2H$（$r = 0, \ldots, 6$ の直接確認）、
$H \times (\{0\} \cup H \cup 2H) = \{0\} \cup H \cup 2H$（Step 2）$\Rightarrow F(n) \bmod 13 \notin 4H$ $\square$

---

## 環境構築

### 依存パッケージ

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install numpy scipy matplotlib mrcfile tqdm
```

Linux/macOS:

```bash
python -m venv .venv && source .venv/bin/activate
pip install numpy scipy matplotlib mrcfile tqdm
```

### ディレクトリ構成

```
.
├── run_b13_fft_v63.py              # メインスクリプト
├── sym_fft.py                      # 対称性ミスマッチ測定 (S4/S5)
├── fib_pell_4h_proof.py            # 4H回避定理検証スクリプト
├── packages/
│   ├── May08_07.30.46.bin.mrc      # クライオEMデータ（50フレーム, 297 MB）をダウンロードして貼り付ける
│   ├── b13phase_extracted/b13phase/
│   └── cone_fft_extracted/
└── data/results/                   # 出力先（自動生成）
```

---

## 実行手順

### クイックスタート（全フレーム平均・推奨）

```powershell
python run_b13_fft_v63.py `
    --input packages\May08_03.05.02.bin.mrc `
            packages\May08_07.30.46.bin.mrc `
    --all-frames `
    --outdir data\results\v63_latest
```

`--all-frames` を指定すると MRC 内の全 50 フレームを逐次積算して平均する（メモリ使用量を抑制しつつ SNR を向上させる）。

### 主なオプション

| オプション | デフォルト | 説明 |
|-----------|-----------|------|
| `--input` | *(必須)* | MRC / NPY（複数指定で平均） |
| `--pixel-size` | `1.34` | ピクセルサイズ [Å/px] |
| `--scale` | `119.0` | ターゲットスケール [Å] |
| `--outer-r` | `60.0` | 粒子外径 [Å] |
| `--inner-r` | `40.0` | 粒子内径 [Å] |
| `--all-frames` | off | 全フレーム平均（推奨） |
| `--no-plot` | off | 可視化スキップ |

### 処理ステップ

```
[1/4] フレーム読み込み    → raw_avg（フレーム平均画像）
[2/4] B13-FFT残渣計算    → residual_2d.npy, edge_product.npy
[3/4] 粒子候補検出        → particle_centers.npy, shell_scores.npy,
                            mean_radial_profile.npy
[4/4] 可視化              → b13_result.png（6パネル）
```

---

## 実験結果（v63）

### 解析パラメータ

| パラメータ | 値 |
|-----------|-----|
| 画像サイズ | 1200 × 1240 px |
| ピクセルサイズ | 1.34 Å/px |
| scale_factor | 89 |
| outer_r / inner_r | 60 Å / 40 Å |
| フレーム数 | 各50（全フレーム平均） |

### 定量結果

| 指標 | 値 |
|------|-----|
| 解析粒子数 | 91 / 101（端部除外後） |
| 殻スコア mean ± std | 0.204 ± 0.121 |
| 殻スコア最大値 | 0.585（cy=939, cx=497） |
| 動径ゼロ交差（内壁） | 41 Å（設計値 40 Å、誤差 1 Å） |
| 動径ゼロ交差（外壁） | 58 Å（設計値 60 Å、誤差 2 Å） |
| 対称性ミスマッチ S4/S5 | 1.325（> 1 → O群優勢） |

### CTFリング仮説の検証結果（v63 棄却）

残渣パワーピーク（54〜60〜115 Å）はアポフェリチンの既知構造寸法と系統的に整合し、CTF零点とは対応しない。B13-FFT は殻の実構造シグナルを保持しており、CTFリング除去フィルタではない。

| パワーピーク | 対応構造 |
|-------------|---------|
| 54〜55 Å（最大） | 殻中央半径 (60+40)/2 ≈ 50 Å |
| 60 Å | 外壁半径（完全一致） |
| 115 Å | B13スケール 119 Å（差 4 Å） |

### 対称性ミスマッチ（S4/S5 = 1.325）

B13格子はIh群（5回対称）を前提に設計されているが、アポフェリチンはO群（4回対称）を持つ。角度方向FFTによる測定：

| 対称成分 | 平均パワー比 ± std |
|---------|-----------------|
| 4回対称（O群特有） | 0.1161 ± 0.0992 |
| 5回対称（Ih群特有） | 0.0876 ± 0.0777 |
| **S4/S5** | **1.325** |

S4/S5 > 1 はB13残渣がO群の4回対称性を保持していることの定量的証拠である。

---

## 出力ファイル

| ファイル | 型 | 内容 |
|---------|-----|------|
| `residual_2d.npy` | float32 (H, W) | B13-FFT残渣マップ（符号付き） |
| `edge_product.npy` | float32 (H, W) | \|row_res\| × \|col_res\| |
| `particle_centers.npy` | int (N, 2) | 検出粒子座標 [(cy, cx)] |
| `shell_scores.npy` | float (N, 3) | [(score, cy, cx)] 降順 |
| `mean_radial_profile.npy` | float (R,) | 全粒子平均残渣プロファイル |
| `b13_result.png` | PNG 120 dpi | 6パネル可視化結果 |

---

## B13パッケージ基本定数

| 定数 | 値 |
|------|-----|
| BASE | 3120 = 2⁴ × 3 × 5 × 13 |
| scale_factor（119 Å） | 89 = round(119 / 1.34) |
| B13ビン数（1240 px） | 100ビン = 8.1% |
| H | {1, 5, 8, 12}（Z₁₃* の部分群） |
| 2H | {2, 3, 10, 11} |
| 4H | {4, 6, 7, 9}（Fib/Pellの禁止域） |

---

## 一行結論

**B13残渣はCTFではなくアポフェリチン殻構造シグナルを選択的に抽出し、O群4回対称性（S4/S5 = 1.325）を保持している。Fib/Pellの4H回避は $M^7 \equiv 8I \pmod{13}$ で代数的に証明された。**
