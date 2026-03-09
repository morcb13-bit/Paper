# VALIDATION PROTOCOL — coneFFT Verification Instrument v0.2.1

第三者が同じ手順で検証を再現できるよう、観測手順・判定基準・解釈ガイドを記述する。

---

## 前提

- `cone_spectrum_demo_v0.2.1.html` をブラウザで開く
- ウィンドウ関数: **Hann**（デフォルト）
- 入力モード: **TEST**（テスト信号）
- N = 1024、SR = 44100Hz

---

## Step 1 — B13 Resonance テスト

**目的:** cone observer が B13構造に再現よく応答するか確認する。

### 手順

1. 信号セレクタを `B13 Resonance` に設定
2. `RESET` ボタンを押してセッションをリセット
3. そのまま **20 frame 以上**（約5〜10秒）観測する
4. 右パネルの以下を記録する：

| 記録項目 | 場所 |
|:--|:--|
| B13再現率 | false positive guard → B13再現率 |
| cone B13 ratio | B13 band energy → cone B13 ratio |
| fft B13 ratio | B13 band energy → fft B13 ratio |
| cone/fft 比 | B13 band energy → cone/fft |
| max Δ dB | DIFF quantification → max Δ dB |
| Φ_d CV | Φ_d session statistics → CV |

### 合格基準（暫定）

| 指標 | 目安 |
|:--|:--|
| B13再現率 | ≥ 80% |
| cone/fft 比 | ≥ 1.5× |
| Φ_d CV | < 0.1 |

---

## Step 2 — White Noise テスト

**目的:** cone observer が任意雑音に対して B13帯域への偽反応を示さないか確認する。

### 手順

1. 信号セレクタを `White Noise` に切り替え（RESETしない）
2. **20 frame 以上**観測する
3. 右パネルの以下を記録する：

| 記録項目 | 場所 |
|:--|:--|
| noise偽陽性率 | false positive guard → noise偽陽性率 |
| ノイズ反応 | false positive guard → ノイズ反応 |
| cone B13 ratio (noise時) | B13 band energy → cone B13 ratio |

### 合格基準（暫定）

| 指標 | 目安 |
|:--|:--|
| noise偽陽性率 | ≤ 5% |
| ノイズ反応 | `✓ 無反応` |

---

## Step 3 — 標準信号との比較

**目的:** B13以外の信号で cone observer が無意味に反応しないか確認する。

### 手順

各信号について 10 frame 以上観測し、cone B13 ratio を記録する：

| 信号 | 期待される cone B13 ratio |
|:--|:--|
| Sine 440Hz | 低い（単純正弦波への反応は限定的） |
| Square 440Hz | 中程度（高調波を含むが B13構造ではない） |
| Fibonacci Chirp | 中〜高（フィボナッチ構造を含む） |
| B13 Resonance | 高い（設計上） |
| White Noise | 低い（偽陽性ガード） |

---

## VERDICT の読み方

右パネル下部の `SESSION VERDICT` は自動判定を表示する。

| 表示 | 意味 |
|:--|:--|
| `✓ B13再現: XX%` | B13信号でのピーク再現率が合格 |
| `△ B13再現: XX%` | 再現率が閾値未満（80%以下） |
| `✓ 偽陽性: XX%` | noise での誤反応率が合格 |
| `⚠ 偽陽性: XX%` | noise での誤反応率が高い |
| `✓ B13帯域 N.NN×` | cone/fft 比が 2.0 以上 |
| `△ B13微優位 N.NN×` | cone/fft 比が 1.2〜2.0 |
| `✓ Φ_d 安定 CV=X.XXX` | 変動係数が 0.1 未満 |

---

## 数値ログの記録方法

現時点では手動記録を推奨する（自動ログ出力は未実装）。

推奨記録フォーマット：

```
日時: YYYY-MM-DD HH:MM
バージョン: v0.2.1
ウィンドウ: Hann

[B13 Resonance / 30frame]
  B13再現率:      XX.X%
  cone B13 ratio: X.XX%
  fft B13 ratio:  X.XX%
  cone/fft:       X.XX×
  max Δ dB:       XX.XX dB @ bin XXX
  Φ_d CV:         0.XXX

[White Noise / 25frame]
  noise偽陽性率:  X.X%
  cone B13 ratio: X.XX%

VERDICT: [OK / WARN / NG]
備考:
```

---

## 既知の制約

- `dftMag()` は `step=4` 間引き計算のため、FFT側との精度差がある
- cone spectrum は DFT ベースのハイブリッド検出器（純粋な R₆ 再帰スペクトルではない）
- B13再現率の「合格」は仮説の証明ではなく、**観測の再現性**を意味する

---

## 改訂履歴

| バージョン | 変更内容 |
|:--|:--|
| v0.2.1 | B13帯域評価を Hz座標系に統一、drawDiff バグ修正、FFTマーカー修正 |
| v0.2 | 右パネル追加（4指標）、Φ_d統計、偽陽性ガード |
| v0.1 | 初版（FFT/cone/DIFF 3段表示） |
