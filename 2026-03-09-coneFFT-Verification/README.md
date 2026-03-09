# coneFFT Verification Instrument v0.2.1

R₆ / N-tree / Φ_d に基づく cone 系応答を、標準FFTと並列比較するためのブラウザ実装。

本実装は「FFTの完全代替」を主張するものではない。  
B13仮説に関わる帯域応答と偽陽性の有無を観測・定量化するための**検証器**として設計されている。

> **本リポジトリは、coneFFT の完成理論を主張するものではなく、観測と検証を行うための実装装置と技術メモを公開するものである。**

---

## これは何か

標準FFTが「どの周波数に何があるか」を測る装置とすれば、  
coneFFT は「どの方向・どの構造に cone 系の応答があるか」を測る観測器である。

実装の核心は3つの演算子：

| 演算子 | 定義 | 役割 |
|:--|:--|:--|
| R₆ | `x(p+d)−x(p-d)+x(p+2d)−x(p+3d)+x(p+5d)−x(p+8d)` | フィボナッチ差分 |
| N(a,b) | `(a+b, a-b)` | 双対ノード（cone展開） |
| Φ_d | `Σ R₆(p)·R₆(p+2d) / n` | 層不変量 |

---

## 何ができるか

- **FFT / coneFFT 並列スペクトル表示**
- **DIFF表示** — cone優位帯域と FFT優位帯域の可視化
- **B13帯域エネルギー評価** — 300Hz / 1000Hz ±80Hz（実注入周波数ベース）
- **Φ_d セッション統計** — 平均・分散・CV・スパークライン
- **偽陽性ガード** — noise信号でのB13帯域誤反応率
- **VERDICT** — セッション判定の自動要約

---

## 何が未解決か

以下は現時点で検証継続中であり、確定を主張しない：

- cone spectrum の出力は `R₆応答 + DFT` の合成（detector 実装）であり、最終的な coneFFT 変換定理そのものではない
- B13仮説（m=3,10 ピークの構造的必然性）は観測段階
- Φ_d の階層接続と多層再帰への拡張は未実装
- `t_* ≈ 48` の構造的必然性は未証明

確定済みの骨格は [TECHNICAL_NOTE.md](./TECHNICAL_NOTE.md) を参照。

---

## 動かし方

```
cone_spectrum_demo_v0.2.1.html をブラウザで開く
```

依存なし。単一HTMLファイル。ローカルで動作。  
マイク入力を使う場合は `https://` または `localhost` 環境が必要。

---

## 検証手順

詳細は [VALIDATION_PROTOCOL.md](./VALIDATION_PROTOCOL.md) を参照。

**最小手順：**

1. 信号を `B13 Resonance` に設定
2. 20 frame 以上観測 → `B13再現率` と `cone/fft` 帯域比を記録
3. 信号を `White Noise` に切り替え
4. 20 frame 以上観測 → `noise偽陽性率` を記録
5. `VERDICT` の判定を確認

---

## ファイル構成

```
README.md                          本ファイル
TECHNICAL_NOTE.md                  理論・実装・限界の詳細
VALIDATION_PROTOCOL.md             検証手順と判定基準
cone_spectrum_demo_v0.2.1.html     単体動作の実装本体
```

---

## 位置づけ

```
確定済みの骨格    → TECHNICAL_NOTE §3
実装済みの観測器  → cone_spectrum_demo_v0.2.1.html
未解決・検証中    → TECHNICAL_NOTE §7, §8
```

本実装の価値は「完全証明の宣言」ではなく、  
**他者が触れて検証できる観測系の提示**にある。
