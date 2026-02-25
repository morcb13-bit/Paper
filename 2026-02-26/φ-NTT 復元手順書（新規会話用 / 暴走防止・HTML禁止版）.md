# φ-NTT 復元手順書（新規会話用 / 暴走防止・HTML禁止版）

---

## 0) 禁止事項（最重要・最初に読む）

- **HTMLを新規作成・改造・提案しない**（phi_conv_demo.html / phi_ntt_demo.html に触れない）
- **新しいデモ・UI・可視化拡張をしない**
- **B=4,5 の再計算・大量走査・新実験をしない**
- **このセッションでやること**：論文 md → `.tex` 化（arXiv最小テンプレ）および文面整形のみ

---

## 1) 正本（freeze point）

```
統合md:     https://github.com/morcb13-bit/Artemis-B13-Archive/tree/main/tutorial
            ファイル: phi-ntt-paper.md @ [PAPER_COMMIT_SHA]

参照実装:   同リポジトリ内 phi_carry_free.py @ [IMPL_COMMIT_SHA]
```

※ SHA は未記入の場合でも作業は進められる。後から埋めれば良い。

---

## 2) 復元に必要な部品（v2 が基礎）

| バージョン | 収録内容 |
|-----------|---------|
| **v2** | 定義・perm10・80^B分解・phi_carry_free.py 全文（数学の核） |
| v4〜v6 | 理論確定：TUT最上位根拠・TT…T射影解釈・UUU≈0・準定理 |
| v7〜v12 | 論文 Section 1〜7 draft（統合md に集約済み） |

v2 の `phi_carry_free.py` さえあれば実装は完全復元可能。論文本文は統合md が正本。

---

## 3) 動作確認（30秒・これだけ）

```python
from phi_carry_free import phi_conv_carry_free
from random import seed, randint
seed(42)
B, N = 3, 1000
x = [(randint(-5,5), randint(-1,1)) for _ in range(N)]
h = [(0,0)]*N; h[0] = (1,0)
y = phi_conv_carry_free(x, h, B)
print("delta conv B=3:", "✓" if all(y[i]==x[i] for i in range(N)) else "✗")
```

✓ が出れば実装は正常。**それ以上の検証・再計測は禁止。**

---

## 4) 次タスク（このセッションでやること）

1. **統合md → `.tex` 変換**（arXiv最小テンプレ）
   - タイトル・著者・abstract・Section 1〜7・Appendix A〜C
   - 数式は既存 draft の範囲のみ（新規導出しない）
   - 図なし、Table 1 / Table 2 のみ（v11 の表をそのまま）
2. **用語統一チェック**（carry-free / Z[φ] / quasi-theorem の表記揺れ）
3. **commit SHA 埋め**（作業完了後にヘッダへ記入）

---

## 5) 論文の現在地（完成状態）

| Section | 内容 | 語数 | 状態 |
|---------|------|------|------|
| 1. Introduction | 背景・貢献・構成 | 310 | ✅ |
| 2. φ-NTT on Z₁₀^B | 群構造・定義・畳み込み定理 | 530 | ✅ |
| 3. Wavelet Interpretation | T/U演算子・二分木・集中メカニズム | 580 | ✅ |
| 4. Quasi-Theorem | 仮定・主張・証明スケッチ・実測整合 | 650 | ✅ |
| 5. Numerical Experiments | Setup・Table1/2・再現精度・計算構造 | 560 | ✅ |
| 6. Discussion | 未解決3点(P1〜P3)・既存手法比較・限界 | 530 | ✅ |
| 7. Conclusion | 貢献3点・未解決・展望 | 230 | ✅ |
| Appendix A〜C | cos5_zphi / perm10 / 80^B分解 | — | ✅ |
| **合計** | | **≈3,390語** | **✅** |

Short paper（8〜12ページ）として適切な分量。

---

## 6) 未解決問題（将来課題・このセッションでは触らない）

| ID | 内容 |
|----|------|
| P1 | T/U 直交性の代数的証明 |
| P2 | 定数 C と ε の閉形式 |
| P3 | carry あり循環畳み込みとの接続（overlap-add） |

---

*作成: 2026-02-26 / 引継書 v12 より*
