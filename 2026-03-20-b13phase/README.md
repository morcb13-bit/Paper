# B13 技術解説書
## 切端20面体の整数フラクタル座標系

---

## 概要

B13 は、切端20面体（サッカーボール型多面体）の60頂点を**浮動小数点演算なし**で扱う整数座標系である。3つの整数 `(k, j, n)` を入力とし、全フィールドが整数の頂点レコード `B13Vertex` を返す。

```python
from b13phase import evaluate

v = evaluate(k=0, j=0, n=0)
# B13Vertex(k=0, j=0, n=0,
#           az_idx=0, cos_az=3120, sin_az=0,
#           cos_theta=2915, sin_theta=1113,
#           r_sq_rat=10, r_sq_phi=9)
```

---

## 本書の構成

| 章 | タイトル | 概要 |
|----|---------|------|
| 第1章 | B13 とは何か | 動作概要・最小動作例・4バンド構造・フラクタル展開の俯瞰 |
| 第2章 | 頂点ラベル m ∈ Z₆₀ | `m = (5k+12j) mod 60` の全単射性・対蹠点・コセット軸 |
| 第3章 | Z₁₃* とコセット | 乗法群の構造・H/2H/4H 分解・Fib/Pell の 4H 回避定理 ※ |
| 第4章 | フラクタル展開 | az_idx 公式・Z[φ] 算術・極角不変性の導出 |
| 第5章 | 大円閉合定理 | sin(θ_B)=2/3 の証明・BASE=3120 の幾何的根拠 ※ |
| 第6章 | evaluator API リファレンス | 全関数の仕様・精度保証・パッケージ構成 |

※ 第3章・第5章は理論背景章であり、初読スキップ可。  
実用優先ルート：**第1章 → 第2章 → 第4章 → 第6章**

---

## 対象読者・前提知識

- **対象**：実装寄りのエンジニア。Python が読める。
- **数学**：高校数学レベル。剰余演算・ベクトルの内積程度。群論・代数的整数論の知識は不要（本書内で必要最小限を補う）。

---

## パッケージのインストール

```bash
# リポジトリ直下に b13phase/ がある場合
pip install -e .

# または PYTHONPATH を通す
export PYTHONPATH=/path/to/b13phase_pkg/b13phase
```

```python
# 動作確認
from b13phase import evaluate
from b13phase.evaluator import verify_all

result = verify_all(n_levels=6)
assert result['az_ok'] and result['r2_ok'] and result['theta_ok'] and result['distinct_ok']
print("OK")
```

---

## B13 の設計原理

B13 の整数表現は以下の4点を保証する。

| フィールド | 精度 |
|-----------|------|
| `az_idx`, `cos_az`, `sin_az` | 厳密（誤差なし） |
| `r_sq_rat`, `r_sq_phi` | Z[φ] 上で厳密（誤差 < 10⁻¹²） |
| `cos_theta`, `sin_theta` | 整数近似（誤差 ≤ 0.02%） |
| 60頂点の識別 | 全レベルで常に区別可能 |

整数化の核心は3つである。

1. **方位角**：長さ 3120 の整数テーブル（`COS_TABLE` / `SIN_TABLE`）をインデックス参照。
2. **半径²**：黄金比の代数的整数環 Z[φ] 上の漸化式 `(a, b) → (a+b, a+2b)`。加算のみ。
3. **極角**：フラクタル展開は球の相似拡大であるため、θ はレベルによらず定数。

---

## パッケージ構成

```
b13phase/
  __init__.py                 公開 API
  constants.py                BASE=3120, Z13_STAR 等
  level0_table.py             COS_TABLE, SIN_TABLE（各3120エントリ）
  evaluator.py                ★ 本番実装（v60 確定版）
  z13_structure.py            Z₁₃ 代数・Fib/Pell 4H 回避
  exact_icosahedron.py        切端20面体 60頂点の Z[φ] 座標
  fractal_vertex_recurrence.py  フラクタル漸化式（参照実装）
  great_circle.py             大円閉合定理
  phase_digits.py             多桁位相算術
  phase_packed_u64.py         u64 パック算術
  evaluator_proto.py          旧プロト（互換性のため残存）
```

---

## 既知の未確定事項

本書（v60 時点）において以下の項目が未確定である。

- Fib/Pell 4H 回避の完全な代数的証明（数値的には n=0〜200 で確認済み）

---
