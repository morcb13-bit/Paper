---
title: "第6章　evaluator API リファレンス"
---

本章は B13 パッケージの完全なAPIリファレンスである。`evaluate()`・`B13Vertex`・補助関数の仕様を網羅する。

---

## 6.1　B13Vertex ── 出力の型

`evaluate()` が返す `B13Vertex` は Python の `NamedTuple` であり、全フィールドが整数である。

```python
class B13Vertex(NamedTuple):
    k: int          # 入力: Z₁₂ 五角形インデックス (0–11)
    j: int          # 入力: Z₅  五角形内頂点番号  (0–4)
    n: int          # 入力: フラクタルレベル      (0以上)
    az_idx:    int  # 方位角インデックス in Z_3120 (厳密)
    cos_az:    int  # COS_TABLE[az_idx]           (厳密)
    sin_az:    int  # SIN_TABLE[az_idx]           (厳密)
    cos_theta: int  # BASE × cos(極角)            (誤差 ≤ 0.02%)
    sin_theta: int  # BASE × sin(極角)            (誤差 ≤ 0.02%)
    r_sq_rat:  int  # 4R²ₙ の有理部  in Z[φ]     (厳密)
    r_sq_phi:  int  # 4R²ₙ の φ係数 in Z[φ]     (厳密)
```

各フィールドのスケールは BASE = 3120 である。浮動小数点座標に変換するには `to_float_coords()` を使う（6.4節）。

![B13Vertex フィールド構成図](/images/fig_ch6_vertex_fields.png)
*図 6-1　B13Vertex の全フィールド。色は精度保証の種別を示す。上段は `evaluate(0, 0, 0)` の実際の出力値。*

---

## 6.2　evaluate() ── メイン関数

### シグネチャ

```python
def evaluate(k: int, j: int, n: int = 0) -> B13Vertex
```

### パラメータ

| パラメータ | 型 | 範囲 | 説明 |
|-----------|-----|------|------|
| `k` | `int` | 0〜11 | 五角形リングのインデックス。z 高さ（極角）を決定する |
| `j` | `int` | 0〜4  | 五角形内の頂点番号。方位角を決定する |
| `n` | `int` | 0以上 | フラクタルレベル。省略時は `n=0` |

範囲外の値を渡すと `ValueError` が送出される。

```python
evaluate(-1, 0, 0)   # ValueError: k must be in 0..11, got -1
evaluate(0, 5, 0)    # ValueError: j must be in 0..4, got 5
evaluate(0, 0, -1)   # ValueError: n must be ≥ 0, got -1
```

### 基本的な使い方

```python
from b13phase import evaluate

# 単一頂点の取得
v = evaluate(k=0, j=0, n=0)
print(v)
# B13Vertex(k=0, j=0, n=0, az_idx=0, cos_az=3120, sin_az=0,
#           cos_theta=2915, sin_theta=1113, r_sq_rat=10, r_sq_phi=9)

# フィールドへのアクセス
print(v.az_idx)     # 0
print(v.cos_theta)  # 2915
print(v.r_sq_rat)   # 10
print(v.r_sq_phi)   # 9
```

### 60頂点の一括取得

```python
# レベル n=0 の全60頂点
all_vertices = [evaluate(k, j, n=0)
                for k in range(12)
                for j in range(5)]

# 頂点ラベル m でアクセスしたい場合
def find_by_m(m_target, n=0):
    for k in range(12):
        for j in range(5):
            if (5*k + 12*j) % 60 == m_target:
                return evaluate(k, j, n)
    raise ValueError(f"m={m_target} not found")

v = find_by_m(m_target=15, n=0)   # m=15 の頂点
```

### 複数レベルの走査

```python
# n=0,1,2 の全頂点（合計180頂点）
multi_level = [evaluate(k, j, n)
               for n in range(3)
               for k in range(12)
               for j in range(5)]

print(len(multi_level))   # 180
```

---

## 6.3　フィールド詳細

### 6.3.1　az_idx と cos_az / sin_az

`az_idx` は方位角を $\mathbb{Z}_{3120}$ 上の整数で表したものである。

$$
\text{az\_idx}(k, j, n) = \bigl((5k + 12j) \times 52 + n \times 312\bigr) \bmod 3120
$$

`cos_az` と `sin_az` はそれぞれ `COS_TABLE[az_idx]` と `SIN_TABLE[az_idx]` の値であり、テーブル参照による厳密整数である。スケールは BASE = 3120。

```python
# az_idx から実際の方位角（度）を計算
v = evaluate(k=1, j=0, n=0)
angle_deg = v.az_idx / 3120 * 360
print(f"az_idx={v.az_idx}, angle={angle_deg:.1f}°")  # 30.0°

# cos/sin の実数値
import math
cos_true = math.cos(math.radians(angle_deg))
print(f"cos_az/BASE = {v.cos_az/3120:.6f}")   # 0.866025
print(f"cos 真値    = {cos_true:.6f}")         # 0.866025
```

### 6.3.2　cos_theta と sin_theta

`cos_theta` と `sin_theta` は極角 $\theta$ の整数近似である（スケール BASE = 3120）。バンドによって4種類の定数値しか取らない。

| バンド | k | `cos_theta` | `sin_theta` | $\cos\theta$ 真値 |
|--------|---|:-----------:|:-----------:|:-----------------:|
| 北極   | 0 | +2915       | 1113        | $\varphi/\sqrt{3} \approx 0.9343$ |
| 上帯   | 1–5 | +1395     | 2791        | $1/\sqrt{5} \approx 0.4472$ |
| 下帯   | 6–10 | −1395    | 2791        | $-1/\sqrt{5} \approx -0.4472$ |
| 南極   | 11 | −2915      | 1113        | $-\varphi/\sqrt{3} \approx -0.9343$ |

`sin_theta` は常に正であることに注意する。

### 6.3.3　r_sq_rat と r_sq_phi

`r_sq_rat` と `r_sq_phi` は半径の2乗 $R^2_n$ を $\mathbb{Z}[\varphi]$ で表したものである。

$$
4R^2_n = \text{r\_sq\_rat} + \text{r\_sq\_phi} \times \varphi
$$

更新式：レベルが1増えるごとに $(\text{rat},\ \text{phi}) \to (\text{rat}+\text{phi},\ \text{rat}+2\times\text{phi})$。

```python
from b13phase.evaluator import r_sq_level
import math

phi = (1 + math.sqrt(5)) / 2

for n in range(5):
    a, b = r_sq_level(n)
    R = math.sqrt((a + b*phi) / 4)
    print(f"n={n}: 4R²={a}+{b}φ,  R={R:.4f}")

# n=0: 4R²=10+9φ,   R=2.4780
# n=1: 4R²=19+28φ,  R=4.0095
# n=2: 4R²=47+75φ,  R=6.4875
# n=3: 4R²=122+197φ, R=10.4971
# n=4: 4R²=319+516φ, R=16.9846
```

---

## 6.4　補助関数

### to_float_coords()

`B13Vertex` を浮動小数点座標 $(x, y, z)$ に変換する。検証・可視化用途向けである。

```python
from b13phase.evaluator import to_float_coords

v = evaluate(k=1, j=2, n=0)
x, y, z = to_float_coords(v)
print(f"x={x:.4f}, y={y:.4f}, z={z:.4f}")
# x=1.9198, y=-0.0617, z=1.1079
```

復元式は以下のとおりである。

$$
R_n = \sqrt{\frac{\text{r\_sq\_rat} + \text{r\_sq\_phi} \cdot \varphi}{4}}, \quad
x = R_n \cdot \frac{\text{sin\_theta}}{\text{BASE}} \cdot \frac{\text{cos\_az}}{\text{BASE}}
$$

$$
y = R_n \cdot \frac{\text{sin\_theta}}{\text{BASE}} \cdot \frac{\text{sin\_az}}{\text{BASE}}, \quad
z = R_n \cdot \frac{\text{cos\_theta}}{\text{BASE}}
$$

### verify_all()

4項目の自己検証を実行し、結果を辞書で返す。

```python
from b13phase.evaluator import verify_all

result = verify_all(n_levels=6)
```

| キー | 型 | 意味 |
|------|-----|------|
| `az_ok` | `bool` | `az_idx` 公式と参照実装の一致 |
| `r2_ok` | `bool` | $R^2$ の $\mathbb{Z}[\varphi]$ 厳密性（誤差 $< 10^{-8}$） |
| `theta_ok` | `bool` | 極角の全レベル不変性（誤差 $< 0.1\%$） |
| `distinct_ok` | `bool` | 全レベルで60頂点が区別可能 |
| `max_r2_err` | `float` | $R^2$ の最大絶対誤差（数値計算誤差） |
| `max_coord_err` | `float` | 極角の最大誤差 |

```python
# 正常な出力例
# {
#   'az_ok':         True,
#   'r2_ok':         True,
#   'theta_ok':      True,
#   'distinct_ok':   True,
#   'max_r2_err':    2.27e-13,
#   'max_coord_err': 1.85e-04
# }
```

### r_sq_level()

レベル $n$ における $4R^2_n$ を $\mathbb{Z}[\varphi]$ の整数ペアで返す。

```python
from b13phase.evaluator import r_sq_level

a, b = r_sq_level(n=3)
print(f"4R²₃ = {a} + {b}φ")   # 4R²₃ = 122 + 197φ
```

### band() と theta_ints()

```python
from b13phase.evaluator import band, theta_ints

print(band(0))          # 0  (北極)
print(band(3))          # 1  (上帯)
print(band(8))          # 2  (下帯)
print(band(11))         # 3  (南極)

print(theta_ints(0))    # (2915, 1113)   (cos_theta, sin_theta)
print(theta_ints(3))    # (1395, 2791)
```

---

## 6.5　精度の保証と限界

![verify_all の検証結果](/images/fig_ch6_verify.png)
*図 6-2　左：$R^2$ の $\mathbb{Z}[\varphi]$ 値と真値の比較（対数軸、誤差 $< 10^{-12}$）。右：各バンドの `cos_theta` 近似誤差（全バンドで 0.02% 以下）。*

### 保証する精度

| フィールド | 精度 | 根拠 |
|-----------|------|------|
| `az_idx` | 厳密（誤差なし） | 整数 mod 演算 |
| `cos_az`, `sin_az` | 厳密（誤差なし） | テーブル参照 |
| `r_sq_rat`, `r_sq_phi` | 厳密（$< 10^{-12}$） | $\mathbb{Z}[\varphi]$ 整数演算 |
| `cos_theta`, `sin_theta` | $\leq 0.02\%$ | 整数丸め |

### to_float_coords() の精度

`to_float_coords()` は内部で浮動小数点演算（平方根・乗算）を使うため、Python の `float` 精度（約 $10^{-15}$）に依存する。`r_sq_rat` / `r_sq_phi` の整数値自体は厳密であり、誤差は変換時にのみ発生する。

---

## 6.6　パッケージ構成

`b13phase` パッケージのファイル構成と各モジュールの役割を示す。

```
b13phase/
  __init__.py               公開 API（evaluate, B13Vertex をエクスポート）
  constants.py              BASE=3120, Z13_STAR, H_SET 等の定数
  level0_table.py           COS_TABLE, SIN_TABLE（各3120エントリ）
  evaluator.py              ★ 本番実装（v60 確定版）
  evaluator_proto.py        旧プロト（互換性のため残存）
  z13_structure.py          Z₁₃ 代数、Fib/Pell 4H 回避
  exact_icosahedron.py      切端20面体 60頂点の Z[φ] 座標
  fractal_vertex_recurrence.py  フラクタル漸化式（参照実装）
  great_circle.py           大円閉合定理
  phase_digits.py           多桁位相算術
  phase_packed_u64.py       u64 パック算術
```

通常の利用では `evaluate` と `B13Vertex` のみをインポートすれば十分である。

---

## 6.7　クイックリファレンス

```python
from b13phase import evaluate, B13Vertex
from b13phase.evaluator import (
    to_float_coords,   # → (x, y, z) float
    verify_all,        # → dict  4項目検証
    r_sq_level,        # → (a, b)  4R²ₙ = a+bφ
    band,              # → 0/1/2/3
    theta_ints,        # → (cos_theta, sin_theta)
    az_idx,            # → int  方位角インデックス
)
from b13phase.constants import BASE   # = 3120

# 典型的なパターン
v = evaluate(k, j, n)          # 頂点レコード取得
x, y, z = to_float_coords(v)   # 浮動小数点変換
result = verify_all(n_levels=6) # 自己検証
assert all(result[k] for k in ['az_ok','r2_ok','theta_ok','distinct_ok'])
```

---

## おわりに

本書では B13 座標系の全体を以下の流れで解説した。

1. **第1章**：`evaluate()` の動作と60頂点の整数表現という全体像。
2. **第2章**：頂点ラベル $m = (5k+12j) \bmod 60$ の全単射性、バンド分類、対蹠点。
3. **第3章**：$\mathbb{Z}_{13}^*$ のコセット構造と Fib/Pell の 4H 回避（理論背景）。
4. **第4章**：`az_idx` 公式の導出、$\mathbb{Z}[\varphi]$ 算術、極角不変性という3つの整数化の核心。
5. **第5章**：大円閉合定理と BASE=3120 の幾何的根拠（理論背景）。
6. **第6章**：本章。API の完全仕様。

B13 の設計思想は一貫している。**角度・半径・極角の3つをすべて整数演算に帰着させ、誤差を定数オーダーに制御する**。この整合性が `verify_all()` によって任意のレベルで確認できる点が、実装の信頼性の核心である。

:::message
**今後の課題**　引継書 v60 時点で未確定の項目として、Fib/Pell 4H 回避の完全な代数的証明が残っている。数値的には $n=0 \sim 200$ で確認済みである。
:::
