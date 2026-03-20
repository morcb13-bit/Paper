---
title: "第4章　フラクタル展開"
---

本章は本書の数学的な核心である。B13 の「フラクタル」とは何か、なぜ整数演算だけで実現できるのかを、`az_idx` 公式・$\mathbb{Z}[\varphi]$ 算術・極角不変性の3点から解説する。

---

## 4.1　BASE = 3120 の設計

すべての議論の出発点として、定数 BASE = 3120 の意味を確認する。

$$
\text{BASE} = 3120 = 2^4 \times 3 \times 5 \times 13
$$

この値は偶然ではなく、次の2つの条件を同時に満たすように選ばれている。

| 条件 | 式 | 値 |
|------|----|----|
| $\mathbb{Z}_{60}$ の1単位 | $\text{BASE} / 60$ | $= 52$（整数） |
| フラクタルレベル1単位（36°） | $\text{BASE} / 10$ | $= 312$（整数） |

因子 13 は $\mathbb{Z}_{13}^*$ の位数（第3章）と対応し、因子 $60 = 4 \times 3 \times 5$ は頂点ラベル $m \in \mathbb{Z}_{60}$ の周期と対応する。これらが BASE の中で整合する最小の値が 3120 である。

```python
from b13phase.constants import BASE

print(BASE)          # 3120
print(BASE // 60)    # 52  = Z₆₀ 1単位のステップ
print(BASE // 10)    # 312 = 1レベルあたりの方位角シフト（36°相当）
```

---

## 4.2　az_idx 公式の導出

方位角インデックス `az_idx` は、頂点の球面上での「経度」を整数で表したものである。

### 4.2.1　公式

$$
\text{az\_idx}(k, j, n) = \bigl(m(k,j) \times 52 + n \times 312\bigr) \bmod 3120
$$

ただし $m(k,j) = (5k + 12j) \bmod 60$ である。2つの項の意味は以下のとおりである。

- **$m \times 52$**：頂点ラベル $m \in \mathbb{Z}_{60}$ をBASEスケールに変換する。1ステップあたり $360° / 60 = 6°$ に対応する。
- **$n \times 312$**：レベルが1増えるごとに方位角を $312/3120 \times 360° = 36°$ だけシフトする。

### 4.2.2　なぜ 36° シフトするのか

切端20面体の5角形は72°（= 360°/5）ごとに並ぶ。フラクタル展開は正10角形（36° = 360°/10）の対称性を持つため、1レベルあたりのシフトが 36° となる。312 = BASE/10 はこの幾何的事実を整数で表現したものである。

```python
from b13phase.evaluator import az_idx

# 同じ頂点 (k=1, j=0) のレベル別 az_idx
for n in range(5):
    ai = az_idx(1, 0, n)
    angle = ai / 3120 * 360
    print(f"n={n}: az_idx={ai:4d}  ({angle:.1f}°)")

# n=0: az_idx= 260  (30.0°)
# n=1: az_idx= 572  (66.0°)   ← +36°
# n=2: az_idx= 884  (102.0°)  ← +36°
# n=3: az_idx=1196  (138.0°)  ← +36°
# n=4: az_idx=1508  (174.0°)  ← +36°
```

![az_idx の構造](/images/fig_ch4_az_idx.png)
*図 4-1　左：$n=0$ での `az_idx` と $m$ の線形関係（傾き 52）。右：同一頂点の `az_idx` がレベルごとに 312 ずつシフトする様子。*

### 4.2.3　整数 cos/sin テーブル

`az_idx` はそのまま `COS_TABLE` / `SIN_TABLE` のインデックスとして使われる。これらは BASE スケールの整数値を格納した長さ 3120 のテーブルである。

$$
\text{COS\_TABLE}[i] \approx \text{BASE} \times \cos\!\left(\frac{2\pi i}{\text{BASE}}\right), \quad i \in \{0, 1, \ldots, 3119\}
$$

テーブル参照は丸め誤差を生まない。`evaluate()` が返す `cos_az` と `sin_az` は厳密な整数値である。

```python
from b13phase.evaluator import evaluate

v = evaluate(k=1, j=0, n=0)
print(f"az_idx = {v.az_idx}")     # 260
print(f"cos_az = {v.cos_az}")     # COS_TABLE[260]（整数）
print(f"sin_az = {v.sin_az}")     # SIN_TABLE[260]（整数）

# 実際の角度と比較
import math
angle = v.az_idx / 3120 * 2 * math.pi
print(f"cos 理論値: {math.cos(angle)*3120:.2f}")  # ≈ cos_az の値
```

---

## 4.3　Z[φ] による半径の厳密表現

### 4.3.1　Z[φ] とは何か

$\mathbb{Z}[\varphi]$ は黄金比 $\varphi = (1+\sqrt{5})/2$ を含む代数的整数環である。任意の元は

$$
a + b\varphi \quad (a, b \in \mathbb{Z})
$$

の形で表される。加算・減算・乗算はすべて整数の演算に帰着する。

特に重要な関係式：

$$
\varphi^2 = \varphi + 1
$$

これにより、$\varphi^2$ による乗算が整数ペアの加算のみで実現できる。

### 4.3.2　R² の初期値

レベル $n=0$ における半径の2乗は

$$
R^2_0 = \frac{9\varphi + 10}{4} \quad \in \mathbb{Z}[\varphi]
$$

である。$\mathbb{Z}[\varphi]$ の元として $(a_0, b_0) = (10, 9)$ と表す（分母 4 は外に出す）。すなわち $4R^2_0 = 10 + 9\varphi$。

### 4.3.3　レベルごとの更新式

レベルが1増えると半径の2乗は $\varphi^2$ 倍になる。

$$
R^2_n = \varphi^{2n} \cdot R^2_0
$$

$\varphi^2 = \varphi + 1$ を使うと、$(a, b)$ で表された $\mathbb{Z}[\varphi]$ 元に $\varphi^2$ を掛ける操作は

$$
(a + b\varphi) \cdot \varphi^2 = (a + b\varphi)(\varphi + 1) = (a+b) + (a+2b)\varphi
$$

すなわち整数の更新式：

$$
\boxed{(a,\, b) \;\longrightarrow\; (a+b,\;\; a+2b)}
$$

加算のみで $\varphi^2$ 倍が実現する。浮動小数点乗算は一切不要である。

```python
from b13phase.evaluator import r_sq_level

for n in range(7):
    a, b = r_sq_level(n)   # 4R²_n = a + bφ
    print(f"n={n}: 4R² = {a:5d} + {b:5d}φ")

# n=0: 4R² =    10 +     9φ
# n=1: 4R² =    19 +    28φ   ← (10+9, 10+18) = (19, 28) ✓
# n=2: 4R² =    47 +    75φ   ← (19+28, 19+56) = (47, 75) ✓
# n=3: 4R² =   122 +   197φ
# n=4: 4R² =   319 +   516φ
# n=5: 4R² =   835 +  1351φ
# n=6: 4R² =  2186 +  3537φ
```

### 4.3.4　係数の増大とフィボナッチ数列の関係

$(a_n, b_n)$ の列をよく観察すると、$(10, 9), (19, 28), (47, 75), (122, 197), \ldots$ となっており、隣接する係数の比が $\varphi$ に近づいていく。これは $\varphi^2$ の乗算を繰り返すと固有ベクトル方向（$\varphi$ の方向）に収束するためであり、フィボナッチ数列の漸化式と同じ構造を持つ。

![Z[φ] 漸化式の可視化](/images/fig_ch4_zphi.png)
*図 4-2　$4R^2_n = a_n + b_n\varphi$ の整数成分（棒グラフ）と実数値 $R^2$（折れ線）の推移。$\varphi^2 \approx 2.618$ 倍で増大する。*

---

## 4.4　極角の全レベル不変性

### 4.4.1　不変性の事実

`cos_theta` と `sin_theta` はレベル $n$ によらず一定である。

```python
from b13phase import evaluate

for n in range(6):
    v = evaluate(k=1, j=0, n=n)
    print(f"n={n}: cos_theta={v.cos_theta}, sin_theta={v.sin_theta}")

# n=0: cos_theta=1395, sin_theta=2791
# n=1: cos_theta=1395, sin_theta=2791  ← 変化なし
# n=2: cos_theta=1395, sin_theta=2791
# ...（全レベルで同一）
```

### 4.4.2　なぜ不変なのか

球面座標 $(r, \theta, \phi_{\text{az}})$ において、極角 $\theta$ は $z/r$ の比で決まる。

$$
\cos\theta = \frac{z}{r}
$$

フラクタル展開では $r$ と $z$ が同じ因子 $\varphi^n$ で拡大する。すなわち

$$
r_n = \varphi^n \cdot r_0, \quad z_n = \varphi^n \cdot z_0
$$

であるから、比 $z_n / r_n = z_0 / r_0$ は $n$ によらず一定である。$\theta$ はレベルに依存しない。

これは幾何的に自然な結果である。フラクタル展開は「球を相似拡大する」操作であり、角度は保存される。

### 4.4.3　整数近似の誤差

`cos_theta` と `sin_theta` の値は $\cos\theta$ と $\sin\theta$ の整数近似であり、BASE でスケールされている。

$$
\text{cos\_theta} = \text{round}(\text{BASE} \times \cos\theta)
$$

近似誤差は $\leq 0.02\%$ であり、`verify_all()` で確認できる。

| バンド | $\cos\theta$ の真値 | `cos_theta` / BASE | 誤差 |
|--------|--------------------|--------------------|------|
| 北極   | $\varphi / \sqrt{3} \approx 0.9342$ | $2915/3120 \approx 0.9342$ | $< 0.01\%$ |
| 上帯   | $1/\sqrt{5} \approx 0.4472$ | $1395/3120 \approx 0.4471$ | $< 0.02\%$ |
| 下帯   | $-1/\sqrt{5}$ | $-1395/3120$ | $< 0.02\%$ |
| 南極   | $-\varphi/\sqrt{3}$ | $-2915/3120$ | $< 0.01\%$ |

---

## 4.5　フラクタル多レベルの全体像

3つのレベルを重ねて描くと、フラクタル構造が視覚的に確認できる。

```python
from b13phase import evaluate
from b13phase.evaluator import to_float_coords
import math

for n in range(3):
    coords = [to_float_coords(evaluate(k, j, n))
              for k in range(12) for j in range(5)]
    rs = [math.sqrt(x**2+y**2+z**2) for x,y,z in coords]
    print(f"n={n}: R = {rs[0]:.4f}  (全60頂点で同一)")

# n=0: R = 2.4780
# n=1: R = 4.0095   ← 2.4780 × φ ≈ 4.010 ✓
# n=2: R = 6.4875   ← 4.0095 × φ ≈ 6.488 ✓
```

![フラクタル多レベル頂点](/images/fig_ch4_multilevels.png)
*図 4-3　$n=0, 1, 2$ の60頂点を重ねた図。角度構造が完全に保存され、半径のみが $\varphi$ 倍で拡大する。*

---

## 4.6　まとめ

本章で示した内容を整理する。

1. **BASE = 3120** は $\mathbb{Z}_{60}$ の周期（÷60 = 52）と10段階のフラクタル構造（÷10 = 312）を同時に整数で表現するための設計値である。
2. **az_idx 公式** $= (m \times 52 + n \times 312) \bmod 3120$ は、頂点ラベルとレベルを方位角インデックスに変換する。レベルあたり 36° のシフトが生じる。
3. **$\mathbb{Z}[\varphi]$ 算術** により、半径の2乗 $R^2_n$ が $(a, b) \to (a+b, a+2b)$ という加算のみの更新式で厳密に計算できる。浮動小数点乗算は不要である。
4. **極角不変性**：フラクタル展開は球の相似拡大であるため、$\theta$ はすべてのレベルで一定である。`cos_theta` / `sin_theta` は定数テーブルから参照するだけでよい。

:::message
**設計の一貫性**　az_idx（整数テーブル参照）・$R^2$（$\mathbb{Z}[\varphi]$ 整数漸化式）・$\theta$（定数参照）の3つがすべて整数演算に帰着することで、`evaluate()` の全出力が整数となる。これが B13 の設計原理である。
:::
