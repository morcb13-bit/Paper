---
title: "第5章　大円閉合定理"
---

本章も初読スキップ可である。BASE = 3120 という値がなぜ「3120」でなければならないのか、その幾何学的な根拠を示す。

:::message
**初読スキップの目安**　BASE = 3120 を所与として使うだけなら本章は不要である。「なぜ3120なのか」という問いが気になるなら読み進めてほしい。
:::

---

## 5.1　大円とは何か

球面上の大円（great circle）とは、球の中心を含む平面と球面の交線であり、球面上の最大円である。切端20面体の表面には、対称性によって定まる複数の大円が存在する。

B13 の観点から重要な大円は、**頂点2個とエッジ中点8個の計10点**を通過するものである。この10点は球面上で正10角形（decagon）に近い配置をとる。

---

## 5.2　10点の座標

エッジ長 $a=2$ の切端20面体において、$x=0$ 断面上に現れる10点の座標は以下のとおりである（$\varphi = (1+\sqrt{5})/2$）。

| 種別 | 座標 $(0, y, z)$ | 半径 $r$ |
|------|-----------------|---------|
| 頂点 $v_0$ | $(0,\ 1,\ 3\varphi)$ | $\sqrt{1+9\varphi^2}$ |
| エッジ中点 $e_0$ | $(0,\ \varphi+2,\ 2\varphi)$ | $\sqrt{(\varphi+2)^2+(2\varphi)^2} = 3\varphi$ |
| エッジ中点 $e_1$ | $(0,\ 3\varphi,\ 0)$ | $3\varphi$ |

この断面を球面に投影すると、各点は大円上に等間隔ではないが規則的に並ぶ。

---

## 5.3　大円閉合定理

大円を1周するとき、弧区間は3種類に分類される。

- $\theta_A$：頂点 ↔ エッジ中点 の弧（8本）
- $\theta_B$：エッジ中点 ↔ エッジ中点 の弧（4本、隣接するもの）
- $\theta_{\text{edge}}$：頂点 ↔ 頂点 の弧（2本、実際のエッジ）

**大円閉合定理**：

$$
4\theta_A + 4\theta_B + 2\theta_{\text{edge}} = 360°
$$

各弧の中心角の数値は以下のとおりである。

| 弧 | 中心角 |
|----|--------|
| $\theta_A$ | $\approx 36.55°$ |
| $\theta_B$ | $\approx 41.81°$ |
| $\theta_{\text{edge}}$ | $\approx 23.28°$ |
| **合計**（$4+4+2$ の重みで） | **$= 360.000°$（厳密）** |

```python
from b13phase.great_circle import great_circle_closure_proof

proof = great_circle_closure_proof()
print(f"θ_A     = {proof['theta_A_deg']:.6f}°")
print(f"θ_B     = {proof['theta_B_deg']:.6f}°")
print(f"θ_edge  = {proof['theta_edge_deg']:.6f}°")
print(f"4θ_A + 4θ_B + 2θ_edge = {proof['full_sum_deg']:.6f}°")
print(f"閉合証明: {proof['closure_proved']}")

# θ_A     = 36.548962°
# θ_B     = 41.810315°
# θ_edge  = 23.281446°
# 4θ_A + 4θ_B + 2θ_edge = 360.000000°
# 閉合証明: True
```

![大円閉合の断面図](/images/fig_ch5_great_circle.png)
*図 5-1　$x=0$ 断面上の10点と3種類の弧区間。青：$\theta_A$、緑：$\theta_B$、橙：$\theta_{\text{edge}}$。*

---

## 5.4　φ 算術による厳密計算

中心角の余弦値は $\mathbb{Z}[\varphi]$ の有理式として厳密に表現できる。

### 5.4.1　cos(θ_B) と sin(θ_B)

エッジ中点 $e_0 = (\varphi+2,\ 2\varphi)$ と $e_1 = (3\varphi,\ 0)$ の半径を計算する。

$$
|e_0|^2 = (\varphi+2)^2 + (2\varphi)^2 = \varphi^2+4\varphi+4 + 4\varphi^2 = 5\varphi^2 + 4\varphi + 4
$$

$\varphi^2 = \varphi+1$ を代入すると：

$$
= 5(\varphi+1) + 4\varphi + 4 = 9\varphi + 9 = 9(\varphi+1) = 9\varphi^2
$$

したがって $|e_0| = 3\varphi$。$|e_1| = 3\varphi$ も同様に確認できる。

内積を計算すると：

$$
e_0 \cdot e_1 = (\varphi+2)(3\varphi) + (2\varphi)(0) = 3\varphi(\varphi+2) = 3\varphi^2 + 6\varphi = 3(\varphi+1)+6\varphi = 9\varphi+3
$$

$$
\cos(\theta_B) = \frac{e_0 \cdot e_1}{|e_0||e_1|} = \frac{9\varphi+3}{9\varphi^2} = \frac{9\varphi+3}{9(\varphi+1)} = \frac{\varphi+2}{3\varphi}
$$

外積の大きさを計算すると：

$$
|e_0 \times e_1| = |(\varphi+2)\cdot 0 - 2\varphi \cdot 3\varphi| = 6\varphi^2
$$

$$
\sin(\theta_B) = \frac{6\varphi^2}{9\varphi^2} = \frac{2}{3}
$$

**$\sin(\theta_B) = 2/3$ は有理数である。** $\varphi$ が完全に約分される。

```python
from b13phase.great_circle import great_circle_closure_sin_b_origin

origin = great_circle_closure_sin_b_origin()
print(f"|r|² = {origin['r_mid_sq']:.6f}")
print(f"9φ²  = {origin['r_mid_phi2']:.6f}")  # 一致を確認
print(f"sin(θ_B) = {origin['sin_B']}  ({origin['sin_B_rational']})")
# sin(θ_B) = 0.6666666666666666  (True)
```

![sin(θ_B)=2/3 の幾何的説明](/images/fig_ch5_sinB.png)
*図 5-2　2つのエッジ中点ベクトルの外積比から $\sin(\theta_B) = 2/3$ が導かれる様子。$\varphi$ が分子・分母で完全に消える。*

### 5.4.2　閉合の代数的証明

定理 $4\theta_A + 4\theta_B + 2\theta_{\text{edge}} = 360°$ は半角に換算すると

$$
2\theta_A + 2\theta_B + \theta_{\text{edge}} = 180°
$$

すなわち $2\theta_A = 180° - (2\theta_B + \theta_{\text{edge}})$ であるから、$\cos(2\theta_A) = -\cos(2\theta_B + \theta_{\text{edge}})$ を示せばよい。

$\sin(\theta_B) = 2/3$ から $\cos(2\theta_B) = 1 - 2\sin^2(\theta_B) = 1 - 8/9 = 1/9$。

加法定理を用いて $\cos(2\theta_B + \theta_{\text{edge}})$ を展開し、整理すると：

$$
\cos(2\theta_B + \theta_{\text{edge}}) = -\frac{70\varphi+55}{252\varphi+171}
$$

一方 $\cos(2\theta_A)$ を $\varphi$ 算術で計算すると同じ値が得られる。これが閉合の代数的証明である。

```python
proof = great_circle_closure_proof()
print(f"cos(2θ_A)         = {proof['cos_2A_phi']:.10f}")
print(f"cos(2θ_B+θ_edge)  = {proof['cos_2B_edge']:.10f}")
print(f"-cos(2θ_B+θ_edge) = {-proof['cos_2B_edge']:.10f}")
print(f"一致: {proof['closure_proved']}")
# 一致: True
```

---

## 5.5　分母 109 と F(7)=13 の関係

$\varphi$ 算術の計算に現れる分母 109 は以下のように分解される。

$$
109 = 8 \times F(7) + F(5) = 8 \times 13 + 5
$$

ここで $F(7) = 13$ はフィボナッチ数であり、第3章で登場した $\mathbb{Z}_{13}^*$ の位数と一致する。分母に 13 が現れることは偶然ではなく、切端20面体の代数的構造が一貫していることを示している。

また $R^2 = (9\varphi+10)/4$ の分子 $9\varphi+10$ を有理化すると：

$$
\frac{1}{9\varphi+10} = \frac{19-9\varphi}{(9\varphi+10)(19-9\varphi)} = \frac{19-9\varphi}{109}
$$

$(9\varphi+10)(19-9\varphi) = 171+90\varphi - 81\varphi^2 - 90\varphi = 171 - 81(\varphi+1) = 90 - 81 = 109$ が成立する（$\varphi^2 = \varphi+1$ を使用）。

---

## 5.6　BASE = 3120 との接続

大円閉合定理が BASE = 3120 を正当化する論理は以下のとおりである。

切端20面体は10回対称性（10個の大円の交点）を持つ。1周 360° を整数で分割するとき、

- **10 分割**（36°/ステップ）が基本単位 → $\text{BASE}/10 = 312$
- **60 分割**（6°/ステップ）が頂点識別単位 → $\text{BASE}/60 = 52$

この2条件を同時に整数で満たす最小値に因子 13（$\mathbb{Z}_{13}^*$ の位数）を掛けると $3120 = \text{lcm}(60, 10) \times 13 \times 4$ となる。

```python
import math
lcm_60_10 = math.lcm(60, 10)  # = 60
print(f"lcm(60,10) = {lcm_60_10}")
print(f"lcm(60,10) × 13 × 4 = {lcm_60_10 * 13 * 4}")  # 3120
```

---

## 5.7　まとめ

本章で示した内容を整理する。

1. 切端20面体の大円は**10点**（頂点2個・エッジ中点8個）を通過し、3種類の弧区間 $\theta_A, \theta_B, \theta_{\text{edge}}$ から構成される。
2. **大円閉合定理** $4\theta_A + 4\theta_B + 2\theta_{\text{edge}} = 360°$ が厳密に成立する。
3. $\sin(\theta_B) = 2/3$ は**有理数**である。エッジ中点ベクトルの比で $\varphi$ が消約するためである。
4. 分母 109 = $8 \times F(7) + F(5)$ を通じて、$\mathbb{Z}_{13}^*$（第3章）・フィボナッチ数列・大円の幾何が同一の代数構造でつながっている。
5. BASE = 3120 は10回対称性と60頂点構造を同時に整数で表現するための最小設計値である。
