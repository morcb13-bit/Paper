---
title: "第2章　頂点ラベル m ∈ Z₆₀"
---

本章では、B13 の頂点識別に用いる整数ラベル $m \in \mathbb{Z}_{60}$ の性質を解説する。$m = (5k + 12j) \bmod 60$ という式が全単射となる理由、バンド（帯）による極角分類との関係、そして対蹠点の代数的表現を順に示す。

---

## 2.1　頂点ラベル m の定義

B13 は $k \in \mathbb{Z}_{12}$ と $j \in \mathbb{Z}_5$ の2つの整数から、次の式によって頂点ラベル $m$ を計算する。

$$
m = (5k + 12j) \bmod 60
$$

$5k$ は $k$ に応じて $0, 5, 10, \ldots, 55$ のいずれかを与え、$12j$ は $j$ に応じて $0, 12, 24, 36, 48$ のいずれかを加算する。結果を $\bmod 60$ で $\mathbb{Z}_{60}$ に収める。

```python
# m = (5k+12j) mod 60 の全写像
m_table = {(k, j): (5*k + 12*j) % 60
           for k in range(12)
           for j in range(5)}

# 値が 0〜59 のすべてを網羅しているか確認
assert set(m_table.values()) == set(range(60))  # 全単射の確認
```

下図は全 `(k, j)` ペアに対する $m$ の値を示したものである。各値がちょうど1回ずつ現れることが確認できる。色はコセット $m \bmod 3$ を示す（詳細は 2.5 節）。

![m = (5k+12j) mod 60 の写像テーブル](/images/fig_ch2_bijection.png)
*図 2-1　$m = (5k + 12j) \bmod 60$ の全60値。各値は1回ずつ現れる（全単射）。*

---

## 2.2　全単射の証明

写像 $(k, j) \mapsto m$ が $\mathbb{Z}_{12} \times \mathbb{Z}_5 \to \mathbb{Z}_{60}$ の全単射であることを示す。定義域と値域の要素数はともに60であるから、単射であることを示せば十分である。

### 2.2.1　gcd(5, 12) = 1 の役割

5 と 12 は互いに素（$\gcd(5, 12) = 1$）である。これにより、$5k + 12j$ の値は $\bmod 60$ において $0$ から $59$ のすべての整数を生成する。

直感的には、$k$ を1増やすと $m$ は $+5$、$j$ を1増やすと $m$ は $+12$ だけ変化する。この2つのステップが互いに素であるため、格子点がちょうど60個の異なる剰余クラスを過不足なく埋め尽くす。

:::message
**数学的背景**　$\mathbb{Z}_n$ 上で $ax + by \equiv r \pmod{n}$ が全ての $r$ に対して解を持つ必要十分条件は $\gcd(a, b, n) \mid r$ である。ここでは $a=5,\; b=12,\; n=60$ かつ $\gcd(5, 12, 60) = 1$ であるため、$r = 0 \sim 59$ の全値が実現する。
:::

### 2.2.2　コードによる確認

```python
def check_bijection():
    seen = set()
    for k in range(12):
        for j in range(5):
            m = (5*k + 12*j) % 60
            assert m not in seen, f"重複: m={m} at k={k}, j={j}"
            seen.add(m)
    assert seen == set(range(60))
    print("全単射を確認: 60個の異なる値")

check_bijection()  # 全単射を確認: 60個の異なる値
```

---

## 2.3　バンド分類と極角の4値

全60頂点は `k` の値のみによって4つの帯（バンド）に分類される。`j` の値はバンドに影響しない。バンドは頂点の z 座標（極角 $\theta$）のクラスを決定する。

### 2.3.1　band(k) の定義

```python
def band(k):
    if k == 0:   return 0  # 北極
    if k <= 5:   return 1  # 上帯
    if k <= 10:  return 2  # 下帯
    return 3               # 南極
```

同じバンド内では $\cos\theta$ と $\sin\theta$ が共通の整数値を持つ（第1章の表を参照）。この整数値は4種類しか存在しない。

### 2.3.2　j の変化が z 座標に影響しない理由

$m = (5k + 12j) \bmod 60$ において、$j$ を変えると `az_idx`（方位角インデックス）は変化するが、`band(k)` は `k` のみで決まるため極角は変化しない。

```python
from b13phase import evaluate

k = 3
for j in range(5):
    v = evaluate(k, j, n=0)
    print(f"j={j}: az_idx={v.az_idx:4d}, cos_theta={v.cos_theta}")

# j=0: az_idx= 780, cos_theta=1395
# j=1: az_idx=1404, cos_theta=1395   <- az_idx は変化
# j=2: az_idx=  28, cos_theta=1395   <- cos_theta は不変
# j=3: az_idx= 652, cos_theta=1395
# j=4: az_idx=1276, cos_theta=1395
```

---

## 2.4　対蹠点：m と m + 30

切端20面体は中心対称性を持つ。任意の頂点 $m$ に対して、原点を挟んで正反対に位置する頂点（対蹠点）が存在し、そのラベルは $(m + 30) \bmod 60$ である。

$$
\text{antipodal}(m) = (m + 30) \bmod 60
$$

これは $\mathbb{Z}_{60}$ の加法的な対称性として表現される。$30 = 60/2$ であり、$\mathbb{Z}_{60}$ 上で「半周」に相当する。

### 2.4.1　コードによる確認

```python
from b13phase import evaluate
from b13phase.evaluator import to_float_coords

def find_vertex_by_m(m_target):
    for k in range(12):
        for j in range(5):
            if (5*k + 12*j) % 60 == m_target:
                return evaluate(k, j, n=0)

# m=0 とその対蹠点 m=30 を検証
v0  = find_vertex_by_m(0)
v30 = find_vertex_by_m(30)
x0, y0, z0   = to_float_coords(v0)
x30, y30, z30 = to_float_coords(v30)

dot = x0*x30 + y0*y30 + z0*z30
print(f"dot product = {dot:.4f}")
# dot product = -6.1414  （= -R²、対蹠点の内積は負の R²）
```

内積が $-R^2$ となる（負値で絶対値が $R^2$ に等しい）ことが、2点が正反対の位置にあることの証拠である。

![対蹠点ペアの3D図](/images/fig_ch2_antipodal.png)
*図 2-2　対蹠点ペア（同色の2点が対蹠関係）。各ペアを結ぶ線分は原点を通過する。*

### 2.4.2　対蹠点のバンド関係

| 頂点のバンド     | 対蹠点のバンド   | 例（m → m+30） |
|-----------------|-----------------|---------------|
| 北極 (k=0)      | 下帯 (k=6〜10)  | 0 → 30        |
| 上帯 (k=1〜5)   | 下帯 (k=6〜10)  | 5 → 35        |
| 下帯 (k=6〜10)  | 上帯・北極      | 30 → 0        |
| 南極 (k=11)     | 上帯 (k=1〜5)   | 55 → 25       |

---

## 2.5　コセット軸：m mod 3

頂点ラベル $m$ を 3 で割った余り（$m \bmod 3$）は、第3章で扱う $\mathbb{Z}_{13}^*$ のコセット分類に対応する。本章では実用上の意味にとどめ、代数的詳細は第3章に委ねる。

$m \bmod 3$ の値は 0・1・2 の3クラスに頂点を等分する（各20頂点）。この分類を「H コセット」「2H コセット」「4H コセット」と呼ぶ。

```python
from b13phase.evaluator import band

for k in range(12):
    b = band(k)
    mod3 = (5*k) % 3  # 12j mod 3 = 0 なので k のみで決まる
    print(f"k={k:2d} (band={b}): m mod 3 = {mod3}  （全 j 共通）")

# k= 0 (band=0): m mod 3 = 0
# k= 1 (band=1): m mod 3 = 2
# k= 2 (band=1): m mod 3 = 1
# k= 3 (band=1): m mod 3 = 0
# ...
```

![バンドとコセットのマトリクス](/images/fig_ch2_band_coset.png)
*図 2-3　バンド（行）とコセット（列）の頂点数分布。同一バンド内では $m \bmod 3$ が一定である。*

重要な観察を2点示す。

1. **バンド内で $m \bmod 3$ は一定**。$12j \bmod 3 = 0$（12 は 3 の倍数）であるため、$m \bmod 3 = 5k \bmod 3$ となり、$j$ に依存しない。
2. **バンドとコセットは独立した2軸**。バンドは $k$ のみで、コセットは $5k \bmod 3$ で決まる。この独立性が evaluator.py の設計の基礎となっている（詳細は第3章）。

---

## 2.6　まとめ

本章で示した内容を整理する。

1. $m = (5k + 12j) \bmod 60$ は $\mathbb{Z}_{12} \times \mathbb{Z}_5 \to \mathbb{Z}_{60}$ の全単射である。$\gcd(5, 12) = 1$ がこれを保証する。
2. バンド $\text{band}(k)$ は $k$ のみで決まり、$j$ の変化に対して不変である。同一バンド内では $\cos\theta$ と $\sin\theta$ が共通の整数値をとる。
3. 対蹠点は $m \mapsto (m + 30) \bmod 60$ で与えられる。$\mathbb{Z}_{60}$ の加法的対称性として表現される。
4. コセット分類 $m \bmod 3$ はバンドとは独立した軸をなす。この代数的構造は第3章（$\mathbb{Z}_{13}^*$）で詳述する。

:::message
**実用上の要点**　`k` を決めればバンド（極角）が定まり、`j` を決めれば方位角が定まる。この2軸の独立性が `evaluator.py` の設計の基礎となっている。
:::
