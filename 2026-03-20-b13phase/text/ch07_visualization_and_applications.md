---
title: "第7章　可視化と応用"
---

本章では `evaluate()` を実際に使う場面を示す。「頂点データをどう取り出し、何ができるか」という視点で、5つのレシピを順に解説する。

---

## 7.1　全60頂点を座標として取り出す

最もシンプルな使い方である。全60頂点を浮動小数点座標のリストに変換する。

```python
from b13phase import evaluate
from b13phase.evaluator import to_float_coords

# 全60頂点の (x, y, z) を生成
points = []
for k in range(12):
    for j in range(5):
        v = evaluate(k, j, n=0)
        x, y, z = to_float_coords(v)
        points.append((x, y, z))

print(len(points))  # 60
```

`to_float_coords()` はデバッグ・可視化専用である。B13 の計算はすべて `B13Vertex` の整数フィールドで行う。

### matplotlib による3D散布図

```python
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

xs = [p[0] for p in points]
ys = [p[1] for p in points]
zs = [p[2] for p in points]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs, ys, zs, s=20)
ax.set_title("B13 — 60 vertices (n=0)")
plt.show()
```

### バンドで色分け

```python
band_colors = {0: 'red', 1: 'orange', 2: 'steelblue', 3: 'navy'}

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for k in range(12):
    for j in range(5):
        v = evaluate(k, j, n=0)
        x, y, z = to_float_coords(v)
        # band(k): 0=北極, 1=上帯, 2=下帯, 3=南極
        b = 0 if k == 0 else (3 if k == 11 else (1 if k <= 5 else 2))
        ax.scatter(x, y, z, color=band_colors[b], s=30)

plt.show()
```

北極（赤）・上帯（橙）・下帯（青）・南極（紺）の4バンドが球面上に分布する様子が視覚的に確認できる。

---

## 7.2　フラクタル展開の可視化

レベル `n` を増やすごとに頂点が球の外側へ螺旋状に展開する。

```python
import math

phi = (1 + math.sqrt(5)) / 2  # 黄金比

print(f"{'n':>3}  {'az(°)':>7}  {'R':>8}  {'4R²':}")
for n in range(6):
    v = evaluate(k=0, j=0, n=n)
    az_deg = v.az_idx / 3120 * 360
    R = math.sqrt((v.r_sq_rat + v.r_sq_phi * phi) / 4)
    print(f"{n:>3}  {az_deg:>7.1f}  {R:>8.4f}  {v.r_sq_rat}+{v.r_sq_phi}φ")
```

出力例：

```
  n    az(°)         R  4R²
  0      0.0    2.4782  10+9φ
  1     36.0    4.0098  19+28φ
  2     72.0    6.4880  47+75φ
  3    108.0   10.4978  122+197φ
  4    144.0   16.9860  319+516φ
  5    180.0   27.4838  835+1351φ
```

**2つの不変量**

1. `cos_theta` はレベル `n` によらず変化しない（極角が一定）。
2. `az_idx` は毎レベル 312（= 36°）ずつ増加する。

**R は φ 倍ずつ拡大する**

`4R²` の更新式は `(a, b) → (a+b, a+2b)` であり、加算だけで計算できる。

```python
from b13phase.evaluator import r_sq_level

for n in range(6):
    a, b = r_sq_level(n)
    print(f"n={n}: 4R² = {a} + {b}φ")
```

### 螺旋プロット

同一 `(k, j)` の頂点を複数レベルにわたって追跡し、螺旋を描く。

```python
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for k in range(0, 12, 3):          # 4系統だけ描画
    xs, ys, zs = [], [], []
    for n in range(8):
        v = evaluate(k, 0, n=n)
        x, y, z = to_float_coords(v)
        xs.append(x); ys.append(y); zs.append(z)
    ax.plot(xs, ys, zs, marker='o', markersize=3)

ax.set_title("B13 フラクタル展開（n=0〜7）")
plt.show()
```

---

## 7.3　対蹠点の整数的確認

B13 では「対蹠点」を浮動小数点比較なしで求められる。

対蹠の規則：`k` を `(k+6) mod 12` に置き換えると対蹠点になる。

```python
def antipodal(k, j, n=0):
    """(k, j, n) の対蹠頂点を返す。"""
    return evaluate((k + 6) % 12, j, n=n)
```

整数フィールドでの確認：

```python
v = evaluate(k=2, j=1, n=0)
v_ant = antipodal(k=2, j=1, n=0)

# az_idx の差が常に BASE/2 = 1560（全ペアで成立）
assert (v_ant.az_idx - v.az_idx) % 3120 == 1560

# r_sq は同じ（同レベルなら）
assert v_ant.r_sq_rat == v.r_sq_rat
assert v_ant.r_sq_phi == v.r_sq_phi

# cos_theta の符号反転は上帯(k=1..5) ↔ 下帯(k=6..10)ペアのみ
# 北極(k=0)の対蹠は k=6（下帯）なのでcos_thetaの絶対値が異なる
```

`az_idx` の差が 1560 になることは全ペアで成立し、整数の等号比較で確認できる。

---

## 7.4　Fib/Pell の 4H 回避を整数で確認する

第3章で述べたとおり、Fibonacci 数列と Pell 数列の項を mod 13 で見ると、4H 帯（`k%3 == 2`）には絶対に入らない。この事実を evaluator と整数比較だけで確認する。

```python
from b13phase.constants import Z13_STAR
from b13phase.z13_structure import coset_of

def fib_mod13(limit=200):
    """F(n) mod 13 の非ゼロ項を (n, val, k) のリストで返す。"""
    a, b = 0, 1
    result = []
    for n in range(limit):
        if a != 0 and a in Z13_STAR:
            k = Z13_STAR.index(a)
            result.append((n, a, k))
        a, b = b, (a + b) % 13
    return result

orbit = fib_mod13(200)

# 4H帯に入るものがないか確認
h4_hits = [(n, val, k) for n, val, k in orbit if k % 3 == 2]
print(f"4H帯への到達: {len(h4_hits)} 件")   # → 0 件
```

evaluator と組み合わせて、軌道上の各頂点のコセットを表示する：

```python
print(f"{'n':>4}  {'F mod 13':>8}  {'k':>3}  coset  az_idx")
for n, val, k in orbit[:10]:
    v = evaluate(k, 0, n=0)
    coset = coset_of(val)
    print(f"{n:>4}  {val:>8}  {k:>3}  {coset:5s}  {v.az_idx:6d}")
```

出力例：

```
   n  F mod 13    k  coset  az_idx
   1         1    0  H         0
   2         1    0  H         0
   3         2    1  2H      260
   4         3    4  2H     1040
   5         5    9  H      2340
   6         8    3  H       780
   8         8    3  H       780
   ...
```

coset 列が `H` か `2H` のみで `4H` が現れない。

---

## 7.5　整数フィールドによる頂点の同一性判定

B13 の利点のひとつは、2頂点が「同じ位置にある」かを整数比較だけで判定できることである。

```python
def same_vertex(v1, v2):
    """整数フィールドだけで頂点の位置的同一性を判定する。"""
    return (v1.az_idx    == v2.az_idx and
            v1.cos_theta == v2.cos_theta and
            v1.r_sq_rat  == v2.r_sq_rat and
            v1.r_sq_phi  == v2.r_sq_phi)
```

使用例：

```python
va = evaluate(0, 0, n=0)
vb = evaluate(0, 0, n=0)   # 同じ頂点
vc = evaluate(1, 0, n=0)   # 別の頂点

print(same_vertex(va, vb))  # True
print(same_vertex(va, vc))  # False
```

`to_float_coords()` を使った浮動小数点比較では、`== 0.000` のような判定が桁落ちで失敗することがある。整数比較は誤りがない。

:::message
**注意**　`same_vertex()` は「同じ `(k, j, n)` から生成された頂点か」を判定するのに十分である。異なる `(k, j, n)` が偶然に同じ位置に来ることはない（全単射性）。
:::

---

## 7.6　頂点ラベル m による直接アクセス

`m = (5k + 12j) mod 60` という頂点ラベルから、対応する `(k, j)` を求めて `evaluate()` を呼ぶ。

```python
def evaluate_by_m(m, n=0):
    """頂点ラベル m ∈ Z₆₀ から B13Vertex を返す。"""
    for k in range(12):
        for j in range(5):
            if (5*k + 12*j) % 60 == m:
                return evaluate(k, j, n)
    raise ValueError(f"m={m} に対応する (k, j) が見つかりません")
```

あるいは中国剰余定理（CRT）で直接求める方が速い。

```python
def m_to_kj(m):
    """m → (k, j) を CRT で直接計算する。"""
    # 5k ≡ m (mod 12)  →  k ≡ 5m (mod 12)   [5^{-1} ≡ 5 mod 12]
    # 12j ≡ m (mod 5)  →  j ≡ 3m (mod 5)    [12 ≡ 2 mod 5, 2^{-1} ≡ 3 mod 5]
    k = (5 * m) % 12
    j = (3 * m) % 5
    assert (5*k + 12*j) % 60 == m   # 検証
    return k, j
```

使用例：

```python
k, j = m_to_kj(15)
v = evaluate(k, j, n=0)
print(f"m=15 → k={k}, j={j}")
# m=15 → k=3, j=0
```

---

## 7.7　まとめ：B13 の使い方の原則

本章で示した5つのレシピを整理する。

| やりたいこと | 使う関数・手法 |
|-------------|--------------|
| 全60頂点の座標リスト | `evaluate(k, j, n)` のループ + `to_float_coords()` |
| フラクタル展開の追跡 | `n` を変えて `evaluate()` + `r_sq_level()` |
| 対蹠点の取得 | `k` → `(k+6)%12` に置き換えて `evaluate()` |
| 4H回避の確認 | `k%3 == 2` で判定 |
| 頂点の同一性判定 | `az_idx`, `cos_theta`, `r_sq_*` の整数比較 |
| m ラベルでのアクセス | `m_to_kj(m)` → `evaluate(k, j, n)` |

**整数フィールドで完結させる**のが B13 の使い方の基本である。`to_float_coords()` は可視化の最終ステップとして使い、比較・計算・判定には常に整数フィールドを使う。

```python
# 典型的な処理フロー
from b13phase import evaluate, B13Vertex
from b13phase.evaluator import to_float_coords

v = evaluate(k, j, n)       # ① 整数レコードを取得
# ↓ 整数のまま処理（比較・フィルタ・集計）
if v.cos_theta > 0:          # ② 整数比較
    pass
# ↓ 最後に変換（可視化・外部ライブラリへの受け渡し）
x, y, z = to_float_coords(v) # ③ 浮動小数点に変換
```
