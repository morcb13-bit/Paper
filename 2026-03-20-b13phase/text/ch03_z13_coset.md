---
title: "第3章　Z₁₃* とコセット"
---

本章は初読スキップ可である。B13 を実用的に使うだけであれば、第2章と第4章の知識で十分である。本章では、頂点ラベル $m$ の分類を支配する代数構造 ── 乗法群 $\mathbb{Z}_{13}^*$ とそのコセット分解 ── を解説する。

:::message
**初読スキップの目安**　$m \bmod 3$ による3クラス分類（H / 2H / 4H）が「そういうものだ」として受け入れられるなら、本章は後回しにして第4章へ進んでよい。
:::

---

## 3.1　Z₁₃* とは何か

$\mathbb{Z}_{13}^*$ は、13 と互いに素な整数からなる乗法群である。13 は素数であるから $\mathbb{Z}_{13}^* = \{1, 2, 3, \ldots, 12\}$ であり、位数は 12 である。

B13 がこの群に注目する理由は、**2 が $\mathbb{Z}_{13}^*$ の原始根**であることにある。すなわち、2 の冪乗がこの群の全要素を生成する。

$$
2^0, 2^1, 2^2, \ldots, 2^{11} \pmod{13} = 1, 2, 4, 8, 3, 6, 12, 11, 9, 5, 10, 7
$$

この列が、第1章で登場した `k` の意味を与える。**`k` は $2^k \bmod 13$ という指数**である。

```python
# Z13* の生成
Z13_STAR = [(2**k) % 13 for k in range(12)]
print(Z13_STAR)
# [1, 2, 4, 8, 3, 6, 12, 11, 9, 5, 10, 7]
```

![Z13* の円環図](/images/fig_ch3_z13star.png)
*図 3-1　$\mathbb{Z}_{13}^*$ の元を円環に配置し、×2 の写像を矢印で示したもの。色はコセット（H / 2H / 4H）を示す。*

---

## 3.2　コセット分解：H / 2H / 4H

$\mathbb{Z}_{13}^*$ は位数 4 の部分群 $H$ を持ち、3 つのコセットに分解される。

$$
\mathbb{Z}_{13}^* = H \cup 2H \cup 4H
$$

各コセットの具体的な要素は以下のとおりである。

| コセット | 要素 | 特徴 |
|---------|------|------|
| $H$     | $\{1, 5, 8, 12\}$ | 部分群（閉じている） |
| $2H$    | $\{2, 3, 10, 11\}$ | $H$ の 2 倍 |
| $4H$    | $\{4, 6, 7, 9\}$  | $H$ の 4 倍（Fib/Pell の禁止域） |

部分群 $H$ は「立方剰余の部分群」であり、$\{1^2, 5^2, 8^2, 12^2\} \bmod 13 = \{1, 12, 12, 1\}$ のような閉包性を持つ。

```python
from b13phase.constants import H_SET, H2_SET, H4_SET

print(f"H  = {sorted(H_SET)}")   # [1, 5, 8, 12]
print(f"2H = {sorted(H2_SET)}")  # [2, 3, 10, 11]
print(f"4H = {sorted(H4_SET)}")  # [4, 6, 7, 9]
```

### 3.2.1　×2 による周期パターン

$\mathbb{Z}_{13}^*$ の元を $k = 0, 1, \ldots, 11$ の順（$2^k \bmod 13$）に並べたとき、コセットのパターンは周期 3 で繰り返す。

$$
H,\ 2H,\ 4H,\ H,\ 2H,\ 4H,\ H,\ 2H,\ 4H,\ H,\ 2H,\ 4H
$$

```python
from b13phase.z13_structure import coset_pattern
print(coset_pattern())
# ['H', '2H', '4H', 'H', '2H', '4H', 'H', '2H', '4H', 'H', '2H', '4H']
```

この周期 3 のパターンが、$m \bmod 3$ による3分類と完全に対応する（3.3 節で詳述）。

---

## 3.3　m mod 3 と Z₁₃* コセットの対応

第2章で述べた $m \bmod 3$ によるコセット分類は、$\mathbb{Z}_{13}^*$ の代数構造と一致する。具体的には次の対応が成立する。

| $m \bmod 3$ | $\mathbb{Z}_{13}^*$ コセット | 対応する $k$ の値 |
|------------|--------------------------|----------------|
| 0          | $H$                      | 0, 3, 6, 9     |
| 1          | $4H$                     | 2, 5, 8, 11    |
| 2          | $2H$                     | 1, 4, 7, 10    |

この対応は偶然ではない。$5k \bmod 3$ は $2k \bmod 3$ と一致し（$5 \equiv 2 \pmod{3}$）、$k \bmod 3$ が 0/1/2 のとき $5k \bmod 3$ は 0/2/1 となる。これが上表の順序を生む。

```python
from b13phase.constants import Z13_STAR
from b13phase.z13_structure import coset_of

for k in range(12):
    g = Z13_STAR[k]
    m_mod3 = (5*k) % 3
    print(f"k={k:2d}: 2^k={g:2d}, coset={coset_of(g):3s}, m mod 3={m_mod3}")

# k= 0: 2^k= 1, coset=H  , m mod 3=0
# k= 1: 2^k= 2, coset=2H , m mod 3=2
# k= 2: 2^k= 4, coset=4H , m mod 3=1
# k= 3: 2^k= 8, coset=H  , m mod 3=0
# ...（周期 3 で繰り返す）
```

![m mod 3 と Z13* コセットの対応表](/images/fig_ch3_coset_map.png)
*図 3-2　各 `k` に対する $m \bmod 3$（列）と $\mathbb{Z}_{13}^*$ コセット（色）の完全対応。*

:::message
**引継書 v60 より確認済み**　$m \bmod 3$ と $\mathbb{Z}_{13}^*$ コセットの完全一致は `verify_all()` のテストスイートで検証されている。
:::

---

## 3.4　Fibonacci / Pell の 4H 回避

$\mathbb{Z}_{13}^*$ の文脈で注目すべき定理がある。**フィボナッチ数列とペル数列は $\bmod 13$ において 4H コセット（$\{4, 6, 7, 9\}$）に決して落ちない**というものである。

### 3.4.1　フィボナッチ数列の 4H 回避

フィボナッチ数列 $F(n)$ を 13 で割った余りの列は周期 7 で繰り返す（Pisano 周期）。

$$
F(n) \bmod 13 : 0, 1, 1, 2, 3, 5, 8, 0, 8, 8, 3, 11, 1, 12, 0, 12, 12, \ldots
$$

非零の項はすべて $H$ または $2H$ に属し、$4H = \{4, 6, 7, 9\}$ には一度も現れない。

### 3.4.2　ペル数列の 4H 回避

ペル数列 $P(n)$（$P(0)=0,\ P(1)=1,\ P(n) = 2P(n-1)+P(n-2)$）も同様に、$\bmod 13$ において $4H$ を回避する。

```python
from b13phase.z13_structure import fib_mod13, pell_mod13, is_forbidden, verify_4h_avoidance

# 個別確認
print("Fib mod 13 の先頭 14 項:")
for n in range(14):
    f = fib_mod13(n)
    print(f"  F({n:2d}) = {f:2d}  {'[4H!]' if is_forbidden(f) else ''}")

# 一括検証
ok = verify_4h_avoidance(n_max=200)
print(f"\nn=0..200 で 4H 回避: {ok}")  # True
```

![Fib/Pell の mod 13 コセット軌跡](/images/fig_ch3_fib_pell.png)
*図 3-3　$F(n) \bmod 13$（左）と $P(n) \bmod 13$（右）のコセット分布。オレンジ（4H）には一度も現れない。*

### 3.4.3　なぜ 4H を回避するのか

直感的な理由を示す。フィボナッチ数列の生成行列

$$
M = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}
$$

に対して、$M^7 \equiv 8I \pmod{13}$ が成立する（$F(7) = 13 \equiv 0$ が根拠）。$8 \in H$ であるから、$M^7$ の固有値的な振る舞いが $H$ に収まる。これが数列の値を $H \cup 2H$ に制限する。

:::message
**代数的証明の状態**　引継書 v60 時点では、この定理の完全な代数的証明は「継続中」の状態である。数値的には $n = 0 \sim 200$ で確認済み。
:::

---

## 3.5　2軸の独立性：band と coset

B13 の設計において重要な事実は、**band（極角の分類）と coset（$m \bmod 3$）が独立した2軸をなす**ことである。

- **band**（`band(k)` = 0/1/2/3）は `k` のみで決まる
- **coset**（$m \bmod 3$）は `k` のみで決まる（`j` に依存しない）

しかし、band と coset は異なる写像によって定まるため、互いに独立である。

```python
from b13phase.evaluator import band as get_band
from b13phase.constants import Z13_STAR
from b13phase.z13_structure import coset_of

print(f"{'k':>3}  {'band':>5}  {'coset':>5}  {'m mod 3':>7}")
print("-" * 28)
for k in range(12):
    b = get_band(k)
    c = coset_of(Z13_STAR[k])
    m3 = (5*k) % 3
    print(f"{k:>3}  {b:>5}  {c:>5}  {m3:>7}")
```

出力から、band=1（上帯）の中に coset が H・2H・4H の全てが混在することが確認できる。band と coset の組み合わせは独立に選べる。

この独立性が `evaluator.py` の実装に直接反映されている。`theta_ints(k)` は band のみを参照し、`az_idx(k, j, n)` はコセット（$m \bmod 3$）を経由した `m` を使う。2つの計算が干渉しない。

---

## 3.6　まとめ

本章で示した内容を整理する。

1. $\mathbb{Z}_{13}^*$ は位数 12 の乗法群であり、2 を原始根とする。`k` は $2^k \bmod 13$ の指数に対応する。
2. $\mathbb{Z}_{13}^*$ は部分群 $H = \{1, 5, 8, 12\}$ によって $H \cup 2H \cup 4H$ に分解される。この周期 3 のパターンが $m \bmod 3$ に対応する。
3. フィボナッチ数列とペル数列は $\bmod 13$ において $4H$ コセットに決して現れない（$n = 0 \sim 200$ で確認済み）。
4. band（極角分類）と coset（$m \bmod 3$）は独立した2軸であり、`evaluator.py` の設計の基礎をなす。
