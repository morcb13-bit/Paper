# 整数軸上の閾値干渉モデル
## √2・√5・フィボナッチ閾値の交差で発火を見る

## 概要

通常の閾値モデルは、

- 残渣 \(R_n\)
- 固定閾値 \(T\)

に対して、

\[
R_n > T
\]

なら発火する、という単純な形をとる。

しかし、今回見たいのはそういう単一閾値モデルではない。

ここでは、整数軸 \(n\) 上で、

- フィボナッチ由来の可変閾値
- √5 系の構造閾値
- √2 系の境界閾値

の **3つが重なりながら交差する** ときに発火するモデルを考える。

---

## 1. 閾値の定義

### フィボナッチ閾値
フィボナッチ数列 \(F_n\) を用いて、

\[
T_F(n) = \frac{F_n}{F_{n+1}}
\]

を可変閾値とする。

これは \(n\) が大きくなると黄金比の逆数

\[
\frac{1}{\varphi}
\]

に収束するが、各ステップではまだ揺らぎを持つ。

### √5 閾値
構造閾値として

\[
T_{\sqrt{5}} = \frac{1}{\varphi}
\]

を置く。

### √2 閾値
別系の境界閾値として

\[
T_{\sqrt{2}} = \frac{1}{\sqrt{2}}
\]

を置く。

---

## 2. 発火の考え方

このモデルでは、発火を単純な

\[
R_n > T
\]

ではなく、**閾値どうしの干渉**として見る。

最小形として、たとえば次のように定義できる。

\[
\mathrm{Fire}(n)=1
\iff
\Big(
R_n > T_F(n)
\Big)
\land
\Big(
\min(|T_F(n)-T_{\sqrt5}|,\ |T_F(n)-T_{\sqrt2}|) < \varepsilon
\Big)
\]

つまり、

1. 残渣がフィボナッチ閾値を超える  
2. その時点でフィボナッチ閾値が √5 または √2 の閾値に十分近い  

という **重なり条件** を満たしたときに発火とみなす。

これは単一閾値モデルではなく、**閾値干渉モデル**である。

---

## 3. 実装例（Python）

```python
import math
import numpy as np
import matplotlib.pyplot as plt

# =========================
# 1. パラメータ
# =========================
N_max = 20
x = np.arange(1, N_max + 1)

# =========================
# 2. フィボナッチ数列
# =========================
F = [0, 1]
for _ in range(2, N_max + 3):
    F.append(F[-1] + F[-2])

# 可変閾値 T_F(n) = F_n / F_{n+1}
T_fib = np.array([F[n] / F[n + 1] for n in range(1, N_max + 1)], dtype=float)

# =========================
# 3. 固定閾値
# =========================
phi = (1 + math.sqrt(5)) / 2
T_sqrt5 = np.full(N_max, 1 / phi)           # ≈ 0.618...
T_sqrt2 = np.full(N_max, 1 / math.sqrt(2))  # ≈ 0.707...

# =========================
# 4. 残渣 R_n
# =========================
# ここではデモ用の残渣を置く。
# 実際には自分のモデルで得た R_n に差し替える。
R = 0.46 + 0.018 * x + 0.045 * np.sin(0.9 * x - 0.7)

# =========================
# 5. 交差判定
# =========================
cross_fib = R > T_fib
cross_sqrt5 = R > T_sqrt5
cross_sqrt2 = R > T_sqrt2

# 発火条件（デモ版）
eps = 0.06
proximity = np.minimum(np.abs(T_fib - T_sqrt5), np.abs(T_fib - T_sqrt2))
fire = (R > T_fib) & (proximity < eps)

# =========================
# 6. プロット
# =========================
plt.figure(figsize=(11, 6))

plt.plot(x, R, marker='o', linewidth=2, label='Residual $R_n$')
plt.step(x, T_fib, where='mid', linewidth=2, label=r'Fibonacci threshold $F_n/F_{n+1}$')
plt.plot(x, T_sqrt5, '--', linewidth=2, label=r'$\sqrt{5}$ threshold $=1/\varphi$')
plt.plot(x, T_sqrt2, '--', linewidth=2, label=r'$\sqrt{2}$ threshold $=1/\sqrt{2}$')

plt.scatter(x[cross_fib], R[cross_fib], marker='x', s=80, label='Cross: Fibonacci')
plt.scatter(x[cross_sqrt5], R[cross_sqrt5], marker='s', s=50, facecolors='none', label='Cross: √5')
plt.scatter(x[cross_sqrt2], R[cross_sqrt2], marker='^', s=60, facecolors='none', label='Cross: √2')
plt.scatter(x[fire], R[fire], marker='*', s=180, label='Fire')

plt.xticks(x)
plt.xlabel('n (integer)')
plt.ylabel('value')
plt.title('Integer-Axis Threshold Interference Model')
plt.grid(True, alpha=0.3)
plt.legend(loc='best')
plt.tight_layout()
plt.show()

print("Fire points:", x[fire].tolist())
