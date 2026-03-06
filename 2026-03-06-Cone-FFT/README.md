# cone-FFT / B13 フラクタル位相ライブラリ

**BASE=3120 整数 FFT における黄金比 φ の必然的出現**

---

## 概要

このリポジトリは、整数フラクタル位相表現（B13）を用いた cone-FFT の数学的基盤を収録しています。

核心的な結果：

$$\mu^* = 2\,\mathrm{arcosh}\!\left(\frac{\sqrt{5}+\sqrt{13}}{4}\right) \approx 1.852, \qquad r(\mu^*) = -\varphi^{-2}$$

黄金比 $\varphi = (1+\sqrt{5})/2$ は仮定でなく、整数演算の構造から**必然的に出てくる結果**です。

---

## 哲学

> 余りにこそ価値がある。整数演算のみ。

- 位相は `integer addition` のみで操作する
- 三角関数評価は二次的
- BASE=3120 の選択は偶然ではない

---

## BASE=3120 の必然性

$$3120 = \mathrm{lcm}(2^4,\; 15,\; 13) = 2 \times 5! \times 13$$

| 因子 | 役割 |
|------|------|
| $2^4 = 16$ | 局所相関 $r_{\text{local}}$ の符号周期 |
| $15 = 3 \times 5$ | 5次対称性 × 3段フラクタル階層 |
| $13$ | B13 格子の周期 |

---

## 主な結果

### 1. 閉じた解析式

$$\mu^* = 2\,\mathrm{arcosh}\!\left(\frac{\sqrt{5}+\sqrt{13}}{4}\right) \approx 1.852$$

### 2. 固定点

$$r(\mu^*) = -\varphi^{-2} \approx -0.38197$$

### 3. 因果鎖

```
BASE = lcm(2⁴, 15, 13) = 3120
        ↓
T/U 偶奇分割 → U·n = T·n + n (mod N)
        ↓
r_local の符号 = n mod 16 で完全決定
        ↓
n=0, N/2 が支配項（直流・ナイキスト集中）
        ↓
ρ_a(μ) ≈ cosh(μ) / cosh(μ/2) = √5
        ↓
μ* = 2 arcosh((√5+√13)/4) ≈ 1.852
        ↓
r(μ*) = −φ⁻²   ← φ が結果として出る
```

### 4. 斜めゲート固定点定理

田問題の中心 $(7,7)$（斜めゲート）が、$m=7$ の二重条件を介して T/U 偶奇分割と接続し、$\sqrt{5}$ 固定点を引き起こす。

$$7 \equiv 2 \pmod{5} \;\Rightarrow\; \mathbb{Q}(\sqrt{5})$$
$$7 \equiv -6 \pmod{13} \;\Rightarrow\; \text{B13 境界}$$

---

## ファイル構成

```
cone_fft_theory.md       # 論文（日本語）
cone_fft_theory_en.md    # 論文（英語）
README.md                # このファイル（日本語）
README_en.md             # README（英語）
```

---

## クイックスタート

```python
import math

PHI = (1 + 5**0.5) / 2
BASE = N = 3120

# μ* の解析値
mu_star = 2 * math.acosh((5**0.5 + 13**0.5) / 4)
print(f"mu* = {mu_star:.6f}")  # → 1.852266

# Level-0 テーブル
COS_TABLE = [round(BASE * math.cos(2*math.pi*k/BASE)) for k in range(BASE)]
SIN_TABLE = [round(BASE * math.sin(2*math.pi*k/BASE)) for k in range(BASE)]

# T/U 分割（全偶数 / 全奇数）
T_list = [k for k in range(BASE) if k % 2 == 0]  # 1560 点
U_list = [k for k in range(BASE) if k % 2 == 1]  # 1560 点

# r_local の符号（代数的に決まる）
def r_local_sign(n):
    if n % 16 == 0: return +1   # 像集合一致
    if n % 8  == 0: return -1   # 半周期ずれ
    return None                  # ゼロ点
```

---

## 理論の状態

| 項目 | 状態 |
|------|------|
| $\mu^*$ の閉じた解析式 | ✅ 完了 |
| $r = -\varphi^{-2}$ の代数的証明 | ✅ 完了 |
| $r_{\text{local}}$ 符号の代数的証明 | ✅ 完了 |
| BASE=3120 の必然性 | ✅ 完了 |
| 田問題との接続 | ✅ 完了 |
| fractal 整数版との差 = 0 | ✅ 完了 |
| 論文体裁への整理 | ✅ 完了（草稿） |
| 実装（cone-FFT 本体） | 🔲 今後 |
| 査読投稿 | 🔲 今後 |

---

## ライセンス

未定（査読前草稿）

---

## 変更履歴

| バージョン | 日付 | 内容 |
|-----------|------|------|
| v22 | 2026-03-04 | μ* ≈ 1.852 数値確定 |
| v23 | 2026-03-05 | 因果鎖完成・BASE の lcm 根拠 |
| v24 | 2026-03-06 | μ* 解析式・a(n) 集中・r_local 符号の代数的証明 |
| v25 | 2026-03-06 | cone 固有値証明・田問題接続・理論完成 |
