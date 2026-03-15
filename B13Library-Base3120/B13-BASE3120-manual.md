# b13phase ライブラリ解説書
**version 0.5.1 / 2026-03-15**  
**coneFFT Verification Instrument — BASE3120**

---

## 目次

1. [概念的背景](#1-概念的背景)
2. [BASE=3120 の構造](#2-base3120-の構造)
3. [Z₁₃の代数構造](#3-z₁₃の代数構造)
4. [切端20面体との対応](#4-切端20面体との対応)
5. [CONEFFTの原点：球と円錐断面](#5-conefffの原点球と円錐断面)
6. [モジュール詳解](#6-モジュール詳解)
7. [クイックスタート](#7-クイックスタート)
8. [確定事項と仮説の区別](#8-確定事項と仮説の区別)
9. [既知の限界](#9-既知の限界)

---

## 1. 概念的背景

### B13とは

B13（Base-13）は、13進数を基底とするフラクタル計算アーキテクチャです。

```
通常の2進数：桁が増えると 2^n の精度
B13：        桁が増えると 13^n の精度
BASE3120：   桁が増えると 3120^n の精度
```

1セルの実効状態空間は Z₁₃（13点）。carry は「保存される状態」ではなく「遷移イベント」です。

### セル設計（確定・変更なし）

```python
DELTA = {0:+1, 1:+2, 2:+3, 3:-1, 4:-2, 5:-3}

def carry_rule(v):
    carry = 0
    while v > +6:  v -= 13;  carry += 1
    while v < -6:  v += 13;  carry -= 1
    return v, carry

def level0_cell(inputs):
    return carry_rule(sum(DELTA[d] for d in inputs))

def level1_cell(r0, c0):    # 解釈B：(r0+c0)合算
    return carry_rule(r0 + c0)
```

**重要**：`level1_cell` の出力 c1 は常に 0。carry c0 は level1 で吸収される。

---

## 2. BASE=3120 の構造

### 素因数分解

```
3120 = 2⁴ × 3 × 5 × 13
     = 16 × 195
     = 60 × 52
     = 切端20面体の頂点数(60) × 活性状態数(52)
```

### 角度単位

BASE=3120 を「1周=3120ステップ」として角度を整数で表現します。

```
3120 / 5   = 624   → 72°  （Z₅の1ステップ、五角形の局所回転）
3120 / 12  = 260   → 30°  （Z₁₂の1ステップ）
3120 / 4   = 780   → 90°  （Hの1ステップ）
3120 / 13  = 240   → 1 Z₁₃ステップ
3120 / 60  = 52    → 6°   （最小単位）
3120 / 360 = 8.67  ← 整数にならない（意図的）
```

360が割り切れないことが B13 の「非自明な角度分割」の根拠です。

### u64への収まり

```
BASE < 2^12 = 4096
→ 1桁 = 12ビット
→ 5桁 = 60ビット → u64に収まる（U64_DIGITS = 5）

5桁の精度：3120⁵ ≈ 2.96 × 10¹⁷
角度分解能：360° / 3120⁵ ≈ 1.22 × 10⁻¹⁵ 度
```

---

## 3. Z₁₃の代数構造

### コセット分割

```
Z₁₃ = {0} ∪ H ∪ 2H ∪ 4H

{0} = {0}              ゼロ元
H   = {1, 5, 8, 12}    3乗余部分群 ≅ Z₄
2H  = {2, 3, 10, 11}   Fibonacci/Pellの到達可能域
4H  = {4, 6, 7, 9}     禁止領域（Fibonacci/Pellが到達しない）
```

### Z₁₃* の巡回順序

生成元 g=2 として：

```
k:  0   1   2   3   4   5   6   7   8   9  10  11
g:  1   2   4   8   3   6  12  11   9   5  10   7
余: H  2H  4H   H  2H  4H   H  2H  4H   H  2H  4H
```

**周期3のパターン：H, 2H, 4H が繰り返す。**

### φ と δₛ の三重対称（確定）

```
φ の余り：8 ∈ H
δₛの余り：5 ∈ H

加法的：8 + 5 = 13 = F(7)
乗法的：8 × 5 ≡ 1 (mod 13)  （互いに逆元）
軌道的：H内で逆向きの同一軌道
  φ軌道：8 → 12 → 5 → 1
  δₛ軌道：5 → 12 → 8 → 1
```

### 4H回避定理（証明済み）

```
定理：F(n) mod 13 は 4H に到達しない。

証明：
  1. M^7 ≡ 8I (mod 13) より F(n+7) ≡ 8·F(n) (mod 13)
  2. 8 ∈ H なので「8倍」は余剰類を保つ
  3. 初期窓 F(0)..F(6) の余剰類 = {0,H,H,2H,2H,H,H} に 4H なし
  4. よって全 n で F(n) ∉ 4H  □

（Pellも同様に証明済み）
```

### 重要な数値事実

```
F(7) = 13 = B13の基底そのもの
P(7) = 169 = 13²
M^7 ≡ 8I (mod 13)
8⁻¹ ≡ 5 (mod 13)
```

---

## 4. 切端20面体との対応

### 基本構造

```
切端20面体：
  頂点：60個
  五角形面：12枚
  六角形面：20枚
  辺：90本

頂点の分解：
  60 = 12（五角形面数）× 5（各五角形の頂点数）
  60 = lcm(|Z₁₃*|, F(5)) = lcm(12, 5)
```

### 頂点写像 f（Proposition D）

```
f: Z₁₃* × Z₅ → V(切端20面体)
f(2^k mod 13, j) = 五角形 P_k の j番頂点

k = 0,...,11  （Z₁₃*の巡回順）
j = 0,...,4   （Z₅の局所番号）
```

### 位相インデックス

```
azimuth_index(k, j) = (k × 260 + j × 624) % 3120
                    = (k × STEP_Z12 + j × STEP_Z5) % BASE
```

### 対面構造（確定）

```
六角形対面：θ = 180°（対蹠点）、ひねりなし（τ=0）
五角形対面：θ ≈ 148°、72°ひねり（τ=1 = Z₅の1ステップ）
```

五角形の対面は36°ずれて向き合い（72°ひねり＝反対側からの展開で正10角形に見える）。

### 余剰類とねじれτの対応

```
2^k ∈ H   → τ=0 遷移を持つ（向き保存可能）
2^k ∈ 2H  → τ=±1 遷移が主
2^k ∈ 4H  → τ=0 が存在しない（常にねじれを伴う）
```

**4H禁止域の幾何的意味：τ=0遷移を持たない五角形層。**

### 52活性頂点

```
52 = Z₁₃ × H = 13 × 4
   = ゼロ点を除く活性状態数
   = BASE / 60 = 3120 / 60
```

---

## 5. CONEFFTの原点：球と円錐断面

### 「球は360°の円錐」

```
円錐の半頂角 α/2 → 90° のとき
断面の離心率 e → 0 （円・球）

α = 180° の円錐 = 球面
= 全方向等距離の定義そのもの

証明：
  e = cos(α/2) / cos(β)  [β = 断面の傾き角]
  e = 0 ⟺ cos(α/2) = 0 ⟺ α = 180°  □
```

### θ分類と円錐断面

球面上の2頂点の中心角 θ：

```
θ < 90°  → 楕円断面（e < 1）← H∪2H（Fib/Pell到達域）
θ = 90°  → 放物線断面（e = 1）← 臨界（carry発火点）
θ > 90°  → 双曲線断面（e > 1）← 4H（禁止域）
θ = 180° → 対蹠点（双曲線極限）← 六角形対面
```

### Observation N（証明済み）

```
命題：切端20面体の全頂点ペアに θ=90° は存在しない。

証明：
  全頂点座標は a+bφ 形式（a,b ∈ Z）の成分を持つ
  任意の2頂点の内積 v_i·v_j = a+bφ
  これがゼロ ⟺ a=b=0（φの無理数性より）
  全ペアで a=b=0 は成立しない  □

系：放物線的臨界は頂点対ではなく辺上にある。
   carry発火の臨界は「通過点」として存在する。
```

### 大円基準線（Observation R）

実物の切端20面体展開図から確認された事実：

```
五角形の頂点2つ＋六角形の1辺を通る大円が存在し
180°で1周する。

1周の構成：h₅×4 + h₆×4 + a×2 （10要素）

φ算術での余弦値（全て分母109で統一）：
  cos(θ_a) = (71+18φ)/109
  cos(θ₅)  = (197-13φ)/218    ← F(7)=13 が出現
  cos(θ₆)  = (54φ-5)/109

109 = 8×F(7)+F(5) = 8×13+5

重要公式：
  (9φ+10)(19-9φ) = 109  （有理化の鍵）
```

### 正10角形断面（Observation P）

```
赤道断面 = 正10角形（10頂点）

対称群 Z₂×Z₅：
  Z₅ ↔ 五角形の局所回転（τの単位）
  Z₂ ↔ carry の符号（±1）

正10角形の10辺 = carry発火の10本の臨界辺
```

### 全方向5回対称（Observation Q）

```
切端20面体の中心Oは全方向に5回回転対称を持つ。

A₅の6本の5回転軸が球面を覆う
→ どの断面方向にも正10角形断面が現れる
→ carry臨界が普遍的構造
→ フラクタル展開の自己相似性の根拠
```

### 重力との接続（仮説）

```
放物線（e=1）= θ=90° = carry発火の臨界
= 束縛軌道（楕円）と脱出軌道（双曲線）の境界

重力 = θ=90°への引力 = carry発火の臨界への引力

切端20面体がθ=90°を頂点として実現しないこと
= B13の系が「重力的束縛」の中に留まる構造的保証
```

---

## 6. モジュール詳解

### constants.py

```python
BASE = 3120              # 基底
STEP_Z5  = 624           # 72°（Z₅の1ステップ）
STEP_Z12 = 260           # 30°（Z₁₂の1ステップ）
STEP_H   = 780           # 90°（Hの1ステップ）
STEP_Z13 = 240           # Z₁₃の1ステップ
STEP_MIN = 52            # 6°（最小単位）

H_SET  = {1, 5, 8, 12}  # 3乗余部分群
H2_SET = {2, 3, 10, 11}
H4_SET = {4, 6, 7, 9}   # 禁止領域

Z13_STAR = [1,2,4,8,3,6,12,11,9,5,10,7]  # 2^k mod 13

F7 = 13                  # F(7) = 13
P7 = 169                 # P(7) = 13²
```

### level0_table.py

```python
# 起動時に math.cos/sin から生成（3120エントリ）
COS_TABLE[i] = round(BASE * cos(2π×i/BASE))
SIN_TABLE[i] = round(BASE * sin(2π×i/BASE))

# 主要値
COS_TABLE[0]    = +3120  # cos(0°)
COS_TABLE[312]  = +2524  # cos(36°) ≈ φ/2 × BASE
COS_TABLE[780]  =     0  # cos(90°)
COS_TABLE[1560] = -3120  # cos(180°)
SIN_TABLE[780]  = +3120  # sin(90°)
```

### phase_digits.py

任意精度のBASE3120桁演算。

```python
digits_from_int(x, n)        # 整数→n桁リスト（MSD→LSD）
digits_to_int(d)             # 桁リスト→整数
digits_add(a, b)             # 加算（carry返却）
digits_inc(d, step=1)        # インクリメント
digits_sub(a, b)             # 減算（borrow返却）
digits_mul_scalar(d, scalar) # スカラー倍
digits_negate(d)             # 符号反転（BASE^n - d）
digits_half(d)               # 1/2（右シフト）
```

### phase_packed_u64.py

5桁u64パック演算（高速・固定精度）。

```python
pack_u64_digits(d)      # 5桁リスト→u64整数
unpack_u64_digits(x)    # u64整数→5桁リスト
add_u64_packed(a, b)    # パック加算（carry返却）
sub_u64_packed(a, b)    # パック減算（borrow返却）
negate_u64_packed(a)    # パック符号反転
```

### evaluator_proto.py

```python
fractal_cos_sin_proto(digits) → (cx, cy)
# digits: MSD→LSD の桁リスト
# 戻り値: (cx, cy) ≈ BASE×(cos θ, sin θ)
# スケール: BASE=3120

norm_sq(cx, cy)         # cx²+cy²（≈ BASE²）
norm_ratio(cx, cy)      # |v|²/BASE²（≈ 1.0）
```

**注意**：protoは近似評価器。本番用は別途実装が必要。

### z13_structure.py

```python
coset_of(g)             # g mod 13 のコセット名 ('0','H','2H','4H')
coset_index(g)          # 0,1,2,3
is_forbidden(g)         # 4Hかどうか
z13_star_index(g)       # 2^k ≡ g の k を返す
coset_pattern()         # ['H','2H','4H','H',...] 周期3
fib_mod13(n)            # F(n) mod 13
pell_mod13(n)           # P(n) mod 13
verify_4h_avoidance(n)  # 4H回避定理の数値検証
phi_residue_orbit()     # [8,12,5,1]（φのH内軌道）
delta_s_orbit()         # [5,12,8,1]（δₛのH内軌道）
```

### truncated_icosahedron.py

```python
R     = 2.478...        # 外接球半径（辺長a=1）
R_SQ  = (9φ+10)/4       # 外接球半径²

vertex_coords_float(k, j) → (x, y, z)
# k=0..11: 五角形インデックス（Z₁₃*巡回順）
# j=0..4:  局所頂点番号（Z₅）

vertex_phase_index(k, j) → int
# = (k×260 + j×624) % 3120

vertex_cos_sin(k, j, n_digits=3) → (cx, cy)
# fractal_cos_sin_protoを使った方位角成分

all_vertices() → List[(k, j, (x,y,z), coset)]
# 全60頂点

inner_product(k1,j1,k2,j2) → float
central_angle(k1,j1,k2,j2) → float  # 度数
classify_pair(k1,j1,k2,j2) → str    # 'elliptic'/'hyperbolic'/'antipodal'
observation_n_verify()      → bool   # θ=90°不在の検証
opposite_face_angles()      → dict   # 対面角度
```

### great_circle.py

```python
great_circle_traversal() → dict
# arc_a, arc_h5, arc_h6, total, theoretical, error_pct

theta_values_degrees() → dict
# θ_a, θ_h5, θ_h6, sum_check（360°に近いか）

phi_arithmetic_verify() → dict
# cos値のφ算術公式の検証

cos_theta_a_phi()  → (71, 18)   # (71+18φ)/109
cos_theta_5_phi()  → (197, -13) # (197-13φ)/218 ← F(7)=13
cos_theta_6_phi()  → (-5, 54)   # (54φ-5)/109

DENOM = 109  # = 8×F(7)+F(5) = 8×13+5
```

---

## 7. クイックスタート

### インストール

```bash
unzip b13phase.zip
cd b13phase_parent_dir
python -m b13phase.demo
```

### 基本的な使い方

```python
from b13phase import (
    BASE, fractal_cos_sin_proto, digits_from_int,
    coset_of, fib_mod13, vertex_coords_float,
    central_angle, great_circle_traversal
)

# 1. 位相評価
digits = digits_from_int(312, 3)   # 36° = 312/3120
cx, cy = fractal_cos_sin_proto(digits)
print(f"cos(36°) ≈ {cx/BASE:.4f}")  # ≈ 0.8090 ≈ φ/2

# 2. Z₁₃コセット判定
for n in range(14):
    f = fib_mod13(n)
    print(f"F({n}) mod 13 = {f:2d}  coset = {coset_of(f)}")

# 3. 頂点座標
x, y, z = vertex_coords_float(k=0, j=0)  # 上極・j=0番頂点
print(f"v(k=0,j=0) = ({x:.4f}, {y:.4f}, {z:.4f})")

# 4. 中心角（θ分類）
theta = central_angle(k1=0, j1=0, k2=6, j2=3)
print(f"θ = {theta:.2f}°")  # > 90° → 双曲線クラス → 4H候補

# 5. 大円解析
gc = great_circle_traversal()
print(f"大円1周誤差: {gc['error_pct']:.4f}%")
print(f"分母109 = {gc['109_factored']}")
```

### 4H回避定理の検証

```python
from b13phase import verify_4h_avoidance
ok = verify_4h_avoidance(n_max=1000)
print("4H回避定理:", "PASS" if ok else "FAIL")
```

### φ算術の確認

```python
from b13phase import phi_arithmetic_verify
result = phi_arithmetic_verify()
for key in ["cos_a", "cos_h5", "cos_h6"]:
    v = result[key]
    print(f"{key}: {v['phi_formula']} → {'✓' if v['match'] else '✗'}")
# cos_h5: (197-13φ)/218  [F(7)=13] → ✓
```

---

## 8. 確定事項と仮説の区別

### 確定事項（証明済み・実測確認済み）

| # | 内容 | 根拠 |
|---|------|------|
| 核1 | F(n),P(n) mod 13 は 4H に到達しない | 数学的証明 |
| 核2 | 8×5≡1(mod 13)、φとδₛの三重対称 | 代数計算 |
| 核3 | θ=90°の頂点ペアは存在しない | Observation N（証明） |
| 核4 | carry はイベントであり保存状態でない | セル設計より |
| 核5 | 六角形対面θ=180°、五角形対面θ≈148° | 幾何計算 |
| 核6 | 大円基準線が h₅×4+h₆×4+a×2 で1周 | 実物展開図による実測 |
| 核7 | 分母109=8×F(7)+F(5)でcos値が統一 | φ算術計算 |
| 核8 | 赤道断面が正10角形 | 幾何計算 |
| 核9 | 中心Oは全方向5回回転対称 | A₅の対称性 |

### 幾何解釈仮説（強い示唆・未証明）

| # | 内容 |
|---|------|
| 仮A | θ<90°↔H∪2H、θ>90°↔4H の完全対応 |
| 仮B | スキルミオン/アンチスキルミオン描像 |
| 仮C | 重力 = θ=90°への引力 = carry発火の臨界への引力 |

### Conjecture D'（予想・検証待ち）

```
臨界断面（大円）はZ₂×Z₅構造を持ち
carry発火の境界を与える

未検証：
  大円が常に10本の辺を横断すること
  10点がcarry臨界点と一致すること
  Z₂×Z₅の群構造が断面に載ること
```

### Proposition D（暫定統合版）

**確定部分**：
1. Fib/Pell は 4H を完全回避（証明済み）
2. θ=90° の頂点ペアは存在しない（証明済み）
3. 放物線的臨界は辺上にある（核3の系）
4. 大円が h₅×4+h₆×4+a×2 で1周（実測確定）

**幾何仮説**：
1. θ<90°（楕円）↔ H∪2H（到達域）
2. θ>90°（双曲線）↔ 4H（禁止域）
3. 臨界断面は正10角形構造を持ちZ₂×Z₅と対応

---

## 9. 既知の限界

### evaluator_proto の近似性

```
fractal_cos_sin_proto は「プロト」実装。
下位桁の補正が近似的で、
完全な数学的精度は保証されない。

本番用には正確な回転行列演算に置き換えること。
```

### vertex_coords_float の近似性

```
頂点座標は float（倍精度）での計算。
θ=90°不在の証明には厳密なφ算術が必要。
（現実装は数値計算による確認のみ）
```

### 大円誤差（1.41%）

```
h₅×4+h₆×4+a×2 の合計が 2πR と
1.41% の誤差がある。

原因：
  - 頂点座標モデルの近似（球面への投影誤差）
  - または通過要素の数え方に漏れがある可能性

実物展開図による実測が正しいので、
理論側の精緻化が次の課題。
```

### f: Z₁₃*×Z₅ → 60頂点 の未完性

```
写像 f の候補は定義されているが：
  - 隣接条件の完全検証が未完
  - twist τ の群論的完全記述が未完
  - 全単射の厳密証明が未完
```

### フラクタル展開（次セッション最優先）

```
level0 → level1 の状態空間拡張と
外接球半径の変化（P(7)=13²との接続）が
未実装・未検証。

漸化式の候補：
  v_{n+1} = φ · R(36°) · v_n
  
cos(36°) = φ/2 がBASE3120で自然に表現できるため
36°回転行列がφ進数の整数係数行列として書ける。
これが確定すれば頂点座標が一意に生成される。
```

---

## 付録：主要な数値事実

```
φ = (1+√5)/2 ≈ 1.6180339887
φ² = φ+1,  φ⁻¹ = φ-1

F(7)  = 13        ← B13の基底
P(7)  = 169 = 13²
M^7 ≡ 8I (mod 13)

R²  = (9φ+10)/4 ≈ 6.1406  （切端20面体外接球半径²）
R   ≈ 2.4780               （辺長a=1のとき）

109 = 8×13+5 = 8×F(7)+F(5) （φ算術の有理化分母）
(9φ+10)(19-9φ) = 109       （有理化の鍵公式）

cos(36°) = φ/2              （BASE3120での表現：2524/3120）
cos(72°) = (φ-1)/2 = 1/(2φ)

3120 = 60×52 = lcm(12,5)×52 = VERTICES × ACTIVE/Z₁₃
```

---

*b13phase v0.5.1 解説書 / 2026-03-15*  
*coneFFT Verification Instrument*
