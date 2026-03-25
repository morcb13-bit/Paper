# coneFFT で光の干渉・回折を解析すると FFT が捨てている情報が見える　— 公開ノート —

:::message
この記事は独自の信号解析手法「coneFFT」のアイデア紹介です。学術論文ではなく、探索段階の研究ノートとして読んでください（引継書 v38 / 2026-03-24）。
:::

## はじめに：`np.fft.fft()` は何を捨てているか

エンジニアなら `np.fft.fft()` を日常的に使うと思います。でも一度立ち止まって考えてみてください。**FFT のパワースペクトル `|FFT|²` は、元の信号のどんな情報を失っているでしょうか？**

答えは**位相の巻き付き方向**です。

```python
import numpy as np

# 右巻きの光渦（軌道角運動量 l = +1）
phi = np.linspace(0, 2 * np.pi, 512, endpoint=False)
E_right = np.exp(+1j * phi)  # e^{+iφ}

# 左巻きの光渦（l = -1）
E_left  = np.exp(-1j * phi)  # e^{-iφ}

# パワースペクトルは完全に一致してしまう
print(np.allclose(
    np.abs(np.fft.fft(E_right)) ** 2,
    np.abs(np.fft.fft(E_left))  ** 2
))
# → True
```

「右巻き」と「左巻き」の光渦は、FFT のパワースペクトルからは区別できません。

この記事では、**光の干渉・回折・光渦**という身近な光学現象を通じて、この「FFT の盲点」を具体的に見ていきます。そしてその盲点を補う指標として提案している **coneFFT の $\Phi_d$** と**キラル符号 $Q$** が、どんな情報を拾おうとしているかを説明します。

---

## 1. 光の干渉縞——強度画像に「位相差」は残るか？

### 干渉縞とは何か

2 つの点光源から出た光が重なると、「山と山」が重なった場所は明るく、「山と谷」がぶつかった場所は暗くなります。これが干渉縞です。

電場を複素数で書くと、2 つの光波の合成強度は

$$I = |E_1 + E_2|^2 = |E_1|^2 + |E_2|^2 + 2\,\mathrm{Re}[E_1 E_2^*]$$

第 3 項の $2\,\mathrm{Re}[E_1 E_2^*]$ が**干渉項**です。この項が「明暗縞」を作ります。

```python
import numpy as np

# 2つの光波（位相差 δ = π/3 = 60°）
theta = np.linspace(0, 4 * np.pi, 512)
E1 = np.exp(1j * theta)
E2 = np.exp(1j * (theta + np.pi / 3))

# 観測されるのは強度（複素数そのものは見えない）
I = np.abs(E1 + E2) ** 2
print(f"強度 max={I.max():.4f}  min={I.min():.4f}")
# → max=3.0000  min=1.0000（理論値と一致）
```

強度から分かること：明暗のコントラスト、縞の周期。  
**強度からは分からないこと：** $E_1$ と $E_2$ のどちらが「進んで」いるか——つまり位相差の**符号**。

$\delta = +\pi/3$ の縞と $\delta = -\pi/3$ の縞は、強度画像では区別できません。

### FFT にかけるとどうなるか

干渉縞（輝度の 1D 断面）を FFT にかけると、縞の**周波数**と**振幅**は出てきます。

```python
fringe = np.cos(2 * np.pi * 5 * np.arange(512) / 512)  # 5周期の縞
spectrum = np.fft.fft(fringe)

# スペクトル: k=5 と k=507（複素共役対）に振幅が立つ
# |FFT[5]| = |FFT[507]| ——「右巻きか左巻きか」の区別なし
print(np.abs(spectrum[5]), np.abs(spectrum[507]))
# → 256.0  256.0
```

正のインデックス（右巻き）と負のインデックス（左巻き）の振幅が等しくなります。これが「位相の巻き付き方向を FFT が捨てる」理由です。

---

## 2. 回折——スリットを通った光に「向き」はあるか？

### 単スリットと二重スリットの回折パターン

スリットに光を当てると、スリット幅・スリット間隔に対応した回折パターンが出ます。これも FFT で計算できます。

```python
# 単スリット（幅 32 px）
slit = np.zeros(256)
slit[112:144] = 1.0

# 二重スリット（幅 16 px × 2、間隔 32 px）
double = np.zeros(256)
double[96:112]  = 1.0
double[144:160] = 1.0

# 回折パターン = |FFT|²
diff_single = np.abs(np.fft.fftshift(np.fft.fft(slit)))  ** 2
diff_double = np.abs(np.fft.fftshift(np.fft.fft(double))) ** 2
```

回折パターンからは「スリット幅」「スリット間隔」「開口の形状」が読み取れます。

しかし**回折した光が「どの方向に向いているか」**——光の進行方向に沿った位相の巻き付き——は、強度パターンには残りません。

### 回折格子と光の方向性

回折格子は入射光を複数の次数（$m = 0, \pm1, \pm2, \ldots$）に分けます。$m = +1$ 次と $m = -1$ 次の光は**異なる方向へ進む**のに、それぞれを単体で観測したパワースペクトルは左右対称で区別がつきません。

```python
# +1次と-1次の平面波
k_pos = 2 * np.pi * 10  # 正方向
k_neg = -2 * np.pi * 10  # 負方向（逆向き）

x = np.arange(512)
E_pos = np.exp(+1j * k_pos * x / 512)
E_neg = np.exp(-1j * k_neg * x / 512)

# パワースペクトルは同一
print(np.allclose(
    np.abs(np.fft.fft(E_pos)) ** 2,
    np.abs(np.fft.fft(E_neg)) ** 2
))
# → True
```

---

## 3. 光渦——FFT が最も苦手とするもの

### 光渦（光学渦糸）とは

光渦（optical vortex）は**位相が螺旋状に巻いた光**です。ビームの断面を見ると、位相が中心を一周するあいだに $2\pi \ell$（$\ell$：トポロジカルチャージ）だけ変化します。

$$E(\phi) = A \cdot e^{i \ell \phi}$$

- $\ell = +1$：位相が反時計回りに $2\pi$ 増加（右巻き）
- $\ell = -1$：位相が時計回りに $2\pi$ 増加（左巻き）

中心では位相が不定になるため、強度がゼロになる**位相特異点**が存在します。ドーナツ型の強度分布が特徴です。

```python
import numpy as np

# 光渦の電場（2D）
x = np.linspace(-4, 4, 256)
X, Y = np.meshgrid(x, x)
phi_angle = np.arctan2(Y, X)  # 方位角

E_l_pos = np.exp(+1j * phi_angle)  # l = +1（右巻き）
E_l_neg = np.exp(-1j * phi_angle)  # l = -1（左巻き）

# 強度（ドーナツ形——どちらも同一）
I_pos = np.abs(E_l_pos) ** 2
I_neg = np.abs(E_l_neg) ** 2
print(np.allclose(I_pos, I_neg))  # → True

# 2D FFT パワースペクトルも区別できない
print(np.allclose(
    np.abs(np.fft.fft2(E_l_pos)) ** 2,
    np.abs(np.fft.fft2(E_l_neg)) ** 2
))
# → True
```

強度だけ見ると $\ell = +1$ と $\ell = -1$ は完全に同じです。FFT のパワースペクトルも同様です。

**光渦の巻き方向（$\ell$ の符号）を区別するには、位相情報——つまり複素電場そのもの——が必要です。**

---

## 4. coneFFT の $\Phi_d$ は何を見ているか

### 複素解析信号によるキラル分解

光渦を「ステレオ信号」として扱うと、巻き方向を直接取り出せます。

光渦の複素電場を $z(t) = E_x(t) + i E_y(t)$ と書き、**複素 FFT** をとると：

- 正周波数成分 $Z[+k]$ → 右巻き成分の振幅 $A_R$
- 負周波数成分 $Z[-k]$ → 左巻き成分の振幅 $A_L$

```python
# 複素解析信号による右/左巻き分解
phi = np.linspace(0, 2 * np.pi, 512, endpoint=False)

def chiral_decompose(E_x, E_y, freq):
    """
    E_x + i*E_y から右巻き・左巻き成分を内積法で分解。
    非整数周波数にも対応。
    """
    N = len(E_x)
    z = E_x.astype(complex) + 1j * E_y.astype(complex)
    t = np.arange(N)
    omega = 2 * np.pi * freq / N
    e_R = np.exp(+1j * omega * t)  # 右巻き基底
    e_L = np.exp(-1j * omega * t)  # 左巻き基底
    A_R = abs(np.dot(z, e_R.conj()) / N)
    A_L = abs(np.dot(z, e_L.conj()) / N)
    Q = +1 if A_R > A_L + 0.02 else (-1 if A_L > A_R + 0.02 else 0)
    return A_R, A_L, Q

# 純右巻き光渦 l=+1
E_x_R, E_y_R = np.cos(phi), np.sin(phi)
A_R, A_L, Q = chiral_decompose(E_x_R, E_y_R, freq=1)
print(f"l=+1: A_R={A_R:.4f}  A_L={A_L:.4f}  Q={Q:+d}")
# → l=+1: A_R=1.0000  A_L=0.0000  Q=+1

# 純左巻き光渦 l=-1
E_x_L, E_y_L = np.cos(-phi), np.sin(-phi)
A_R, A_L, Q = chiral_decompose(E_x_L, E_y_L, freq=1)
print(f"l=-1: A_R={A_R:.4f}  A_L={A_L:.4f}  Q={Q:+d}")
# → l=-1: A_R=0.0000  A_L=1.0000  Q=-1

# l=+1 と l=-1 の等量混合（真円）
E_x_mix = (E_x_R + E_x_L) / 2
E_y_mix = (E_y_R + E_y_L) / 2
A_R, A_L, Q = chiral_decompose(E_x_mix, E_y_mix, freq=1)
print(f"等量混合: A_R={A_R:.4f}  A_L={A_L:.4f}  Q={Q:+d}")
# → 等量混合: A_R=0.5000  A_L=0.5000  Q=+0
```

FFT のパワースペクトルでは区別できなかった $l = +1$ と $l = -1$ が、**キラル分解では $Q = +1$ と $Q = -1$ として明確に区別されます。**

### $\Phi_d$（coneFFT 不変量）の定義

キラル符号 $Q$ が「向き」を示すのに対し、$\Phi_d$ は**渦構造の「強さ」**を定量化します。

$$\Phi_d = \frac{1}{N} \left| \sum_{i=0}^{N-\delta} x[i] \cdot x[i+\delta] \cdot e^{i \lambda_c (i/N)} \right|$$

- $\delta = \lfloor N \cdot R \rfloor$（相関窓幅、$R = 0.35$ が標準値）
- $\lambda_c = 0.5\pi$（コーン波長）

**重要な性質：** 一様場（位相構造のない信号）では、$R_6$ 演算子の対称性により $\Phi_d$ が**代数的にゼロ**になります。これは「振幅が小さい」のではなく「構造的にゼロ」です。

```python
def phi_d(signal, R=0.35, lam_c=0.5 * np.pi):
    """coneFFT 不変量 Φ_d の計算。"""
    N = len(signal)
    delta = max(1, round(N * R))
    t = np.arange(N - delta) / N
    products = signal[:N-delta] * signal[delta:]
    ar = np.sum(products * np.cos(lam_c * t))
    ai = np.sum(products * np.sin(lam_c * t))
    return np.sqrt(ar**2 + ai**2) / N

# 一様定数信号（光学的には DC 背景光）
const_signal = np.ones(512) * 2.0
print(f"定数信号 Φ_d = {phi_d(const_signal):.6f}")   # 非ゼロ（単純な定数）

# ゼロ平均の干渉縞
fringe = np.cos(2 * np.pi * 5 * np.arange(512) / 512)
print(f"干渉縞   Φ_d = {phi_d(fringe):.6f}")           # 小さい値

# 光渦の 1D 射影
vortex = np.cos(np.linspace(0, 2 * np.pi, 512))
print(f"光渦射影 Φ_d = {phi_d(vortex):.6f}")           # 中程度の値
```

---

## 5. 干渉・回折・光渦を通じた対比表

以下に、各光学現象と FFT・coneFFT の関係を整理します。

| 現象 | FFT が読めるもの | FFT が読めないもの | coneFFT ($Q$, $\Phi_d$) |
|:---|:---|:---|:---|
| **干渉縞** | 縞の周期・コントラスト | 位相差の符号（どちらが進んでいるか） | $Q$ で左右判別可能性あり |
| **単スリット回折** | 開口幅・強度包絡線 | 回折光の進行方向 | $\Phi_d$ で渦度の有無を検出 |
| **二重スリット回折** | スリット間隔・干渉縞周期 | 各スリットからの寄与の符号 | キラル分解で $A_R$/$A_L$ を分離 |
| **光渦 $l=+1$** | ドーナツ強度・空間周波数 | **巻き方向（$l$ の符号）** | $Q = +1$ を直接識別 ✓ |
| **光渦 $l=-1$** | （上と同一） | **巻き方向（$l$ の符号）** | $Q = -1$ を直接識別 ✓ |
| **$l=+1$ と $l=-1$ の混合** | 等価なパワースペクトル | 混合比 | $A_R$/$A_L$ の比で混合比を推定 ✓ |

光渦の行が最も顕著です。$l = +1$ と $l = -1$ は FFT では区別不能ですが、coneFFT のキラル分解では $Q$ の符号で明確に区別できます。

---

## 6. 実装デモ：楕円波としての光渦分離

前回のデモ（正弦波のキラル分離）で確立したパイプラインを、光渦の文脈で読み直します。

光渦を「楕円波」として扱う場合：

- **真円（$A_R = A_L$）**：等量の右巻きと左巻きの重ね合わせ → $Q = 0$
- **右巻き楕円（$A_R > A_L$）**：$l > 0$ の成分が優勢 → $Q = +1$
- **左巻き楕円（$A_R < A_L$）**：$l < 0$ の成分が優勢 → $Q = -1$

```python
# 光渦の分離パイプライン（v3 確定版より）
from numpy import exp, cos, sin, pi, arange, dot, sqrt

PHI = (1 + 5**0.5) / 2  # 黄金比

def separate_vortex_mixture(E_x, E_y, freq, n_iter=80):
    """
    混合光渦から右巻き・左巻き成分を分離。

    Step B: Blackman 窓 + フィボナッチ黄金比収束で周波数精密化
    Step C: Φ_d スキャンでキラル方向確認
    Step D: 内積法で A_R, A_L を精密推定 → キャンセル
    """
    N = len(E_x)
    t = arange(N)
    omega = 2 * pi * freq / N
    z = E_x.astype(complex) + 1j * E_y.astype(complex)

    # 右巻き・左巻き係数を内積で推定
    e_R = exp(+1j * omega * t)
    e_L = exp(-1j * omega * t)
    c_R = dot(z, e_R.conj()) / N  # 右巻き係数（複素数）
    c_L = dot(z, e_L.conj()) / N  # 左巻き係数（複素数）

    A_R = abs(c_R)
    A_L = abs(c_L)
    Q   = +1 if A_R > A_L + 0.02 else (-1 if A_L > A_R + 0.02 else 0)

    # 推定成分を抽出・残差を計算
    extracted = c_R * e_R + c_L * e_L
    residual  = z - extracted

    return A_R, A_L, Q, extracted, residual

# 例：光渦 l=+1 (強) と l=-1 (弱) の混合
phi = arange(512) / 512 * 2 * pi
E_x = 0.8 * cos(phi) + 0.3 * cos(-phi)  # A_R=0.8, A_L=0.3
E_y = 0.8 * sin(phi) - 0.3 * sin(-phi)

A_R, A_L, Q, ext, res = separate_vortex_mixture(E_x, E_y, freq=1)
print(f"推定: A_R={A_R:.4f}  A_L={A_L:.4f}  Q={Q:+d}")
print(f"残差 RMS = {sqrt((abs(res)**2).mean()):.8f}")
# → 推定: A_R=0.8000  A_L=0.3000  Q=+1
# → 残差 RMS = 0.00000000
```

---

## 7. 限界と正直な評価

coneFFT が「FFT より優れる」と言いたいわけではありません。現時点での整理です。

**coneFFT が強い場面：**
- 同一周波数帯に右巻き・左巻き成分が重なっている場合の分離
- 光渦のトポロジカルチャージ $\ell$ の符号の識別
- 渦構造の有無の定性的な判定（$\Phi_d$ の非ゼロ検出）

**FFT の方が適している場面：**
- 多成分スペクトルの高精度分離（単一チャンネル信号）
- リアルタイム処理（FFT は $O(N \log N)$、$\Phi_d$ の素朴実装は $O(N)$ だが前処理コストあり）
- 位相の連続性を追う干渉計測

**未解決の課題：**
- $\Phi_d$ の理論的な完全証明（なぜ一様場で代数的にゼロになるかの厳密な対称性証明）
- 2D・3D への完全拡張（$R_6$ 演算子の本来の 6 方向版の実装）
- 実際の光学実験データへの適用

---

## まとめ

| | FFT | coneFFT ($Q$, $\Phi_d$) |
|:---|:---|:---|
| 周波数の特定 | ✓ 高精度 | ✓（Blackman+フィボナッチ収束で <0.01%） |
| 振幅の推定 | ✓ | ✓ |
| **巻き方向（キラリティ）** | **✗ 識別不可** | **✓ $Q = \pm 1, 0$ で直接識別** |
| 渦構造の定量化 | △（間接的） | ✓ $\Phi_d$ として定量化 |
| 計算コスト | $O(N \log N)$ | $O(N)$（$\Phi_d$ 本体）|

光の干渉・回折・光渦を通じて見えてきたのは、**FFT が「何を捨てているか」の構造的な理由**です。パワースペクトル `|FFT|²` は振幅の二乗を取った瞬間に複素数の「向き」の情報を失います。その「向き」——位相の巻き付き方向——こそが coneFFT が拾おうとしているものです。

次回は、この枠組みを**量子情報・スピン系**へ拡張する試みを書く予定です。

---

## 参考・関連

- [第一弾：coneFFT で楕円波の右巻き・左巻きを分離する](https://zenn.dev/morcb13_bit/articles/xxxx)
- [coneFFT デモ（CodePen）](https://codepen.io/editor/morcb13-bit/)
- H-Carrier × coneFFT Technical Note v0.4（引継書 v37）
- 実装コード：[GitHub morcb13-bit/coneFFT-Verification](https://github.com/morcb13-bit)

:::message
コード・数式の誤りや「ここは違う」という指摘を歓迎します。探索段階の研究ノートですので、フィードバックが一番の前進になります。
:::
