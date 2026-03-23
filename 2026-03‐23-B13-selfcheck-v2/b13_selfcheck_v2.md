# B13 図読み取り自己チェック指示書
**貼るだけで動く・完全自己完結版**

---

## Claudeへの指示

この指示書が会話に含まれている場合、図・画像の分析を求められたとき、**言葉での説明の前に必ず以下のコードを実行すること。**  
説明だけで終わることは禁止。数値を出してから解釈する。

---

## Step 1：b13phaseをセットアップ（毎回最初に実行）

```python
# === B13 SETUP ===
import sys, os, zipfile, urllib.request
import numpy as np
from PIL import Image
from collections import defaultdict

B13_PATH = '/tmp/b13phase_env'
ZIP_URL  = 'https://github.com/morcb13-bit/Paper/raw/main/2026-03-23-b13phase-v066/b13phase_v066_optimized.zip'
ZIP_LOCAL = '/tmp/b13phase_v066_optimized.zip'

if not os.path.exists(os.path.join(B13_PATH, 'b13phase')):
    os.makedirs(B13_PATH, exist_ok=True)

    # アップロード済みzipを優先、なければGitHubから取得
    upload_zip = '/mnt/user-data/uploads/b13phase_v066_optimized.zip'
    if os.path.exists(upload_zip):
        src = upload_zip
        print(f"📦 アップロード済みzipを使用")
    else:
        print(f"📥 GitHubからダウンロード中...")
        urllib.request.urlretrieve(ZIP_URL, ZIP_LOCAL)
        src = ZIP_LOCAL

    with zipfile.ZipFile(src, 'r') as z:
        z.extractall(B13_PATH)
    print(f"✅ 展開完了")

sys.path.insert(0, B13_PATH)
from b13phase import COS_TABLE, SIN_TABLE, BASE, STEP_MIN
from b13phase.z13_structure import coset_of
print(f"✅ b13phase ready. BASE={BASE}, STEP_MIN={STEP_MIN}")
```

---

## Step 2：巻き取って読む（パスと半径だけ変える）

```python
# === B13 WINDING ANALYSIS ===

def b13_winding_analysis(image_path, cx=None, cy=None, R=80, label=""):
    """
    image_path : 分析する画像パス
    cx, cy     : 関心領域の中心ピクセル（Noneで画像中央）
    R          : 動径プロファイルの半径（ピクセル）
    label      : 表示用ラベル
    """
    img = Image.open(image_path).convert('L')
    arr = np.array(img).astype(float)
    arr -= arr.mean()

    h, w = arr.shape
    if cx is None: cx = w // 2
    if cy is None: cy = h // 2

    # 動径プロファイル
    Y, X = np.mgrid[-(cy):(h-cy), -(cx):(w-cx)]
    D = np.sqrt(X**2 + Y**2)
    profile = np.array([
        arr[(D >= r) & (D < r+1)].mean() if ((D >= r) & (D < r+1)).sum() > 0 else 0.0
        for r in range(R)
    ])

    # B13巻き取り
    N = len(profile)
    results = []
    for freq_idx in range(1, 61):
        phases = np.array([(freq_idx * r * STEP_MIN) % BASE for r in range(N)])
        cos_vals = np.array([COS_TABLE[int(p)] for p in phases]) / BASE
        sin_vals = np.array([SIN_TABLE[int(p)] for p in phases]) / BASE
        re = np.dot(profile, cos_vals)
        im = np.dot(profile, sin_vals)
        power = re**2 + im**2
        coset = coset_of(freq_idx % 13) if freq_idx % 13 != 0 else '0'
        results.append((freq_idx, power, coset))

    # コセット別集計
    coset_power = defaultdict(float)
    for _, power, coset in results:
        coset_power[coset] += power
    total = sum(coset_power.values())

    h_r  = 100 * coset_power['H']  / total
    h2_r = 100 * coset_power['2H'] / total
    h4_r = 100 * coset_power['4H'] / total
    struct = h_r + h2_r

    verdict = (
        "✅ 本物の構造シグナル" if struct >= 70 else
        "⚠️  ノイズ混入・要注意"  if struct >= 50 else
        "❌ ノイズ優勢"
    )

    print(f"\n{'='*45}")
    print(f"📊 B13巻き取り結果: {label or image_path}")
    print(f"{'='*45}")
    print(f"  H  （構造許容域）: {h_r:5.1f}%")
    print(f"  2H （構造許容域）: {h2_r:5.1f}%")
    print(f"  4H （禁止域）    : {h4_r:5.1f}%")
    print(f"  ───────────────────────────")
    print(f"  構造シグナル合計 : {struct:5.1f}%")
    print(f"  判定             : {verdict}")
    results.sort(key=lambda x: -x[1])
    print(f"\n  TOP5周波数:")
    for freq_idx, power, coset in results[:5]:
        print(f"    freq={freq_idx:3d}  coset={coset:>3}  power={power:,.0f}")
    return struct, h4_r, verdict

# ===== ここだけ変えて使う =====
b13_winding_analysis(
    image_path='YOUR_IMAGE.png',  # ← 画像パス
    cx=None,                       # ← 中心X（Noneで画像中央）
    cy=None,                       # ← 中心Y（Noneで画像中央）
    R=80,                          # ← 動径半径（ピクセル）
    label="解析対象"
)
```

---

## 判定基準

| H+2H合計 | 判定 |
|:---:|:---:|
| ≥ 70% | ✅ 本物の構造シグナル |
| 50〜70% | ⚠️ ノイズ混入 |
| ≤ 50% | ❌ ノイズ優勢 |

**Fib/Pellの4H回避定理（代数証明済み）：本物の構造は4Hを踏まない。**

---

## Claudeの自己チェックリスト

```
□ コードを実行したか？（説明だけで終わっていないか）
□ H / 2H / 4H のパワーを数値で出したか？
□ 判定（✅⚠️❌）を確認したか？
□ 数値を踏まえて図の解釈を述べたか？
```

---

*「巻き取って残渣を見る」― アリでもやっていることだ。*  
*自分で作った図は、自分で数値として読め。*
