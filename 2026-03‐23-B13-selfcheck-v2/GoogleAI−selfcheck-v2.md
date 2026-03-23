import numpy as np
from PIL import Image
import math
from collections import defaultdict

def b13_logic_audit(image_path, R=80):
    """
    画像を周波数空間で「巻き取り」、その代数的整合性をスコアリングする
    """
    # 1. 信号抽出：中心からの動径平均プロファイルを取得
    img = Image.open(image_path).convert('L')
    arr = np.array(img).astype(float)
    arr -= arr.mean() # DC成分の除去
    h, w = arr.shape
    cx, cy = w // 2, h // 2

    Y, X = np.mgrid[-cy:h-cy, -cx:w-cx]
    D = np.sqrt(X**2 + Y**2)
    profile = np.array([
        arr[(D >= r) & (D < r+1)].mean() if ((D >= r) & (D < r+1)).sum() > 0 else 0.0
        for r in range(R)
    ])

    # 2. B13 Winding Engine：剰余類解析
    # H = {1, 3, 4, 9, 10, 12}, 2H = {2, 5, 6, 7, 8, 11} (mod 13)
    results = []
    h_coset = {1, 3, 4, 9, 10, 12}
    h2_coset = {2, 5, 6, 7, 8, 11}

    for f in range(1, 61):
        # 複素平面上での回転（巻き取り）によるパワー算出
        theta = 2 * np.pi * f * np.arange(len(profile)) / 13
        re = np.dot(profile, np.cos(theta))
        im = np.dot(profile, np.sin(theta))
        power = re**2 + im**2
        
        coset = 'H' if (f % 13) in h_coset else '2H' if (f % 13) in h2_coset else '4H'
        results.append((f, power, coset))

    # 3. 統計的評価
    stats = defaultdict(float)
    for _, p, c in results: stats[c] += p
    total = sum(stats.values()) or 1.0
    
    struct_sig = 100 * (stats['H'] + stats['2H']) / total
    forbidden_sig = 100 * stats['4H'] / total

    # 4. 判定ロジック：構造か、ただのコピペ（ノイズ）か
    print(f"--- B13 Audit Report ---")
    print(f"Structure Strength (H+2H): {struct_sig:.2f}%")
    print(f"Forbidden Zone (4H):       {forbidden_sig:.2f}%")
    
    if struct_sig >= 70:
        return "✅ PURE STRUCTURE (数学的整合性を検知。論理マップとして詳細解析せよ)"
    elif struct_sig >= 50:
        return "⚠️ MIXED SIGNAL (構造とノイズが混在。図表・加工画像の可能性あり)"
    else:
        return "❌ REAL-WORLD NOISE (代数的構造なし。実写またはコピペ画像と断定し、通常検索へ移行せよ)"

# audit_result = b13_logic_audit('input.png')
# print(audit_result)
