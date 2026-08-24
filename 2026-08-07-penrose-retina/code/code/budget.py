#  検定FL15  セルを何枚まで減らせるか
#  検定FL16  何枚に一度読めばよいか
#
#      監視機の代金は電力とデータ量で決まる。装置側で効くのはこの二つだけ。
#      OK なら：画素の数%のセル・数分の1の頻度でも、ずれが対象の幅以内に収まる
#      NG なら：細かい担体と全枚読みでしか成立しない＝監視機の利点が無い

import math, random, json
import numpy as np
from PIL import Image

W_IMG, H_IMG, NFR = 270, 480, 275
TAU2 = 2 * math.pi
NPIX = W_IMG * H_IMG

def penrose_cells(spacing, w, h, seed=3):
    rng = random.Random(seed)
    e = [(math.cos(TAU2 * j / 5), math.sin(TAU2 * j / 5)) for j in range(5)]
    g = [rng.uniform(0.1, 0.9) for _ in range(4)]
    g.append(-sum(g))
    R = math.hypot(w, h) / spacing * 0.75
    K = int(R) + 2
    V = set()
    for j in range(5):
        for l in range(j + 1, 5):
            det = e[j][0] * e[l][1] - e[j][1] * e[l][0]
            for kj in range(-K, K + 1):
                for kl in range(-K, K + 1):
                    a, b = kj - g[j], kl - g[l]
                    px = (a * e[l][1] - b * e[j][1]) / det
                    py = (b * e[j][0] - a * e[l][0]) / det
                    if px * px + py * py > R * R:
                        continue
                    idx = []
                    for p in range(5):
                        if p == j:   idx.append(kj)
                        elif p == l: idx.append(kl)
                        else:        idx.append(math.ceil(px * e[p][0] + py * e[p][1] + g[p]))
                    V.add((round(sum(idx[p] * e[p][0] for p in range(5)), 5),
                           round(sum(idx[p] * e[p][1] for p in range(5)), 5)))
    return [(x * spacing + w / 2, y * spacing + h / 2) for x, y in V
            if 0 <= x * spacing + w / 2 < w and 0 <= y * spacing + h / 2 < h]

GRAY = np.array([np.asarray(Image.open("fr/%04d.png" % t).convert("L")).astype(np.float32)
                 for t in range(1, NFR + 1)])
T = json.load(open("truth.json"))
WHITE = T["white"]

def build(spacing):
    C = np.array(penrose_cells(spacing, W_IMG, H_IMG))
    lum = GRAY[:, C[:, 1].astype(int), C[:, 0].astype(int)]
    return C, lum

def run(C, EXC, tau, gamma, start, sig, k=200.0, sub=32, skip=1):
    NC = len(C)
    S = np.zeros(NC)
    x, y = start
    vx = vy = 0.0
    dt = 1.0 / sub
    inv = 1.0 / (sig * sig)
    cut = 3.0 * sig
    f = math.exp(-gamma * dt)
    path, idx = [], list(range(0, NFR, skip))
    for t in idx:
        L = EXC[t].astype(np.float64)
        Wt = L - S
        sel = np.abs(Wt) > 0.05
        pw, pp = Wt[sel], C[sel]
        for _ in range(sub):
            dx = x - pp[:, 0]; dy = y - pp[:, 1]
            r2 = dx * dx + dy * dy
            m = r2 < cut * cut
            if m.any():
                w = pw[m] * np.exp(-0.5 * r2[m] * inv) * inv
                gx, gy = float((w * dx[m]).sum()), float((w * dy[m]).sum())
            else:
                gx = gy = 0.0
            vx = vx * f - k * gx * dt
            vy = vy * f - k * gy * dt
            x = min(max(x + dt * vx, 0.0), W_IMG - 1.0)
            y = min(max(y + dt * vy, 0.0), H_IMG - 1.0)
        path.append((x, y))
        S += (L - S) / tau
    d = [math.hypot(p[0] - WHITE[i][0], p[1] - WHITE[i][1])
         for p, i in zip(path, idx) if WHITE[i]]
    return sum(d) / len(d), max(d), sum(1 for v in d if v < 20) / len(d)

print("検定FL15  セルを何枚まで減らせるか（しきい値140・τ=300・γ=30）")
print("%-8s %8s %10s %10s %10s %10s" %
      ("間隔px", "セル数", "画素比", "平均ずれ", "最大ずれ", "20px以内"))
KEEP = {}
for sp in (3.2, 5.0, 8.0, 12.0, 18.0):
    C, lum = build(sp)
    KEEP[sp] = (C, lum)
    sig = max(12.0, 2.0 * sp)
    m, mx, fr = run(C, lum > 140.0, 300.0, 30.0, tuple(WHITE[0]), sig)
    print("%-8.1f %8d %9.1f%% %10.1f %10.1f %9.0f%%"
          % (sp, len(C), 100.0 * len(C) / NPIX, m, mx, 100 * fr))
print("  アヒルの幅は約40px。担体の間隔より小さいずれなら、セルの刻みより細かく決まっている")

print()
print("検定FL16  何枚に一度読めばよいか（間隔5.0px・25fps を間引く）")
C, lum = KEEP[5.0]
EXC = lum > 140.0
print("%-10s %8s %10s %10s %10s" % ("間引き", "実効fps", "平均ずれ", "最大ずれ", "20px以内"))
for skip in (1, 2, 5, 10, 25):
    m, mx, fr = run(C, EXC, 300.0 / skip, 30.0, tuple(WHITE[0]), 12.0, skip=skip)
    print("%-10s %8.1f %10.1f %10.1f %9.0f%%"
          % ("1/%d" % skip, 25.0 / skip, m, mx, 100 * fr))
print("  τ は読む枚数で数えるので、間引きに合わせて τ/間引き にしてある")
