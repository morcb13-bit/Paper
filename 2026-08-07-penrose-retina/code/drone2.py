#  検定FL17  カメラ自身が動く場合（ドローン）
#
#      膜は動かない荷を均す。カメラが対象を追尾すると、対象は画面の中で止まり
#      背景が流れる。すると均されるのは対象のほうで、井戸を作るのは背景になる。
#
#      OK（＝この見立てが正しい）なら：
#          追尾するカメラでは視線が対象から外れ、背景側に捕まる
#      NG なら：見立てが違う。カメラが動いても対象を追える
#
#      対照：同じ動画・同じ切り出し寸法で、カメラを止めた場合。
#      これが通ることは検定FL12 で出ている。

import math, random, json
import numpy as np
from PIL import Image

NFR = 275
CW, CH = 160, 240                      # 切り出し（ドローンの画角）
TAU2 = 2 * math.pi

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
WHITE = json.load(open("truth.json"))["white"]

C = np.array(penrose_cells(5.0, CW, CH))
NC = len(C)
CXi, CYi = C[:, 0].astype(int), C[:, 1].astype(int)

def clampc(cx, cy):
    cx = min(max(cx, CW / 2), 270 - CW / 2)
    cy = min(max(cy, CH / 2), 480 - CH / 2)
    return cx, cy

def make(mode, lag=0):
    """(切り出した明るさの列, 切り出しの中での対象の位置, カメラ中心の列)"""
    lum, tru, cen = [], [], []
    sx, sy = clampc(WHITE[0][0], WHITE[0][1])
    for t in range(NFR):
        if mode == "静止":
            cx, cy = clampc(WHITE[0][0], WHITE[0][1])
        elif mode == "追尾":
            cx, cy = clampc(WHITE[t][0], WHITE[t][1])
        elif mode == "追尾（遅れ）":
            a = 1.0 / max(lag, 1)
            sx += a * (WHITE[t][0] - sx); sy += a * (WHITE[t][1] - sy)
            cx, cy = clampc(sx, sy)
        elif mode == "パン":
            f = t / (NFR - 1)
            cx, cy = clampc(95 + 60 * f, 300 - 130 * f)
        x0, y0 = int(cx - CW / 2), int(cy - CH / 2)
        lum.append(GRAY[t, y0:y0 + CH, x0:x0 + CW][CYi, CXi])
        tru.append((WHITE[t][0] - x0, WHITE[t][1] - y0))
        cen.append((cx, cy))
    return np.array(lum), tru, cen

def run(EXC, tau, gamma, start, sig=12.0, k=200.0, sub=32):
    S = np.zeros(NC)
    x, y = start
    vx = vy = 0.0
    dt = 1.0 / sub
    inv = 1.0 / (sig * sig); cut = 3.0 * sig
    f = math.exp(-gamma * dt)
    path = []
    for t in range(NFR):
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
            x = min(max(x + dt * vx, 0.0), CW - 1.0)
            y = min(max(y + dt * vy, 0.0), CH - 1.0)
        path.append((x, y))
        S += (L - S) / tau
    return path

print("切り出し %dx%d／担体 %d セル／アヒルの幅は約40px" % (CW, CH, NC))
print()
print("判別法2：膜を止めたら結果が消えるか（τ=止め は S≡0＝生の明るさだけで井戸を作る）")
print("%-16s %-5s %-8s %9s %9s %10s" % ("カメラ", "しきい", "τ", "平均ずれ", "最大ずれ", "20px以内"))
for mode, lag in (("静止", 0), ("追尾", 0), ("パン", 0)):
    lum, tru, cen = make(mode, lag)
    for thr in (140.0, 70.0):
        for tau, nm in ((10.0, "10"), (30.0, "30"), (100.0, "100"), (300.0, "300"), (1e9, "止め")):
            p = run(lum > thr, tau, 30.0, tru[0])
            d = [math.hypot(a[0]-b[0], a[1]-b[1]) for a, b in zip(p, tru)]
            print("%-16s %-5.0f %-8s %9.1f %9.1f %9.0f%%"
                  % (mode, thr, nm, sum(d)/len(d), max(d),
                     100*sum(1 for v in d if v < 20)/len(d)))
        print()
