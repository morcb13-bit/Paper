#  検定FL18  静止した明るい物（日に温まった岩に相当）があるとき
#
#      分担：ホバリングして「気づく」段は膜つき（動く熱源だけを井戸にする）。
#            掛かったら「追う」段は膜を止める（生の井戸で足りる。カメラが動いてよい）。
#
#      OK なら：膜つきでは静止した明るい物が井戸にならず、動くものだけに掛かる。
#               膜を止めると静止した物に捕まる（＝二段に分ける理由がある）
#      NG なら：どちらも同じ。分ける理由が無い
#
#      必ず落ちる設定：動くものを消し、静止した明るい物だけ置いた画面では
#                     膜つきの視線が行き先を持たないこと。

import math, random, json
import numpy as np
from PIL import Image

W_IMG, H_IMG, NFR = 270, 480, 275
TAU2 = 2 * math.pi

def penrose_cells(spacing, w, h, seed=3):
    rng = random.Random(seed)
    e = [(math.cos(TAU2 * j / 5), math.sin(TAU2 * j / 5)) for j in range(5)]
    g = [rng.uniform(0.1, 0.9) for _ in range(4)]; g.append(-sum(g))
    R = math.hypot(w, h) / spacing * 0.75; K = int(R) + 2
    V = set()
    for j in range(5):
        for l in range(j + 1, 5):
            det = e[j][0] * e[l][1] - e[j][1] * e[l][0]
            for kj in range(-K, K + 1):
                for kl in range(-K, K + 1):
                    a, b = kj - g[j], kl - g[l]
                    px = (a * e[l][1] - b * e[j][1]) / det
                    py = (b * e[j][0] - a * e[l][0]) / det
                    if px * px + py * py > R * R: continue
                    idx = []
                    for p in range(5):
                        if p == j:   idx.append(kj)
                        elif p == l: idx.append(kl)
                        else:        idx.append(math.ceil(px*e[p][0] + py*e[p][1] + g[p]))
                    V.add((round(sum(idx[p]*e[p][0] for p in range(5)), 5),
                           round(sum(idx[p]*e[p][1] for p in range(5)), 5)))
    return [(x*spacing + w/2, y*spacing + h/2) for x, y in V
            if 0 <= x*spacing + w/2 < w and 0 <= y*spacing + h/2 < h]

GRAY = np.array([np.asarray(Image.open("fr/%04d.png" % t).convert("L")).astype(np.float32)
                 for t in range(1, NFR + 1)])
WHITE = json.load(open("truth.json"))["white"]

ROCK = (140.0, 205.0)          # 静止した明るい物（アヒルの通り道の脇）
RR = 13.0

def with_rock(hide_duck=False):
    G = GRAY.copy()
    yy, xx = np.mgrid[0:H_IMG, 0:W_IMG]
    m = (xx - ROCK[0]) ** 2 + (yy - ROCK[1]) ** 2 < RR * RR
    if hide_duck:
        for t in range(NFR):                      # アヒルを水の明るさで塗り潰す
            dx, dy = WHITE[t]
            md = (xx - dx) ** 2 + (yy - dy) ** 2 < 26 ** 2
            G[t][md] = 55.0
    G[:, m] = 225.0
    return G

C = np.array(penrose_cells(5.0, W_IMG, H_IMG))
NC = len(C)
CXi, CYi = C[:, 0].astype(int), C[:, 1].astype(int)

def run(EXC, tau, gamma, start, sig=12.0, k=200.0, sub=32):
    S = np.zeros(NC); x, y = start; vx = vy = 0.0
    dt = 1.0 / sub; inv = 1.0 / (sig*sig); cut = 3.0*sig
    f = math.exp(-gamma*dt); path = []
    for t in range(NFR):
        L = EXC[t].astype(np.float64); Wt = L - S
        sel = np.abs(Wt) > 0.05; pw, pp = Wt[sel], C[sel]
        for _ in range(sub):
            dx = x - pp[:, 0]; dy = y - pp[:, 1]
            r2 = dx*dx + dy*dy; m = r2 < cut*cut
            if m.any():
                w = pw[m] * np.exp(-0.5*r2[m]*inv) * inv
                gx, gy = float((w*dx[m]).sum()), float((w*dy[m]).sum())
            else: gx = gy = 0.0
            vx = vx*f - k*gx*dt; vy = vy*f - k*gy*dt
            x = min(max(x + dt*vx, 0.0), W_IMG-1.0)
            y = min(max(y + dt*vy, 0.0), H_IMG-1.0)
        path.append((x, y)); S += (L - S) / tau
    return path

def stat(p):
    dd = [math.hypot(a[0]-b[0], a[1]-b[1]) for a, b in zip(p, WHITE)]
    dr = [math.hypot(a[0]-ROCK[0], a[1]-ROCK[1]) for a in p]
    return (sum(dd)/len(dd), sum(dr)/len(dr),
            sum(1 for i in range(NFR) if dd[i] < dr[i]) / NFR)

GR = with_rock()
print("静止した明るい物を (%.0f,%.0f) 半径%.0f に置く。アヒルは (%.0f,%.0f)→(%.0f,%.0f)"
      % (ROCK[0], ROCK[1], RR, WHITE[0][0], WHITE[0][1], WHITE[-1][0], WHITE[-1][1]))
print("担体 %d セル（間隔5.0px）／しきい値140" % NC)
print()
print("検定FL18a  アヒルの上から始める")
print("%-10s %12s %12s %14s" % ("膜", "アヒルとのずれ", "岩とのずれ", "アヒル側にいた枚"))
EXC = GR[:, CYi, CXi] > 140.0
for tau, nm in ((10.0, "τ=10"), (30.0, "τ=30"), (100.0, "τ=100"), (1e9, "止め")):
    a, b, f = stat(run(EXC, tau, 30.0, tuple(WHITE[0])))
    print("%-10s %12.1f %12.1f %13.0f%%" % (nm, a, b, 100*f))

print()
print("検定FL18b  岩の上から始める（気づく段。どちらへ行くか）")
for tau, nm in ((10.0, "τ=10"), (30.0, "τ=30"), (1e9, "止め")):
    p = run(EXC, tau, 30.0, ROCK)
    a, b, f = stat(p)
    print("%-10s 275枚後 (%.0f,%.0f)  アヒルとのずれ %.1f  岩とのずれ %.1f"
          % (nm, p[-1][0], p[-1][1], a, b))

print()
print("検定FL18c  必ず落ちる設定（アヒルを消し、静止した明るい物だけ）")
EXC2 = with_rock(hide_duck=True)[:, CYi, CXi] > 140.0
for tau, nm in ((10.0, "τ=10"), (30.0, "τ=30"), (1e9, "止め")):
    for st in ((ROCK[0]+25, ROCK[1]+25), (90, 300)):
        p = run(EXC2, tau, 30.0, st)
        path = sum(math.hypot(p[i][0]-p[i-1][0], p[i][1]-p[i-1][1]) for i in range(1, NFR))
        print("  %-6s 始点(%3.0f,%3.0f) → (%3.0f,%3.0f)  道のり %6.1f px  岩まで %.0f px"
              % (nm, st[0], st[1], p[-1][0], p[-1][1], path,
                 math.hypot(p[-1][0]-ROCK[0], p[-1][1]-ROCK[1])))
