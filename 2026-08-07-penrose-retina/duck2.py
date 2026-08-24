#  検定FL12  実動画でアヒルを追う（照合を羽ごとに直した版）
#  検定FL13  必ず落ちる設定
#
#      FL12  OK なら：視線が一羽の上に乗り、275枚のあいだ離れない。
#                     羽ごとの照合とのずれがアヒルの幅（約40px）以内に収まる
#            NG なら：波紋に捕まる／別の羽へ飛び移る／画面外へ流れる
#      FL13  水面だけの場所から始めると、視線は行き先を持たないこと。
#            持ってしまうなら、追えたのは対象がいたからではない。

import math, random, json
import numpy as np
from PIL import Image, ImageDraw

W_IMG, H_IMG, NFR = 270, 480, 275
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

SPACING, THR = 3.2, 140.0
CELLS = penrose_cells(SPACING, W_IMG, H_IMG)
NC = len(CELLS)
CPOS = np.array(CELLS)
CX = CPOS[:, 0].astype(int); CY = CPOS[:, 1].astype(int)

EXC = np.array([np.asarray(Image.open("fr/%04d.png" % t).convert("L"))
                .astype(np.float32)[CY, CX] > THR for t in range(1, NFR + 1)])
T = json.load(open("truth.json"))
WHITE, GULL = T["white"], T["gull"]
print("担体 %d セル（最近接 %.1f px）／%d 枚／一枚あたりの励起 %.0f セル"
      % (NC, SPACING, NFR, EXC.sum(1).mean()))

SIG = 12.0
CUT = 3.0 * SIG

def run(tau, gamma, start, k=200.0, sub=32):
    S = np.zeros(NC)
    x, y = start
    vx = vy = 0.0
    dt = 1.0 / sub
    inv = 1.0 / (SIG * SIG)
    f = math.exp(-gamma * dt)
    path = []
    for t in range(NFR):
        L = EXC[t].astype(np.float64)
        Wt = L - S
        sel = np.abs(Wt) > 0.05
        pw, pp = Wt[sel], CPOS[sel]
        for _ in range(sub):
            dx = x - pp[:, 0]; dy = y - pp[:, 1]
            r2 = dx * dx + dy * dy
            m = r2 < CUT * CUT
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
    return path

def score(path, truth):
    d = [math.hypot(p[0] - q[0], p[1] - q[1]) for p, q in zip(path, truth) if q]
    return sum(d) / len(d), max(d), sum(1 for v in d if v < 20) / len(d)

print()
print("検定FL12  白いアヒルの上から始める（照合は白いアヒルのみ）")
print("%-6s %-7s %9s %9s %11s %12s" % ("τ", "γ", "平均ずれ", "最大ずれ", "20px以内", "最後の位置"))
best = None
for tau in (30.0, 100.0, 300.0):
    for gamma in (2.0, 8.0, 30.0):
        p = run(tau, gamma, tuple(WHITE[0]))
        m, mx, fr = score(p, WHITE)
        print("%-6.0f %-7.1f %9.1f %9.1f %10.0f%% %12s"
              % (tau, gamma, m, mx, 100 * fr, "(%.0f,%.0f)" % p[-1]))
        if best is None or m < best[0]:
            best = (m, tau, gamma, p)
m, tau, gamma, P = best
print("いちばん良い組：τ=%.0f γ=%.1f 平均ずれ %.1f px（アヒルの幅は約40px）" % (tau, gamma, m))

print()
print("同じ組で別の始点から（τ=%.0f γ=%.1f）" % (tau, gamma))
pg = run(tau, gamma, tuple(GULL[0]))
print("  上のカモメから始める：カモメとのずれ 平均 %.1f / 最大 %.1f、白いアヒルとのずれ 平均 %.1f"
      % (score(pg, GULL)[0], score(pg, GULL)[1], score(pg, WHITE)[0]))

print()
print("検定FL13  水面だけの場所から始める（どの羽からも遠い）")
for st in ((60, 400), (200, 420), (40, 120)):
    q = run(tau, gamma, st)
    disp = math.hypot(q[-1][0] - st[0], q[-1][1] - st[1])
    path = sum(math.hypot(q[i][0] - q[i-1][0], q[i][1] - q[i-1][1]) for i in range(1, NFR))
    dw = min(math.hypot(q[-1][0] - w[0], q[-1][1] - w[1]) for w in (WHITE[-1], GULL[-1]))
    print("  始点(%3d,%3d) → 275枚後 (%3.0f,%3.0f)  始点から %5.1f px / 道のり %6.1f px"
          "  一番近い羽まで %.0f px" % (st[0], st[1], q[-1][0], q[-1][1], disp, path, dw))
print("  （アヒルの上から始めた場合の道のり %.1f px）"
      % sum(math.hypot(P[i][0]-P[i-1][0], P[i][1]-P[i-1][1]) for i in range(1, NFR)))

json.dump({"path": P, "gull_path": pg, "tau": tau, "gamma": gamma}, open("path2.json", "w"))
