#  検定GS3  跳びの類（実長×残渣）と、静止した輪郭の埋まり方
#
#  対象を一行で（判別法9）：視線2回（0 と d）で静止した輪郭を読んだとき、
#      輪郭の取りこぼしが d の何で決まるか。実長か、残渣か、その積か。
#
#  積は Z[φ] の分岐で二種しかない（b13-verify「実空間と残渣」）
#      五芒星→円環（φ³・φ⁵）… 積 1（単数）
#      連続接続・2つ飛ばし    … 積 √5
#
#  測る量：輪郭の各点から最も近い標本点（世界座標）までの距離。最大と平均。閾値なし。
#  OK なら：積で層に分かれず、実長で並ぶ（＝跳びの類は静止物の読みに効かない）
#  NG なら：積の類で層に分かれる（＝静止物でも跳びの使い分けが要る）
#  必ず落ちる設定：d=0 が改善 0 で最下位に来ること

import math
from collections import defaultdict
import b13_chain_units as U

rows, place, offs = U.build_stack()
CELLS = list(U.fits(sum(place, [])))
XY = {q: U.xy(q) for q in CELLS}
CTR = min(CELLS, key=lambda q: (XY[q][0]) ** 2 + (XY[q][1] - 38.0) ** 2)
OX, OY = XY[CTR]

H = {}
for q in CELLS:
    x, y = XY[q]
    H.setdefault((int(x // 3), int(y // 3)), []).append(q)

def land(px, py):
    gx, gy = int(px // 3), int(py // 3)
    best, bd = None, None
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for q in H.get((gx + dx, gy + dy), ()):
                x, y = XY[q]
                d = (x - px) ** 2 + (y - py) ** 2
                if bd is None or d < bd:
                    bd, best = d, q
    return best

def resid(v):
    re = im = 0.0
    for k in range(4):
        th = math.pi * (3 * k) / 5.0
        re += v[k] * math.cos(th); im += v[k] * math.sin(th)
    return math.hypot(re, im)

def realen(v):
    x, y = U.xy(v)
    return math.hypot(x, y)

def arc(R, n):
    return [(OX + R * math.cos(2 * math.pi * i / n), OY + R * math.sin(2 * math.pi * i / n))
            for i in range(n)]

def poly(R, k, n=40, turn=0.0):
    V = [(OX + R * math.cos(math.radians(90 + turn + 360 * i / k)),
          OY + R * math.sin(math.radians(90 + turn + 360 * i / k))) for i in range(k)]
    out = []
    for i in range(k):
        for j in range(n):
            t = j / n
            out.append((V[i][0] + (V[(i + 1) % k][0] - V[i][0]) * t,
                        V[i][1] + (V[(i + 1) % k][1] - V[i][1]) * t))
    return out

SH = {"丸 R=6.5": arc(6.5, 260), "三角の輪": poly(6.0, 3)}

def samples(shape, gazes):
    pts = set()
    for g in gazes:
        gx, gy = U.xy(g)
        for px, py in shape:
            q = land(px - gx, py - gy)
            pts.add((XY[q][0] + gx, XY[q][1] + gy))
    return list(pts)

def miss(shape, pts):
    mx = 0.0; s = 0.0
    for px, py in shape:
        d = min((px - x) ** 2 + (py - y) ** 2 for x, y in pts) ** 0.5
        mx = max(mx, d); s += d
    return mx, s / len(shape)

ZERO = (0, 0, 0, 0)
CAND = []
seen = set()
for a0 in range(-2, 3):
    for a1 in range(-2, 3):
        for a2 in range(-2, 3):
            for a3 in range(-2, 3):
                v = (a0, a1, a2, a3)
                L = realen(v)
                if L > 6.5:
                    continue
                key = (round(L, 6), round(resid(v), 6), round(U.xy(v)[0], 6))
                if key in seen:
                    continue
                seen.add(key); CAND.append(v)
print("掃く跳び %d 通り（実長6.5以内・重複除く）" % len(CAND))

for nm, sh in SH.items():
    base = miss(sh, samples(sh, [ZERO]))
    print("\n■ %s   視線1回：最大 %.3f / 平均 %.3f" % (nm, base[0], base[1]))
    res = []
    for v in CAND:
        mx, av = miss(sh, samples(sh, [ZERO, v]))
        res.append((realen(v), resid(v), realen(v) * resid(v), mx, av, v))
    res.sort(key=lambda r: r[4])
    print("   最も埋まる10通り")
    print("   %8s %8s %8s %7s %7s  %s" % ("実長", "残渣", "積", "最大", "平均", "跳び"))
    for r in res[:10]:
        print("   %8.4f %8.4f %8.4f %7.3f %7.3f   %s" % (r[0], r[1], r[2], r[3], r[4], r[5]))
    print("   最も埋まらない5通り")
    for r in res[-5:]:
        print("   %8.4f %8.4f %8.4f %7.3f %7.3f   %s" % (r[0], r[1], r[2], r[3], r[4], r[5]))

    # 積の類ごと
    grp = defaultdict(list)
    for L, R, P, mx, av, v in res:
        if L == 0:
            continue
        k = "積 1（単数）" if abs(P - 1.0) < 1e-6 else ("積 √5" if abs(P - 5 ** 0.5) < 1e-6 else "その他 %.2f" % P)
        grp[k].append((L, mx, av))
    print("   積の類ごと（平均取りこぼしの範囲）")
    for k in sorted(grp, key=lambda s: (s.startswith("その他"), s)):
        g = grp[k]
        if k.startswith("その他") and len(g) < 6:
            continue
        print("      %-14s %3d通り  実長 %5.2f〜%5.2f  平均 %.3f〜%.3f"
              % (k, len(g), min(x[0] for x in g), max(x[0] for x in g),
                 min(x[2] for x in g), max(x[2] for x in g)))

    # 実長の帯ごと
    print("   実長の帯ごと（平均取りこぼしの中央）")
    for lo, hi in ((0.0, 0.001), (0.001, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 4.5), (4.5, 5.5), (5.5, 6.5)):
        g = sorted(r[4] for r in res if lo <= r[0] < hi)
        if not g:
            continue
        print("      %4.1f〜%4.1f  %3d通り  中央 %.3f  最小 %.3f  最大 %.3f"
              % (lo, hi, len(g), g[len(g) // 2], g[0], g[-1]))
