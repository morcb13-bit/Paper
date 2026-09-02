#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
動く図で見るB13講座(その3) の検定。
2+1 次元ミンコフスキー (ct, x, y) で、双曲回転（ローレンツ変換）の合成と、
光円錐の断面の離心率を確かめる。標準ライブラリのみ。

    python3 lorentz.py         検定1〜11 をすべて走らせる
"""

import math

# ---------------------------------------------------------------- 行列（順は t, x, y）

def boost(phi, alpha=0.0):
    c, s = math.cosh(phi), math.sinh(phi)
    nx, ny = math.cos(alpha), math.sin(alpha)
    d = c - 1.0
    return [c,     s * nx,          s * ny,
            s * nx, 1 + d * nx * nx, d * nx * ny,
            s * ny, d * nx * ny,     1 + d * ny * ny]

def galilei(beta):
    """光速を無限大とみなした場合の変換。負の対照に使う。"""
    return [1.0,    0.0, 0.0,
            -beta,  1.0, 0.0,
            0.0,    0.0, 1.0]

def mul(A, B):
    M = [0.0] * 9
    for i in range(3):
        for j in range(3):
            M[i * 3 + j] = sum(A[i * 3 + k] * B[k * 3 + j] for k in range(3))
    return M

def ap(M, v):
    return [M[0] * v[0] + M[1] * v[1] + M[2] * v[2],
            M[3] * v[0] + M[4] * v[1] + M[5] * v[2],
            M[6] * v[0] + M[7] * v[1] + M[8] * v[2]]

def form(v):
    return v[0] * v[0] - v[1] * v[1] - v[2] * v[2]

def decompose(M):
    """M = boost(phi, alpha) * rotation(theta) に分ける。"""
    phi = math.acosh(max(1.0, M[0]))
    sh = math.sinh(phi)
    al = 0.0 if sh < 1e-12 else math.atan2(M[6], M[3])
    R = mul(boost(-phi, al), M)
    return phi, al, math.atan2(R[7], R[4])

# ---------------------------------------------------------------- 断面の離心率

def section_points(beta, a=1.0, n=400, xr=8.0):
    """光円錐 c^2t^2 = x^2+y^2 を平面 ct = beta*x + a で切った点列（未来側）。"""
    pts = []
    for i in range(n + 1):
        x = -xr + 2 * xr * i / n
        D = (beta * beta - 1.0) * x * x + 2 * beta * x + a * a
        ct = beta * x + a
        if D >= 0.0 and ct > 0.0:
            y = math.sqrt(D)
            pts.append((ct, x, y))
            if y > 1e-12:
                pts.append((ct, x, -y))
    return pts

def fit_conic(pts):
    """平面内の直交座標に落とし、Ax^2+Bxy+Cy^2+Dx+Ey+F=0 を最小二乗で当てる。"""
    # 平面 ct = beta*x + a の中の正規直交基底
    e1 = None
    # 基底は呼び出し側で渡さず、点群から主成分的に作る
    cx = sum(p[1] for p in pts) / len(pts)
    cy = sum(p[2] for p in pts) / len(pts)
    ct0 = sum(p[0] for p in pts) / len(pts)
    v1 = None
    for p in pts:
        d = (p[0] - ct0, p[1] - cx, p[2] - cy)
        if d[0] * d[0] + d[1] * d[1] + d[2] * d[2] > 1e-6:
            v1 = d
            break
    nrm = math.sqrt(sum(c * c for c in v1))
    e1 = tuple(c / nrm for c in v1)
    v2 = None
    for p in pts:
        d = (p[0] - ct0, p[1] - cx, p[2] - cy)
        dot = sum(d[i] * e1[i] for i in range(3))
        w = tuple(d[i] - dot * e1[i] for i in range(3))
        if sum(c * c for c in w) > 1e-6:
            v2 = w
            break
    nrm = math.sqrt(sum(c * c for c in v2))
    e2 = tuple(c / nrm for c in v2)
    uv = []
    for p in pts:
        d = (p[0] - ct0, p[1] - cx, p[2] - cy)
        uv.append((sum(d[i] * e1[i] for i in range(3)),
                   sum(d[i] * e2[i] for i in range(3))))
    # F = -1 に固定して 5 元の正規方程式を解く
    rows = [[u * u, u * v, v * v, u, v] for (u, v) in uv]
    N = [[sum(r[i] * r[j] for r in rows) for j in range(5)] for i in range(5)]
    b = [sum(r[i] for r in rows) for i in range(5)]
    for i in range(5):                       # ガウス消去
        piv = max(range(i, 5), key=lambda k: abs(N[k][i]))
        N[i], N[piv] = N[piv], N[i]
        b[i], b[piv] = b[piv], b[i]
        for k in range(i + 1, 5):
            f = N[k][i] / N[i][i]
            for j in range(i, 5):
                N[k][j] -= f * N[i][j]
            b[k] -= f * b[i]
    s = [0.0] * 5
    for i in range(4, -1, -1):
        s[i] = (b[i] - sum(N[i][j] * s[j] for j in range(i + 1, 5))) / N[i][i]
    A, B, C, D, E = s
    return A, B, C, D, E, -1.0

def eccentricity(coef):
    A, B, C, D, E, F = coef
    det = 4 * A * C - B * B
    if abs(det) < 1e-7 * max(1.0, abs(A) + abs(C)):
        return 1.0                                    # 放物線
    x0 = (B * E - 2 * C * D) / det
    y0 = (B * D - 2 * A * E) / det
    F2 = F + (D * x0 + E * y0) / 2
    tr, dd = A + C, A * C - B * B / 4
    disc = math.sqrt(max(0.0, tr * tr / 4 - dd))
    l1, l2 = tr / 2 + disc, tr / 2 - disc
    s1, s2 = -F2 / l1, -F2 / l2                       # 半軸の二乗
    if s1 > 0 and s2 > 0:
        a2, b2 = max(s1, s2), min(s1, s2)
        return math.sqrt(1 - b2 / a2)
    a2 = s1 if s1 > 0 else s2
    b2 = -(s2 if s1 > 0 else s1)
    return math.sqrt(1 + b2 / a2)

def kind(e):
    if e < 1e-6:  return "円"
    if e < 1 - 1e-4:  return "楕円"
    if e < 1 + 1e-4:  return "放物線"
    return "双曲線"

# ---------------------------------------------------------------- 検定

def main():
    ok = lambda b: "OK" if b else "NG"
    print("検定1  不変量 (ct)^2-x^2-y^2 が一歩で変わらない")
    v = [math.cosh(0.7), math.sinh(0.7) * 0.6, math.sinh(0.7) * 0.8]
    worst = max(abs(form(ap(boost(0.13 * i, 0.37 * i), v)) - form(v)) for i in range(1, 40))
    print("       最大のずれ %.3e  →  %s" % (worst, ok(worst < 1e-9)))

    print("検定2  向きが同じなら 合成の歩幅 = φ1+φ2、残る回転 = 0")
    worst = wr = 0.0
    for i in range(1, 40):
        p1, p2, a = 0.06 * i, 0.05 * i, 0.16 * i
        phi, _, rot = decompose(mul(boost(p2, a), boost(p1, a)))
        worst = max(worst, abs(phi - (p1 + p2)))
        wr = max(wr, abs(rot))
    print("       歩幅のずれ %.3e / 回転 %.3e  →  %s" % (worst, wr, ok(worst < 1e-9 and wr < 1e-9)))

    print("検定3  (v1+v2)/(1+v1v2) = tanh(φ1+φ2)")
    worst = max(abs((math.tanh(.07 * i) + math.tanh(.05 * i)) /
                    (1 + math.tanh(.07 * i) * math.tanh(.05 * i)) - math.tanh(.12 * i))
                for i in range(1, 30))
    print("       最大のずれ %.3e  →  %s" % (worst, ok(worst < 1e-12)))

    print("検定4  k = e^φ = √((1+β)/(1−β))")
    worst = max(abs(math.exp(.08 * i) - math.sqrt((1 + math.tanh(.08 * i)) / (1 - math.tanh(.08 * i))))
                for i in range(1, 30))
    print("       最大のずれ %.3e  →  %s" % (worst, ok(worst < 1e-9)))

    print("検定5  光円錐座標で u=ct+x が k 倍、w=ct−x が 1/k 倍")
    p = 0.9
    u = ap(boost(p), [1, 0, 0])
    print("       u %.12f (e^φ %.12f) / w %.12f (e^-φ %.12f)  →  %s"
          % (u[0] + u[1], math.exp(p), u[0] - u[1], math.exp(-p),
             ok(abs((u[0] + u[1]) - math.exp(p)) < 1e-12 and abs((u[0] - u[1]) - math.exp(-p)) < 1e-12)))

    print("検定6  同じ一歩を n 回 → φ = n·φ1")
    v, B1, good = [1, 0, 0], boost(0.8), True
    for k in range(1, 8):
        v = ap(B1, v)
        good = good and abs(math.acosh(v[0]) - 0.8 * k) < 1e-9
    print("       7段まで  →  %s" % ok(good))

    print("検定7  双曲扇形の面積 = φ/2、双曲線に沿った歩幅 = φ")
    def sector(phi, N=200000):
        return sum(0.5 * (math.cosh(phi * (i + .5) / N) ** 2 - math.sinh(phi * (i + .5) / N) ** 2)
                   * (phi / N) for i in range(N))
    def arclen(phi, N=200000):
        s = 0.0
        for i in range(N):
            t = phi * (i + .5) / N
            dt = phi / N
            s += math.sqrt(abs((math.sinh(t) * dt) ** 2 - (math.cosh(t) * dt) ** 2))
        return s
    a17, l17 = sector(1.7), arclen(1.7)
    print("       φ=1.7: 面積 %.8f (φ/2 = %.8f) / 歩幅 %.8f  →  %s"
          % (a17, 0.85, l17, ok(abs(a17 - .85) < 1e-7 and abs(l17 - 1.7) < 1e-7)))

    print("検定8  断面の離心率（点列から当てた円錐曲線） e = β√(2/(1+β²))")
    good = True
    for beta in (0.0, 0.5, 0.9, 1.0, 1.2, 3.0):
        pts = section_points(beta)
        e = eccentricity(fit_conic(pts))
        pred = beta * math.sqrt(2.0 / (1.0 + beta * beta))
        good = good and abs(e - pred) < 2e-4
        print("       β=%4.1f  e=%.6f  式 %.6f  %s" % (beta, e, pred, kind(e)))
    print("       →  %s" % ok(good))

    print("検定9  e=1 になるのは β=1 のときだけ")
    xs = [(b / 1000.0, b / 1000.0 * math.sqrt(2.0 / (1 + (b / 1000.0) ** 2))) for b in range(0, 3001)]
    hits = [b for (b, e) in xs if abs(e - 1.0) < 1e-6]
    print("       e=1 に当たる β: %s  →  %s" % (hits, ok(hits == [1.0])))

    print("検定9b 軸に平行な面 y=a で切ると、切り口は一歩の軌道そのもの")
    a = 1.0
    pts = []
    for i in range(801):
        t = -2.0 + 4.0 * i / 800
        pts.append((a * math.cosh(t), a * math.sinh(t), a))
    e = eccentricity(fit_conic(pts))
    inv = [form([a * math.cosh(t / 10.0), a * math.sinh(t / 10.0), a]) for t in range(-20, 21)]
    print("       e=%.6f（√2 = %.6f） / 面 y=%.1f の上で (ct)²-x²-y² は %s"
          % (e, math.sqrt(2), a, ok(max(inv) - min(inv) < 1e-9)))
    print("       →  %s" % ok(abs(e - math.sqrt(2)) < 2e-4))

    print("検定10 負の対照A：向きを 90° ずらす（φ1=φ2=1）")
    phi, _, rot = decompose(mul(boost(1.0, math.pi / 2), boost(1.0, 0.0)))
    print("       合成の歩幅 %.6f（和は 2.000000）/ 残る回転 %.4f°  →  足し算は崩れる"
          % (phi, rot * 180 / math.pi))

    print("検定11 負の対照B：速度をそのまま足す変換（光速を無限大とみなす）")
    v = [1.0, 0.0, 0.0]
    w = ap(galilei(0.5), v)
    n1 = [1.0, 1.0, 0.0]                     # 光円錐の上の点
    n2 = ap(galilei(0.5), n1)
    print("       不変量 %.4f → %.4f / 光円錐の点 %.4f → %.4f  →  境界が動く"
          % (form(v), form(w), form(n1), form(n2)))

if __name__ == "__main__":
    main()
