#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kabai_ring.py ── 菱形三十面体（RT）を面どうしで繋いで閉じた輪を作る。
                 Kabai–Bérczi 2009 p.100 の再現。

型抜きは使わない。**輪が閉じるかどうかは6次元の整数の足し算で決まる。**

  射影は 6次元 → 3次元（20面体対称）。RT は生成子6本のゾノトープ＝受容窓
  そのもの。RT を面で繋ぐ一歩は、6次元では ±1 を4つ持つ整数ベクトルになる。

  事前登録（判別法6）

  検定1 5回回転が格子の上にあるか
      OK なら：72°回転が6次元で符号つき置換行列 g になり、g^5 = I（整数のまま）
      NG なら：軸の取り違え。以降は全部無効（誤り122の型）

  検定2 一歩が面接触か
      OK なら：30本の歩みそれぞれで、その面の法線だけ支持値ちょうど（＝接する）、
               他14本は支持値未満（＝食い込まない）
      NG なら：一歩の作り方が違う

  検定3 輪が閉じるか
      OK なら：歩みの総和が6次元でちょうど 0（浮動小数を使わない判定）
      NG なら：閉じない

  検定4 重なりゼロ
      OK なら：どの2個も、中心差が 2·RT の内部に入らない
      NG なら：立体が食い込んでいる（dodecahedron-cluster の「代数は閉じるが
               立体が通らない」と同じ型）

  検定5 5回対称
      OK なら：輪全体が g で自分に写る
      NG なら：対称を仮定した組み方が壊れている

  対照   閉じない歩みの組・重なる歩みの組が実際に存在すること（＝検定3・4が
         落ちうること）を数える
"""

import math
import itertools
import json

PHI = (1.0 + math.sqrt(5.0)) / 2.0

# 生成子6本（正20面体の頂点方向。どれも5回対称軸）
GEN = [(1.0, PHI, 0.0), (-1.0, PHI, 0.0),
       (0.0, 1.0, PHI), (0.0, -1.0, PHI),
       (PHI, 0.0, 1.0), (PHI, 0.0, -1.0)]
S = math.sqrt(1.0 + PHI*PHI)
GEN = [tuple(c/S for c in v) for v in GEN]        # 辺長1

def dot(u, v): return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def cross(u, v): return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
def norm(u): return math.sqrt(dot(u, u))
def add3(u, v): return (u[0]+v[0], u[1]+v[1], u[2]+v[2])
def mul3(u, t): return (u[0]*t, u[1]*t, u[2]*t)

def phys(n):
    p = (0.0, 0.0, 0.0)
    for k in range(6):
        if n[k]:
            p = add3(p, mul3(GEN[k], n[k]))
    return p

# ---------------------------------------------------------------- RT の面

def facets():
    """RT の面法線15方向と支持値。面を張る生成子の対も返す。"""
    out = []
    seen = []
    for i, j in itertools.combinations(range(6), 2):
        nv = cross(GEN[i], GEN[j])
        L = norm(nv)
        nv = mul3(nv, 1.0/L)
        if any(abs(abs(dot(nv, m)) - 1.0) < 1e-9 for m in seen):
            continue
        seen.append(nv)
        h = 0.5*sum(abs(dot(GEN[k], nv)) for k in range(6))
        out.append((nv, h, i, j))
    return out

FAC = facets()

def support_ok(c, scale, eps=1e-9):
    """c が scale·RT の内部にあるか（内部なら True＝重なり）。"""
    for nv, h, _, _ in FAC:
        if abs(dot(c, nv)) > scale*h - eps:
            return False
    return True

# ---------------------------------------------------------------- 面で繋ぐ一歩

def steps():
    """RT を面どうしで繋ぐ一歩。6次元の整数ベクトル（±1 が4つ）。"""
    out = []
    for nv, h, i, j in FAC:
        for sgn in (+1, -1):
            v = [0]*6
            for k in range(6):
                if k == i or k == j:
                    continue
                d = dot(GEN[k], nv)
                v[k] = sgn * (1 if d > 0 else -1)
            out.append((tuple(v), (nv, h, i, j), sgn))
    return out

STEP = steps()

# ---------------------------------------------------------------- 5回回転

def rot_about(axis, ang):
    a = mul3(axis, 1.0/norm(axis))
    c, s = math.cos(ang), math.sin(ang)
    def f(p):
        return add3(add3(mul3(p, c), mul3(cross(a, p), s)),
                    mul3(a, dot(a, p)*(1-c)))
    return f

def signed_perm(axis, ang):
    """回転が生成子の集合を符号つきで写すなら、6×6整数行列を返す。"""
    R = rot_about(axis, ang)
    M = [[0]*6 for _ in range(6)]
    for k in range(6):
        q = R(GEN[k])
        hit = None
        for m in range(6):
            if norm((q[0]-GEN[m][0], q[1]-GEN[m][1], q[2]-GEN[m][2])) < 1e-9:
                hit = (m, +1); break
            if norm((q[0]+GEN[m][0], q[1]+GEN[m][1], q[2]+GEN[m][2])) < 1e-9:
                hit = (m, -1); break
        if hit is None:
            return None
        M[hit[0]][k] = hit[1]
    return M

def matvec(M, v):
    return tuple(sum(M[r][c]*v[c] for c in range(6)) for r in range(6))

def matmul(A, B):
    return [[sum(A[r][k]*B[k][c] for k in range(6)) for c in range(6)]
            for r in range(6)]

# ---------------------------------------------------------------- 描画

def rot_matrix(az, el):
    ca, sa = math.cos(az), math.sin(az)
    ce, se = math.cos(el), math.sin(el)
    def f(p):
        x = p[0]*ca - p[1]*sa
        y = p[0]*sa + p[1]*ca
        return (x, y*ce - p[2]*se, y*se + p[2]*ce)
    return f

# 色は「明るさ＝光の当たり方」だけで決める。裏返しに見えないよう、
# 面の向きに対して明るさが単調になるようにする（暗い紺を作らない）。
HUES = [212, 224, 236, 248, 260, 272, 284, 296,
        308, 320, 332, 344, 356, 8, 20]

def hsl(h, s, l):
    c = (1 - abs(2*l - 1)) * s
    x = c * (1 - abs((h/60.0) % 2 - 1))
    m = l - c/2
    r, g, b = [(c, x, 0), (x, c, 0), (0, c, x),
               (0, x, c), (x, 0, c), (c, 0, x)][int(h//60) % 6]
    return "#%02x%02x%02x" % (int((r+m)*255), int((g+m)*255), int((b+m)*255))

def rt_faces(center):
    """RT 1個の30面。各面は (4頂点, 外向き法線, 対の番号)。"""
    out = []
    for idx, (nv, h, i, j) in enumerate(FAC):
        for sgn in (+1, -1):
            n = mul3(nv, sgn)
            f = (0.0, 0.0, 0.0)
            for k in range(6):
                d = dot(GEN[k], n)
                if abs(d) < 1e-9:
                    continue
                f = add3(f, mul3(GEN[k], 0.5 if d > 0 else -0.5))
            f = add3(f, center)
            u, v = mul3(GEN[i], 0.5), mul3(GEN[j], 0.5)
            quad = [add3(add3(f, mul3(u, -1)), mul3(v, -1)),
                    add3(add3(f, u), mul3(v, -1)),
                    add3(add3(f, u), v),
                    add3(add3(f, mul3(u, -1)), v)]
            out.append((quad, n, idx))
    return out

def frame_maps():
    w = GEN[2]
    p = GEN[0]
    d = dot(p, w); p = add3(p, mul3(w, -d)); p = mul3(p, 1.0/norm(p))
    q = cross(w, p)
    def f(v):
        return (dot(v, p), dot(v, q), dot(v, w))
    return f

FRAME = frame_maps()

def draw(centers, path, az, el, shift=(0.0, 0.0, 0.0)):
    R = rot_matrix(az, el)
    # 光は視点のやや左上。裏面は捨ててあるので、見えている面は必ず光を受ける
    light = (-0.28, -0.34, 0.90)
    light = mul3(light, 1.0/norm(light))
    prep = []
    for c in centers:
        for quad, n, idx in rt_faces(c):
            n = FRAME(n)
            quad = [add3(FRAME(p), shift) for p in quad]
            rn = R(n)
            if rn[2] <= 0.02:            # 裏面を捨てる（RT は凸）
                continue
            lam = max(0.0, dot(rn, light))          # 絶対値を取らない
            rq = [R(p) for p in quad]
            depth = sum(p[2] for p in rq)/4.0
            col = hsl(HUES[idx % len(HUES)], 0.30, 0.60 + 0.26*lam)
            prep.append((depth, [(p[0], p[1]) for p in rq], col))
    prep.sort(key=lambda r: r[0])
    xs = [p[0] for _, q, _ in prep for p in q]
    ys = [p[1] for _, q, _ in prep for p in q]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    pad = 0.05*max(x1-x0, y1-y0)
    W = 1100
    H = int(W*(y1-y0+2*pad)/(x1-x0+2*pad))
    # 幅・高さは書かない。viewBox だけにして、置いた枠に合わせて縮む形にする
    L = ['<svg xmlns="http://www.w3.org/2000/svg" '
         'viewBox="%.4f %.4f %.4f %.4f" preserveAspectRatio="xMidYMid meet" '
         'style="width:100%%;height:auto;max-height:80vh">'
         % (x0-pad, -(y1+pad), x1-x0+2*pad, y1-y0+2*pad),
         '<g stroke="#5b5b6a" stroke-width="0.012" stroke-linejoin="round">']
    for _, q, col in prep:
        L.append('<polygon points="%s" fill="%s"/>'
                 % (" ".join("%.4f,%.4f" % (p[0], -p[1]) for p in q), col))
    L.append('</g></svg>')
    open(path, "w").write("\n".join(L))
    return len(prep), W, H

# ---------------------------------------------------------------- 本体

def main():
    print("== 検定1 5回回転が格子の上にあるか ==")
    g = signed_perm(GEN[2], 2*math.pi/5)
    ok1 = g is not None
    if ok1:
        P = [[1 if r == c else 0 for c in range(6)] for r in range(6)]
        for _ in range(5):
            P = matmul(g, P)
        ok1 = all(P[r][c] == (1 if r == c else 0) for r in range(6) for c in range(6))
    print("   g は符号つき置換／g^5 = I → %s" % ("OK" if ok1 else "NG"))
    if not ok1:
        return

    print("== 検定2 一歩が面接触か ==")
    bad = 0
    for v, (nv, h, i, j), sgn in STEP:
        c = phys(v)
        touch = abs(abs(dot(c, nv)) - 2*h) < 1e-9
        inside = all(abs(dot(c, m)) < 2*hm + 1e-9 for m, hm, _, _ in FAC)
        if not (touch and inside):
            bad += 1
    print("   歩み %d 本／不合格 %d → %s" % (len(STEP), bad, "OK" if bad == 0 else "NG"))

    # ---- 輪を探す。1周を5区画に分け、区画あたり m 歩。全体で 5m 個の RT。
    print("\n== 輪の探索（5回対称・区画あたり1〜3歩） ==")
    found = []
    tried = closed_but_overlap = 0
    for m in (1, 2, 3):
        for combo in itertools.product(range(len(STEP)), repeat=m):
            tried += 1
            path = [STEP[k][0] for k in combo]
            # 5区画を g で回して並べる
            cents6 = []
            cur = (0,)*6
            ok = True
            for r in range(5):
                for v in path:
                    w = v
                    for _ in range(r):
                        w = matvec(g, w)
                    cents6.append(cur)
                    cur = tuple(cur[t] + w[t] for t in range(6))
            if any(cur[t] != 0 for t in range(6)):      # 検定3：閉じるか
                continue
            if len(set(cents6)) != 5*m:
                continue
            cs = [phys(c) for c in cents6]
            over = False
            for a in range(len(cs)):
                for b in range(a+1, len(cs)):
                    d = (cs[a][0]-cs[b][0], cs[a][1]-cs[b][1], cs[a][2]-cs[b][2])
                    if support_ok(d, 2.0):              # 検定4：重なり
                        over = True; break
                if over: break
            if over:
                closed_but_overlap += 1
                continue
            found.append((m, combo, cents6, cs))
        print("   区画%d歩まで：試行 %d／閉じたが重なる %d／成立 %d"
              % (m, tried, closed_but_overlap, len(found)))

    print("\n== 検定3・4 ==")
    print("   閉じる歩みの組は存在し、そのうち重なるものも %d 件ある"
          "（＝検定3・4はどちらも落ちうる）" % closed_but_overlap)
    if not found:
        print("   成立 0 件 → NG")
        return

    # 一番小さい輪と、次に大きい輪を採る
    found.sort(key=lambda r: (r[0], r[1]))
    best = {}
    for m, combo, c6, cs in found:
        best.setdefault(5*m, (combo, c6, cs))
    print("   成立した輪の個数（RT の個数ごと）: %s"
          % {k: sum(1 for f in found if 5*f[0] == k) for k in sorted(best)})

    out = {}
    for N in sorted(best):
        combo, c6, cs = best[N]
        # 検定5：輪全体が g で自分に写るか
        st = set(c6)
        P = [0]*6
        for k in combo:
            P = [P[t] + STEP[k][0][t] for t in range(6)]
        P = tuple(P)
        ok5 = all(tuple(matvec(g, c)[t] + P[t] for t in range(6)) in st
                  for c in c6)
        rr = [math.hypot(*[dot(p, GEN[2])*0 + p[t] for t in (0, 1)]) for p in cs]
        print("\n   ── RT %d 個の輪  歩み %s  検定5 %s"
              % (N, combo, "OK" if ok5 else "NG"))
        print("      中心間の最小距離 %.6f（面接触の一歩は %.6f）"
              % (min(norm((cs[a][0]-cs[b][0], cs[a][1]-cs[b][1], cs[a][2]-cs[b][2]))
                     for a in range(len(cs)) for b in range(a+1, len(cs))),
                 norm(phys(STEP[combo[0]][0]))))
        out[N] = {"steps": list(combo), "test5": ok5,
                  "centers6": [list(c) for c in c6]}

    # 描画：成立した輪をすべて
    for N in sorted(best):
      combo, c6, cs = best[N]
      mid = FRAME((sum(c[0] for c in cs)/len(cs),
                   sum(c[1] for c in cs)/len(cs),
                   sum(c[2] for c in cs)/len(cs)))
      shift = mul3(mid, -1.0)

      for tag, az, el in (("iso", 0.55, 1.05), ("top", 0.00, 0.00),
                          ("side", 0.20, 1.5708)):
          name = "kabai_rt%02d_%s.svg" % (N, tag)
          n, W, H = draw(cs, name, az, el, shift)
          print("  %-22s 多角形 %d 枚  %dx%d" % (name, n, W, H))

    pro = obl = 0
    for i, j, k in itertools.combinations(range(6), 3):
        u, v, w2 = GEN[i], GEN[j], GEN[k]
        d = abs(dot(u, cross(v, w2)))
        if d > 0.6:
            pro += 1
        else:
            obl += 1
    print("\n   RT 1個の内訳：尖った菱形六面体 %d ＋ 平たい %d ＝ %d"
          % (pro, obl, pro+obl))

    json.dump({"rings": out}, open("kabai_ring.json", "w"), indent=1)

if __name__ == "__main__":
    main()
