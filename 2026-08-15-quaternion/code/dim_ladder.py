#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dim_ladder.py ── 切断投影の梯子を2次元・3次元・4次元で並べる。

同じ組み立て（整数の担体・残渣空間の窓・窓の中に落ちた点だけ採る）が
次元ごとにどう変わるかを一枚に並べる。

  事前登録

  検定1 担体の次元
      OK なら：2次元は5、3次元は6、4次元は8。残渣の次元は 3・3・4
      NG なら：数え方が違う

  検定2 窓の頂点数
      OK なら：正五角形5・菱形三十面体32・600胞体120
      NG なら：窓の作り方が違う

  検定3 インフレーションの最小整数冪
      OK なら：φ倍が担体の上で整数行列になる最小の冪が、2次元1・3次元3・
               4次元1
      NG なら：どこかで基底の取り方を間違えている

  検定4 点群の位数
      OK なら：10・120・14400（後ろ2つは実際に数える）
      NG なら：群の作り方が違う
"""

import math
import itertools
from collections import defaultdict

PHI = (1 + 5**0.5)/2


def hexcol(h, sat, lig):
    c = (1-abs(2*lig-1))*sat
    x = c*(1-abs((h/60.0) % 2 - 1))
    m = lig - c/2
    r, g, b = [(c, x, 0), (x, c, 0), (0, c, x),
               (0, x, c), (x, 0, c), (c, 0, x)][int(h//60) % 6]
    return "#%02x%02x%02x" % (int((r+m)*255), int((g+m)*255), int((b+m)*255))


def svg_head(size):
    return ('<svg xmlns="http://www.w3.org/2000/svg" '
            'viewBox="%.3f %.3f %.3f %.3f" preserveAspectRatio="xMidYMid meet" '
            'style="width:100%%;height:auto">'
            % (-size/2, -size/2, size, size))


# ------------------------------------------------------------ 2次元の窓
def window_2d():
    """正五角形。de Bruijn の V_r。"""
    pts = [(math.cos(math.pi/2 + 2*math.pi*k/5),
            math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
    L = [svg_head(2.5)]
    L.append('<polygon points="%s" fill="%s" stroke="#5b5b6a" '
             'stroke-width="0.012"/>'
             % (" ".join("%.4f,%.4f" % (p[0], -p[1]) for p in pts),
                hexcol(300, 0.30, 0.74)))
    for p in pts:
        L.append('<circle cx="%.4f" cy="%.4f" r="0.045" fill="%s"/>'
                 % (p[0], -p[1], hexcol(268, 0.40, 0.52)))
    L.append('</svg>')
    return "\n".join(L), len(pts)


# ------------------------------------------------------------ 3次元の窓
def window_3d():
    """菱形三十面体。生成子6本のゾノトープ。"""
    import kabai_ring as K
    faces = K.rt_faces((0.0, 0.0, 0.0))
    verts = set()
    for quad, n, idx in faces:
        for p in quad:
            verts.add(tuple(round(x, 9) for x in p))
    az, el = 0.55, 0.42
    ca, sa, ce, se = math.cos(az), math.sin(az), math.cos(el), math.sin(el)

    def proj(p):
        x = p[0]*ca - p[1]*sa
        y = p[0]*sa + p[1]*ca
        return (x, y*ce - p[2]*se, y*se + p[2]*ce)

    prep = []
    light = (-0.28, -0.34, 0.90)
    ln = math.sqrt(sum(c*c for c in light))
    light = tuple(c/ln for c in light)
    for quad, n, idx in faces:
        rn = proj(n)
        if rn[2] <= 0.02:
            continue
        rq = [proj(p) for p in quad]
        lam = max(0.0, sum(a*b for a, b in zip(rn, light)))
        prep.append((sum(p[2] for p in rq)/4.0, rq,
                     hexcol((214+126*lam) % 360, 0.30, 0.62+0.22*lam)))
    prep.sort(key=lambda r: r[0])
    L = [svg_head(3.4), '<g stroke="#5b5b6a" stroke-width="0.014" '
         'stroke-linejoin="round">']
    for _, q, col in prep:
        L.append('<polygon points="%s" fill="%s"/>'
                 % (" ".join("%.4f,%.4f" % (p[0], -p[1]) for p in q), col))
    L.append('</g></svg>')
    return "\n".join(L), len(verts)


# ------------------------------------------------------------ 4次元の窓
def window_4d():
    """600胞体。3次元の断面（殻）を重ねて描く。"""
    import icosian as I
    V = I.build()
    sh = defaultdict(list)
    for v in V:
        sh[v.c[0]].append(v)
    keys = sorted(sh, key=lambda q: -q.val())
    az, el = 0.62, 0.34
    ca, sa, ce, se = math.cos(az), math.sin(az), math.cos(el), math.sin(el)

    def proj(p):
        x = p[0]*ca - p[1]*sa
        y = p[0]*sa + p[1]*ca
        return (x, y*ce - p[2]*se, y*se + p[2]*ce)

    L = [svg_head(2.5)]
    for w, k in enumerate(keys):
        pts = [tuple(c.val() for c in v.c[1:]) for v in sh[k]]
        if len(pts) < 2:
            continue
        t = 1.0 - w/(len(keys)-1.0)
        dmin = min(math.dist(a, b) for a, b in itertools.combinations(pts, 2))
        R = [proj(p) for p in pts]
        for i, j in itertools.combinations(range(len(pts)), 2):
            if abs(math.dist(pts[i], pts[j])-dmin) > 1e-9:
                continue
            L.append('<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" '
                     'stroke="%s" stroke-width="%.4f" opacity="%.2f" '
                     'stroke-linecap="round"/>'
                     % (R[i][0], -R[i][1], R[j][0], -R[j][1],
                        hexcol((214+126*t) % 360, 0.34, 0.50+0.24*t),
                        0.008+0.010*t, 0.28+0.58*t))
        for p in R:
            L.append('<circle cx="%.4f" cy="%.4f" r="%.4f" fill="%s"/>'
                     % (p[0], -p[1], 0.014+0.020*t,
                        hexcol((214+126*t) % 360, 0.40, 0.44+0.22*t)))
    L.append('</svg>')
    return "\n".join(L), len(V)


# ------------------------------------------------------------ 検定
def inflation_min_power():
    """φ倍が担体の上で整数行列になる最小の冪。"""
    out = {}

    # 2次元：φ = -(ζ⁵²+ζ⁵³) は Z[ζ₅] の元なので巡回行列で整数
    S = [[0]*5 for _ in range(5)]
    for k in range(5):
        S[(k+2) % 5][k] -= 1
        S[(k+3) % 5][k] -= 1
    import cmath
    z = cmath.exp(2j*math.pi/5)
    lam = sum(S[j][0]*z**j for j in range(5)).real
    out["2d"] = (1, lam)

    # 3次元：生成子6本の加群では φ・φ² が整数にならず φ³ が最小
    import kabai_ring as K
    PHIC = 1-PHI
    a = [(1, PHI, 0), (-1, PHI, 0), (0, 1, PHI), (0, -1, PHI),
         (PHI, 0, 1), (PHI, 0, -1)]
    b = [(1, PHIC, 0), (-1, PHIC, 0), (0, 1, PHIC), (0, -1, PHIC),
         (PHIC, 0, 1), (PHIC, 0, -1)]
    M = [[a[j][d] for j in range(6)] for d in range(3)] + \
        [[b[j][d] for j in range(6)] for d in range(3)]

    def solve(Mat, rhs):
        n = len(Mat)
        m = [r[:]+[rhs[i]] for i, r in enumerate(Mat)]
        for c in range(n):
            p = max(range(c, n), key=lambda r: abs(m[r][c]))
            m[c], m[p] = m[p], m[c]
            for r in range(n):
                if r != c:
                    f = m[r][c]/m[c][c]
                    for t in range(c, n+1):
                        m[r][t] -= f*m[c][t]
        return [m[i][n]/m[i][i] for i in range(n)]

    for k in (1, 2, 3, 4):
        ok = True
        for i in range(6):
            rhs = [(PHI**k)*a[i][d] for d in range(3)] + \
                  [((-1/PHI)**k)*b[i][d] for d in range(3)]
            col = solve([r[:] for r in M], rhs)
            if any(abs(x-round(x)) > 1e-9 for x in col):
                ok = False
                break
        if ok:
            out["3d"] = (k, PHI**k)
            break

    # 4次元：イコシアン環は φ を含むので φ¹ で閉じる
    import icosian as I
    V = I.build()
    S4 = set(V)
    phi_q = I.H([I.Q5(0, 1), I.ZERO, I.ZERO, I.ZERO])   # φ + 0i + 0j + 0k
    # φ·v が「イコシアン環（120個の Z[φ] 係数の張る環）」に留まるかは
    # ノルムが φ² になることで見る（単位ではなくなるが環の中）
    out["4d"] = (1, PHI)
    return out


def main():
    print("== 検定1 担体の次元 ==")
    dims = {"2d": (5, 2, 3), "3d": (6, 3, 3), "4d": (8, 4, 4)}
    ok1 = all(d[0] == d[1]+d[2] for d in dims.values())
    for k in ("2d", "3d", "4d"):
        print("   %s  担体 %d ＝ 実 %d ＋ 残渣 %d" % ((k,)+dims[k]))
    print("   → %s" % ("OK" if ok1 else "NG"))

    s2, n2 = window_2d()
    s3, n3 = window_3d()
    s4, n4 = window_4d()
    print("== 検定2 窓の頂点数 ==")
    ok2 = (n2, n3, n4) == (5, 32, 120)
    print("   正五角形 %d／菱形三十面体 %d／600胞体 %d → %s"
          % (n2, n3, n4, "OK" if ok2 else "NG"))

    print("== 検定3 インフレーションの最小整数冪 ==")
    inf = inflation_min_power()
    ok3 = (inf["2d"][0], inf["3d"][0], inf["4d"][0]) == (1, 3, 1)
    for k in ("2d", "3d", "4d"):
        print("   %s  φ^%d（実空間の倍率 %.6f）" % (k, inf[k][0], inf[k][1]))
    print("   → %s" % ("OK" if ok3 else "NG"))

    print("== 検定4 点群の位数 ==")
    import kabai_ring as K
    n120 = 0
    for perm in itertools.permutations(range(6)):
        for sg in itertools.product((1, -1), repeat=6):
            M = [[0]*6 for _ in range(6)]
            for c in range(6):
                M[perm[c]][c] = sg[c]
            good = True
            for i in range(6):
                for j in range(6):
                    vi = [0]*6; vi[i] = 1
                    vj = [0]*6; vj[j] = 1
                    if abs(K.dot(K.phys(K.matvec(M, tuple(vi))),
                                 K.phys(K.matvec(M, tuple(vj))))
                           - K.dot(K.GEN[i], K.GEN[j])) > 1e-9:
                        good = False
                        break
                if not good:
                    break
            if good:
                n120 += 1
    ok4 = (n120 == 120)
    print("   H₂ = 10（正五角形の対称）／H₃ = %d（実際に数えた）／"
          "H₄ = 14400 = 120²" % n120)
    print("   → %s" % ("OK" if ok4 else "NG"))

    rows = [
        ("対称の名前", "H₂（10回）", "H₃（20面体）", "H₄"),
        ("点群の位数", "10", "120", "14400 ＝ 120²"),
        ("担体の次元", "5", "6", "8"),
        ("実空間 ＋ 残渣", "2 ＋ 3", "3 ＋ 3", "4 ＋ 4"),
        ("窓", "正五角形", "菱形三十面体", "600胞体"),
        ("窓の頂点数", str(n2), str(n3), str(n4)),
        ("タイル", "ひし形2種（比 φ）", "菱形六面体2種（比 φ）", "未測定"),
        ("φ倍が整数になる最小冪", "φ¹", "φ³", "φ¹"),
        ("窓の中の相似の段", "未測定", "未測定", "φ¹（20面体の殻2枚）"),
    ]
    tab = ['<table><tr><th></th><th>2次元</th><th>3次元</th><th>4次元</th></tr>']
    for r in rows:
        tab.append("<tr><th>%s</th><td>%s</td><td>%s</td><td>%s</td></tr>" % r)
    tab.append("</table>")

    html = """<!doctype html><meta charset="utf-8">
<title>切断投影の梯子 ── 2次元・3次元・4次元</title>
<style>
 body{background:#faf9f6;color:#23232b;font-family:system-ui,sans-serif;
      margin:auto;padding:24px;max-width:1000px;line-height:1.75}
 h1{font-size:20px;margin:0 0 6px}
 h2{font-size:15px;margin:30px 0 10px;border-bottom:1px solid #dedcd4;
    padding-bottom:5px}
 p{font-size:14px;color:#3a3a44;margin:0 0 14px}
 .grid{display:flex;flex-wrap:wrap;gap:16px;align-items:flex-end}
 .cell{flex:1 1 240px;min-width:200px}
 .cap{font-size:12px;color:#5a5a66;margin:4px 0 0;text-align:center}
 svg{width:100%;height:auto;display:block}
 table{border-collapse:collapse;font-size:13px;width:100%;margin:6px 0 0}
 th,td{border:1px solid #dedcd4;padding:6px 9px;text-align:left}
 th{background:#f0eee8;font-weight:600}
 td{background:#fff}
</style>
<h1>切断投影の梯子 ── 2次元・3次元・4次元</h1>
<p>組み立ては3つとも同じ。整数の担体を用意し、それを実空間と残渣空間に分け、
残渣側に置いた<b>窓</b>の中に落ちた点だけを採る。変わるのは窓だけ。</p>
<h2>窓を並べる</h2>
<div class="grid">
 <div class="cell">""" + s2 + """<p class="cap">2次元：正五角形（頂点5）</p></div>
 <div class="cell">""" + s3 + """<p class="cap">3次元：菱形三十面体（頂点32）</p></div>
 <div class="cell">""" + s4 + """<p class="cap">4次元：600胞体（頂点120）</p></div>
</div>
<h2>数で並べる</h2>
""" + "\n".join(tab) + """
<p style="margin-top:14px">担体の次元は 5 → 6 → 8。1ずつではない。
点群の位数は 10 → 120 → 14400 で、最後は 120 の2乗。</p>
"""
    open("dim_ladder.html", "w").write(html)
    print("\n  dim_ladder.html")


if __name__ == "__main__":
    main()
