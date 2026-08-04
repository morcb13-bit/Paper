#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
draw_fivefold.py ── 91環の三角形を頂点まわりに5回回転対称で並べ、
                    隙間に合同な図形が収まるかを描く。

事前登録（判別法6）
  検査A 5枚（72°刻み・頂点共有）が置けるか
      OK なら：頂点を共有したまま5回回転対称の図形が成立する
      NG なら：衝突している五角形の位置と個数を出す（どこが止めているか）
  検査B 隙間に合同な図形（36°回転）が収まるか
      OK なら：10枚で360°が閉じる
      NG なら：衝突の半径分布を出す。1か所なら局所、辺全体なら系統的
"""
import math, sys
import b13_chain_units as B

rows, place, offs = B.build_stack()
allc = sum(place, [])

def cellmap(cs):
    d = {}
    for c in cs:
        for q, k in B.ring_cells(c):
            d.setdefault(q, k)
    return d

def conflicts(A, C):
    out = []
    for q in A:
        pa = B.xy(q)
        for q2 in C:
            if q2 == q: continue
            if math.dist(pa, B.xy(q2)) < 1.1755 - 1e-9:
                out.append((q, q2))
    return out

def poly(q, a):
    return [B.xy(B.zadd(q, B.zt(a + 2 * i))) for i in range(5)]

# ── 5枚（72°刻み）＋ 隙間に1枚（36°回転） ──────────────────
copies = [(k, cellmap([B.zrot(c, k) for c in allc])) for k in range(0, 10, 2)]
ghost  = cellmap([B.zrot(c, 1) for c in allc])

bad = []
for i in range(5):
    for j in range(i + 1, 5):
        bad += conflicts(copies[i][1], copies[j][1])
bad_ghost = conflicts(copies[0][1], ghost) + conflicts(copies[1][1], ghost)

print("検査A 5枚（72°刻み・頂点共有）… 衝突 %d 対" % len(bad))
for q, q2 in bad[:4]:
    print("   半径 %.3f / %.3f  距離 %.4f" %
          (math.hypot(*B.xy(q)), math.hypot(*B.xy(q2)), math.dist(B.xy(q), B.xy(q2))))
print("検査B 隙間に合同な1枚（36°回転）… 衝突 %d 対" % len(bad_ghost))
rr = sorted(math.hypot(*B.xy(q)) for q, _ in bad_ghost)
print("   衝突の半径:", " ".join("%.1f" % r for r in rr))

# ── 描画 ───────────────────────────────────────────────────
PC = ["#7c8798", "#6f7b8c", "#7c8798", "#6f7b8c", "#7c8798"]
pts = [p for _, cm in copies for q, a in cm.items() for p in poly(q, a)]
pts += [p for q, a in ghost.items() for p in poly(q, a)]
x0 = min(p[0] for p in pts) - 3; x1 = max(p[0] for p in pts) + 3
y0 = min(p[1] for p in pts) - 3; y1 = max(p[1] for p in pts) + 3

o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="%.2f %.2f %.2f %.2f">'
     % (x0, y0, x1 - x0, y1 - y0),
     '<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="#0e1014"/>'
     % (x0, y0, x1 - x0, y1 - y0)]

# 隙間に置いた合同な図形（成立しないので破線の輪郭だけ）
o.append('<g fill="none" stroke="#e0603a" stroke-width="0.16" stroke-dasharray="0.9 0.7" opacity="0.85">')
for q, a in ghost.items():
    o.append('<polygon points="' + " ".join("%.4f,%.4f" % p for p in poly(q, a)) + '"/>')
o.append('</g>')

# 5枚
for idx, (k, cm) in enumerate(copies):
    o.append('<g fill="none" stroke="%s" stroke-width="0.13">' % PC[idx])
    for q, a in cm.items():
        o.append('<polygon points="' + " ".join("%.4f,%.4f" % p for p in poly(q, a)) + '"/>')
    o.append('</g>')

# 衝突点
o.append('<g fill="#f2b544" opacity="0.95">')
for q, q2 in bad:
    p = B.xy(q); o.append('<circle cx="%.4f" cy="%.4f" r="1.5"/>' % p)
o.append('</g><g fill="#e0603a" opacity="0.95">')
for q, q2 in bad_ghost:
    p = B.xy(q); o.append('<circle cx="%.4f" cy="%.4f" r="1.5"/>' % p)
o.append('</g>')
o.append('</svg>')
open(sys.argv[1] if len(sys.argv) > 1 else "b13_fivefold.svg", "w").write("\n".join(o))
print("→", sys.argv[1] if len(sys.argv) > 1 else "b13_fivefold.svg")
