#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
draw_star10.py ── 五芒星を中心に、91環の三角形を10枚組む
============================================================
内側5枚は円環中心を半径 φ³、外側5枚は φ⁵ に置く（差は φ⁴）。
円環の内側だけを 2,3 の篩の色で塗る。

事前登録（判別法6）
  検定α 10枚が重なりなく置けるか
      OK なら：中心を五芒星に取った監督の見立てが成立する
      NG なら：φ³/φ⁵ の取り方が違う
  検定β 隙間は4種（細ひし形・五芒星・舟・正十角形）だけか
      OK なら：局所配置がペンローズの4種から出ない
      NG なら：新しい隙間が出ている＝タイリングとして破綻
  検定γ 中心の隙間は五芒星か
  検定δ 図全体は ζ²(72°) で閉じ、ζ(36°) では閉じないか
      OK なら：対称群は D₅（10枚でも10回対称にはならない）
"""
import math, sys
from collections import Counter, defaultdict
import b13_chain_units as B

rows, place, offs = B.build_stack()
T = (-2, 2, 0, 3)     # |T| = φ³  画面方位90°
S = (-5, 5, 0, 8)     # |S| = φ⁵  画面方位90°

ringrow = {}
U = []; seen = set()
for k in range(10):
    t = B.zrot(T if k % 2 == 0 else S, k)
    for ri, r in enumerate(place):
        for c in r:
            c2 = B.zadd(B.zrot(c, k), t)
            if c2 not in seen:
                seen.add(c2); U.append(c2); ringrow[c2] = ri
cells = B.fits(U)

print("検定α 10枚 重なり0        %s  環%d 五角形%d 共有%d"
      % ("OK" if cells else "NG", len(U), len(cells) if cells else 0, 6280 - (len(cells) if cells else 0)))
fs = B.gaps(cells)
kinds = Counter(round(a, 4) for a, _ in fs)
print("検定β 隙間の種類          %s  %s"
      % ("OK" if set(kinds) == set(B.GAP_NAME) else "NG",
         {B.GAP_NAME.get(k, "不明%.4f" % k): v for k, v in sorted(kinds.items())}))

inner = []; ctr = None
for a, cyc in fs:
    P = [B.xy(p) for p in cyc]
    g = (sum(p[0] for p in P) / len(P), sum(p[1] for p in P) / len(P))
    c = min(ringrow, key=lambda c: math.dist(g, B.xy(c)))
    if math.dist(g, B.xy(c)) < 2.0:
        inner.append((ringrow[c], round(a, 4), cyc))
    if math.hypot(*g) < 1e-6:
        ctr = round(a, 4)
print("検定γ 中心の隙間          %s  %s" % ("OK" if ctr == 2.9389 else "NG", B.GAP_NAME.get(ctr, "?")))
Sset = set(U)
c1 = all(B.zrot(c, 1) in Sset for c in U); c2 = all(B.zrot(c, 2) in Sset for c in U)
print("検定δ ζ²で閉じ ζでは閉じず %s  (ζ=%s, ζ²=%s)" % ("OK" if (c2 and not c1) else "NG", c1, c2))

# ── 描画 ───────────────────────────────────────────────────
def poly(q, a): return [B.xy(B.zadd(q, B.zt(a + 2 * i))) for i in range(5)]
pts = [p for q, a in cells.items() for p in poly(q, a)]
x0 = min(p[0] for p in pts) - 2; x1 = max(p[0] for p in pts) + 2
y0 = min(p[1] for p in pts) - 2; y1 = max(p[1] for p in pts) + 2
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="%.2f %.2f %.2f %.2f">' % (x0, y0, x1 - x0, y1 - y0),
     '<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="#0e1014"/>' % (x0, y0, x1 - x0, y1 - y0),
     '<g stroke="#3a414c" stroke-width="0.07" fill="none">']
for q, a in cells.items():
    o.append('<polygon points="' + " ".join("%.4f,%.4f" % p for p in poly(q, a)) + '"/>')
o.append('</g><g stroke="none">')
for ri, _, cyc in inner:
    d = "M" + "L".join("%.4f %.4f" % p for p in (B.xy(p) for p in cyc)) + "Z"
    o.append('<path d="%s" fill="%s"/>' % (d, B.COLOR[B.sieve_class(B.ROW_ORDER[ri])]))
o.append('</g></svg>')
out = sys.argv[1] if len(sys.argv) > 1 else "b13_star10.svg"
open(out, "w").write("\n".join(o))
print("→", out)
