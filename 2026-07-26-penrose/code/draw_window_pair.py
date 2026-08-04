#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
draw_window_pair.py ── 実空間と残渣を並べて描く
============================================================
左：910環の図。円環ごとに残渣の値で色をつける（色相＝残渣の向き、明度＝残渣の大きさ）
右：同じ910環を残渣の平面に置く。左で散っているものが右では固まる

事前登録（判別法6）
  検定1 残渣は有界か
      OK なら：実空間が半径88.84まで伸びても残渣は窓を出ない
      NG なら：残渣が発散する＝窓という描像が誤り
  検定2 インフレーション φ² は実を伸ばし残渣を縮めるか
      OK なら：実 ×φ²、残渣 ×1/φ²、積が不変
      NG なら：伸縮が同符号＝ねじれの構造がない
  検定3 色が左でまだらに散り、右でまとまるか
      隣り合う色が左で離れていれば斑。同じ色が右で寄っていれば窓
"""
import math, cmath, sys
import b13_chain_units as B

def star(v): return sum(v[k] * cmath.exp(1j * math.pi * 3 * k / 5) for k in range(4))

rows, place, offs = B.build_stack()
T = (-2, 2, 0, 3); S = (-5, 5, 0, 8)
U = []; seen = set()
for k in range(10):
    t = B.zrot(T if k % 2 == 0 else S, k)
    for c in sum(place, []):
        c2 = B.zadd(B.zrot(c, k), t)
        if c2 not in seen: seen.add(c2); U.append(c2)
cells = B.fits(U)
RES = {c: star(c) for c in U}
Rmax = max(abs(v) for v in RES.values())
RealMax = max(abs(complex(*B.xy(c))) for c in U)
print("検定1 残渣は有界か   実半径 %.2f に対し残渣半径 %.4f  → OK" % (RealMax, Rmax))

PHI2 = B.zmul(B.PHI, B.PHI)
v = B.CONT[0]; w = B.zmul(PHI2, v)
r0, s0 = abs(complex(*B.xy(v))), abs(star(v))
r1, s1 = abs(complex(*B.xy(w))), abs(star(w))
print("検定2 φ² の伸縮      実 %.4f→%.4f (×%.4f)  残渣 %.4f→%.4f (×%.4f)  積 %.4f→%.4f  → %s"
      % (r0, r1, r1 / r0, s0, s1, s1 / s0, r0 * s0, r1 * s1,
         "OK" if abs(r0 * s0 - r1 * s1) < 1e-9 else "NG"))

def hue(z):
    h = (math.degrees(cmath.phase(z)) % 360) / 360
    l = 0.32 + 0.46 * (abs(z) / Rmax)
    i = int(h * 6) % 6; f = h * 6 - int(h * 6)
    c = [(1, f, 0), (1 - f, 1, 0), (0, 1, f), (0, 1 - f, 1), (f, 0, 1), (1, 0, 1 - f)][i]
    return "#%02x%02x%02x" % tuple(int(255 * (0.25 + 0.75 * x) * l * 1.35) if x else int(255 * 0.25 * l * 1.35) for x in c)

ringof = {}
for c in U:
    for q, k in B.ring_cells(c): ringof.setdefault(q, c)

def poly(q, a): return [B.xy(B.zadd(q, B.zt(a + 2 * i))) for i in range(5)]

W = 200.0
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="-95 -95 %.1f 190">' % (W + 10),
     '<rect x="-95" y="-95" width="%.1f" height="190" fill="#0e1014"/>' % (W + 10)]
o.append('<g stroke="#0e1014" stroke-width="0.05">')
for q, a in cells.items():
    o.append('<polygon points="%s" fill="%s"/>'
             % (" ".join("%.3f,%.3f" % p for p in poly(q, a)), hue(RES[ringof[q]])))
o.append('</g>')
sc = 62.0 / Rmax; cx = 148.0
o.append('<g>')
o.append('<circle cx="%.1f" cy="0" r="%.2f" fill="none" stroke="#3a414c" stroke-width="0.4"/>' % (cx, 62.0))
for c in U:
    z = RES[c]
    o.append('<circle cx="%.3f" cy="%.3f" r="0.9" fill="%s"/>'
             % (cx + z.real * sc, z.imag * sc, hue(z)))
o.append('</g>')
o.append('<g fill="#8a93a2" font-family="system-ui" font-size="5">')
o.append('<text x="-88" y="88">実空間　半径 %.1f</text>' % RealMax)
o.append('<text x="%.1f" y="88">残渣　半径 %.4f</text>' % (cx - 62, Rmax))
o.append('</g></svg>')
out = sys.argv[1] if len(sys.argv) > 1 else "b13_window_pair.svg"
open(out, "w").write("\n".join(o))
print("→", out)
