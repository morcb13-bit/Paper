#  図：中心をずらすと窓になる
#  対象を一行で（判別法9）：残渣の平面で、窓 W と、ずらした窓 W+σ(d) と、
#      その重なり（＝新しい窓）。および窓が閉じる条件の図形 W+(−W)。
#
#  検定G1 担体の整合  ── OK なら：重なり0・2つ飛ばし対が全て等半径。NG なら描かない
#  検定G2 図中の数値  ── OK なら：図に書く面積比が shift_window.py と一致する
#  検定G3 版面の重なり ── OK なら：二つのパネルの外接矩形が交わらない。NG なら描かない
#
#  使い方: python3 draw_shift_window.py fig5_shift_window.svg

import sys
from fractions import Fraction as F
import b13_two_tilings as T
from qphi import Qp, zmul, zconj, zsigma, zre

OUT = sys.argv[1] if len(sys.argv) > 1 else "fig5_shift_window.svg"
SIN36 = 0.5877852522924731

FIG = T.build(); V = T.verify(FIG)
G1 = bool(V["重なり0"]) and bool(V["全て等半径"])
print(f"検定G1 担体の整合: {'OK' if G1 else 'NG'}")
if not G1: sys.exit("担体が整合しないので描かない")
rings = [tuple(F(x) for x in v) for v in FIG.rings]

def zim(a): return Qp(a[1], a[2] + a[3])
def zsub(a, b): return tuple(x - y for x, y in zip(a, b))
def norm2(a): return zre(zmul(a, zconj(a)))
def res2(a):
    s = zsigma(a); return zre(zmul(s, zconj(s)))
def XY(v):
    s = zsigma(v); return (zre(s), zim(s))
def cross(o, a, b):
    return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
def hull(pts):
    pts = sorted(set(pts), key=lambda p: (p[0].val(), p[1].val()))
    lo = []
    for p in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p).sign() <= 0: lo.pop()
        lo.append(p)
    up = []
    for p in reversed(pts):
        while len(up) >= 2 and cross(up[-2], up[-1], p).sign() <= 0: up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]
def area2(poly):
    s = Qp(0, 0)
    for i in range(len(poly)):
        a, b = poly[i], poly[(i+1) % len(poly)]
        s = s + (a[0]*b[1] - b[0]*a[1])
    return s
def clip(poly, a, b):
    out = []; n = len(poly)
    for i in range(n):
        p, q = poly[i], poly[(i+1) % n]
        sp, sq = cross(a, b, p).sign(), cross(a, b, q).sign()
        if sp >= 0: out.append(p)
        if (sp > 0 and sq < 0) or (sp < 0 and sq > 0):
            d1, d2 = cross(a, b, p), cross(a, b, q)
            t = d1 / (d1 - d2)
            out.append((p[0]+(q[0]-p[0])*t, p[1]+(q[1]-p[1])*t))
    return out
def inter(A, B):
    poly = A
    for i in range(len(B)):
        if not poly: return []
        poly = clip(poly, B[i], B[(i+1) % len(B)])
    return poly
def sh(poly, s): return [(p[0]+s[0], p[1]+s[1]) for p in poly]
def pt(p): return (p[0].val(), p[1].val()*SIN36)

W = hull([XY(v) for v in rings]); AW = area2(W)

# ずらしは「クラスタ間の一歩」（実長 8.057480）
target = None
for a in rings:
    for b in rings:
        d = zsub(a, b)
        if norm2(d) == Qp(18, 29): target = d; break
    if target: break
s = XY(target)
I = inter(W, sh(W, s))
ratio = (area2(I)/AW).val()
G2 = abs(ratio - 0.5483) < 5e-4
print(f"検定G2 図中の面積比 {ratio:.4f}: {'OK' if G2 else 'NG'}")

mW = [(-p[0], -p[1]) for p in W]
DB = hull([(a[0]+b[0], a[1]+b[1]) for a in W for b in mW])

# ── 版面 ────────────────────────────────────────────────────
Wpx, Hpx = 880, 470
PA = (30, 60, 400, 380)      # x, y, w, h
PB = (470, 60, 380, 380)
G3 = PA[0]+PA[2] < PB[0]
print(f"検定G3 版面の重なり: {'OK' if G3 else 'NG'}")
if not (G2 and G3): sys.exit("検定が通らないので描かない")

def mk(panel, scale, cx=0.0, cy=0.0):
    x0, y0, w, h = panel
    ox, oy = x0 + w/2, y0 + h/2
    return lambda p: (ox + (p[0]-cx)*scale, oy - (p[1]-cy)*scale)
def poly_s(poly, m):
    return " ".join(f"{m(pt(p))[0]:.2f},{m(pt(p))[1]:.2f}" for p in poly)

mA = mk(PA, 235, cx=pt(s)[0]/2, cy=pt(s)[1]/2)
mB = mk(PB, 185)

L = []
L.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {Wpx} {Hpx}" '
         f'width="{Wpx}" height="{Hpx}" font-family="Noto Sans CJK JP, Hiragino Sans, '
         f'Yu Gothic, Meiryo, sans-serif">')
L.append(f'<rect width="{Wpx}" height="{Hpx}" fill="#ffffff"/>')
L.append('<text x="30" y="30" font-size="17" fill="#111">'
         '中心をずらすと窓が細くなる ── 残渣の平面で見たところ</text>')
L.append('<text x="30" y="50" font-size="12" fill="#555">'
         'ずらし＝クラスタ間の一歩（実長 8.057480／残渣 0.277515）</text>')

# パネルA
L.append(f'<rect x="{PA[0]}" y="{PA[1]}" width="{PA[2]}" height="{PA[3]}" '
         f'fill="none" stroke="#eeeeee"/>')
L.append(f'<polygon points="{poly_s(W, mA)}" fill="#2f6fa8" fill-opacity="0.10" '
         f'stroke="#2f6fa8" stroke-width="1.6"/>')
L.append(f'<polygon points="{poly_s(sh(W, s), mA)}" fill="#c8632f" fill-opacity="0.10" '
         f'stroke="#c8632f" stroke-width="1.6" stroke-dasharray="5 3"/>')
L.append(f'<polygon points="{poly_s(I, mA)}" fill="#4f9d69" fill-opacity="0.45" '
         f'stroke="#3a7a51" stroke-width="2"/>')
ax, ay = mA(pt((Qp(0,0), Qp(0,0))))
bx, by = mA(pt(s))
L.append(f'<line x1="{ax:.1f}" y1="{ay:.1f}" x2="{bx:.1f}" y2="{by:.1f}" '
         f'stroke="#333" stroke-width="1.2"/>')
for X, Y, t in ((ax, ay, ''), (bx, by, '')):
    L.append(f'<circle cx="{X:.1f}" cy="{Y:.1f}" r="3" fill="#333"/>')
L.append(f'<text x="{(ax+bx)/2:.0f}" y="{(ay+by)/2-8:.0f}" font-size="12" '
         f'text-anchor="middle" fill="#333">σ(d)</text>')
L.append(f'<text x="{PA[0]+10}" y="{PA[1]+PA[3]-46}" font-size="12.5" fill="#2f6fa8">'
         f'■ 窓 W（正五角形）</text>')
L.append(f'<text x="{PA[0]+10}" y="{PA[1]+PA[3]-28}" font-size="12.5" fill="#c8632f">'
         f'■ ずらした窓 W+σ(d)</text>')
L.append(f'<text x="{PA[0]+10}" y="{PA[1]+PA[3]-10}" font-size="12.5" fill="#3a7a51">'
         f'■ 新しい窓（重なり）面積比 {ratio:.4f}</text>')

# パネルB
L.append(f'<rect x="{PB[0]}" y="{PB[1]}" width="{PB[2]}" height="{PB[3]}" '
         f'fill="none" stroke="#eeeeee"/>')
L.append(f'<polygon points="{poly_s(DB, mB)}" fill="#7a4fa3" fill-opacity="0.08" '
         f'stroke="#7a4fa3" stroke-width="1.6"/>')
L.append(f'<polygon points="{poly_s(W, mB)}" fill="none" stroke="#2f6fa8" '
         f'stroke-width="1.2" stroke-dasharray="4 3"/>')
# ずらしの候補は担体の差ではなく格子 L そのものから取る。
# 担体の差だけを打つと σ(d)=σ(a)−σ(b) が構造上必ず W+(−W) の中に入り、
# どんな入力でも「全部開く」と出る（誤り86 の型）。
# L = { d : d0+d2 ≡ d1+d3 (mod 5) }（担体の差から求めた合同条件）
nin = nout = 0
B = 6
for a0 in range(-B, B+1):
 for a1 in range(-B, B+1):
  for a2 in range(-B, B+1):
   for a3 in range(-B, B+1):
    if (a0+a2-a1-a3) % 5: continue
    if a0 == a1 == a2 == a3 == 0: continue
    d = (F(a0), F(a1), F(a2), F(a3))
    q = XY(d)
    if max(abs(q[0].val()), abs(q[1].val()*SIN36)) > 1.05: continue
    px, py = mB(pt(q))
    ins = all(cross(DB[i], DB[(i+1) % len(DB)], q).sign() >= 0
              for i in range(len(DB)))
    nin += ins; nout += (not ins)
    L.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="1.7" '
             f'fill="{"#3a7a51" if ins else "#cfcfcf"}"/>')
L.append(f'<text x="{PB[0]+10}" y="{PB[1]+PB[3]-28}" font-size="12.5" fill="#7a4fa3">'
         f'■ W+(−W)（正十角形）── この中なら窓が開く</text>')
L.append(f'<text x="{PB[0]+10}" y="{PB[1]+PB[3]-10}" font-size="12.5" fill="#555">'
         f'点＝ずらしの残渣 σ(d)　開く {nin} ／ 閉じる {nout}</text>')
L.append('</svg>')

open(OUT, "w").write("\n".join(L))
print(f"→ {OUT}   （σ(d) の点 開く {nin} ／ 閉じる {nout}）")
