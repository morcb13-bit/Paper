#  図  六角とペンローズの見え方の違い／動きを追うとは何か
#
#      対象：同じ円板・同じ程度のセル数で、同じ像を標本したときの励起パターン。
#            および、像を一歩ずらしたときにパターンがどう変わるか。
#
#  検定F1  二つの網膜が比べられる状態にあるか
#      測る量：セル数と、覆っている円板の半径
#      OK なら：セル数の差が 5% 以内・同じ半径 → 以後の比較は並べ方だけの違い
#      NG なら：数か広さが違う → 比較にならない
#
#  検定F2  格子ベクトルで像をずらすと、励起パターンは自分に重なるか（台地）
#      測る量：ずらす前と後の励起セル集合の一致率（ずらした分を戻して比べる）
#      予想：六角は格子ベクトルのずれで 100%（模様がそのまま平行移動する）
#            ペンローズは 100% にならない（ずれが模様に出る）
#      OK なら：ペンローズ側でだけ、ずれが模様の変化として読める
#      NG なら：ペンローズも 100% → 動きを模様から読む筋が立たない
#
#  検定F3  五芒星の符牒はずれで変わるか
#      測る量：ずらす前と後で符牒（周囲5枚の励起の5ビット）が変わった五芒星の個数
#      必ず落ちる設定：ずれ 0 のとき、変わった五芒星が 0 個であること
#
#  検定F4  版面の重なり
#      測る量：各パネルの外接円が版面からはみ出さず、隣と重ならないこと
#      NG なら：図を出さない

exec(open('wind_core.py').read())

import math

# ── 担体（ペンローズ側） ────────────────────────────────
PXY = {q: tuple(float(t) for t in U.xy(q)) for q in CELLS}
cx = sum(p[0] for p in PXY.values()) / len(PXY)
cy = sum(p[1] for p in PXY.values()) / len(PXY)
PXY = {q: (p[0] - cx, p[1] - cy) for q, p in PXY.items()}
RAD = max(math.hypot(*p) for p in PXY.values())
NP = len(PXY)

# ── 共通の領域＝担体の凸包（円板ではない） ───────────────
def hull(P):
    P = sorted(P)
    def half(P):
        h = []
        for p in P:
            while len(h) >= 2 and ((h[-1][0]-h[-2][0])*(p[1]-h[-2][1])
                                   - (h[-1][1]-h[-2][1])*(p[0]-h[-2][0])) <= 0:
                h.pop()
            h.append(p)
        return h
    return half(P)[:-1] + half(P[::-1])[:-1]

HULL = hull(list(PXY.values()))

def inside(p):
    n = len(HULL)
    for i in range(n):
        a, b = HULL[i], HULL[(i + 1) % n]
        if (b[0]-a[0])*(p[1]-a[1]) - (b[1]-a[1])*(p[0]-a[0]) < -1e-9:
            return False
    return True

# ── 六角格子（同じ領域・同じ程度のセル数） ──────────────
def hexcells(s):
    out = []
    n = int(RAD / s) + 3
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            x = s * (i + j * 0.5)
            y = s * j * math.sqrt(3) / 2
            if inside((x, y)):
                out.append(((i, j), (x, y)))
    return out

lo, hi = 0.5, 6.0
for _ in range(60):
    mid = (lo + hi) / 2
    if len(hexcells(mid)) > NP:
        lo = mid
    else:
        hi = mid
HS = (lo + hi) / 2
HEX = dict(hexcells(HS))
NH = len(HEX)

print(f"ペンローズ  五角形 {NP} 枚 / 凸包の頂点 {len(HULL)} / 外接半径 {RAD:.2f}")
print(f"六角格子    セル {NH} 個 / 同じ凸包の中 / 間隔 {HS:.3f}")
f1 = abs(NP - NH) / NP <= 0.05
print(f"検定F1 セル数の差 {abs(NP-NH)/NP*100:.1f}%   {'OK' if f1 else 'NG'}")
if not f1:
    raise SystemExit("検定F1 NG。図を出さない")

# ── 像：数字の 3（太さを持つ折れ線） ──────────────────
def digit3(scale, ox, oy):
    pts = [(-0.45, 0.80), (0.35, 0.80), (-0.10, 0.18), (0.30, 0.18),
           (0.45, -0.05), (0.40, -0.55), (0.00, -0.80), (-0.42, -0.62)]
    return [(x * scale + ox, y * scale + oy) for x, y in pts]

def seglist(pts):
    return list(zip(pts[:-1], pts[1:]))

def near(p, segs, w):
    x, y = p
    for (x1, y1), (x2, y2) in segs:
        dx, dy = x2 - x1, y2 - y1
        L = dx * dx + dy * dy
        t = 0.0 if L == 0 else max(0.0, min(1.0, ((x - x1) * dx + (y - y1) * dy) / L))
        if math.hypot(x - (x1 + t * dx), y - (y1 + t * dy)) <= w:
            return True
    return False

SCALE = RAD * 0.46
WID = RAD * 0.045

def excite(cells, segs):
    return {k for k, p in cells.items() if near(p, segs, WID)}

OY = RAD * 0.26
segs0 = seglist(digit3(SCALE, 0, OY))
EP0 = excite(PXY, segs0)
EH0 = excite(HEX, segs0)
print(f"像「3」   ペンローズ励起 {len(EP0)} / 六角励起 {len(EH0)}")

# ── 検定F2 格子ベクトルでずらす ────────────────────────
DX, DY = HS * 1.0, 0.0                      # 六角の格子ベクトル 1 本
segs1 = seglist(digit3(SCALE, DX, OY + DY))
EH1 = excite(HEX, segs1)
EP1 = excite(PXY, segs1)

EH1b = {(i - 1, j) for (i, j) in EH1}        # ずらした分を戻す

# 縁で出入りするセルは退化として除外する（除外件数を出す）
def margin(p, m):
    n = len(HULL)
    for i in range(n):
        a, b = HULL[i], HULL[(i + 1) % n]
        L = math.hypot(b[0]-a[0], b[1]-a[1])
        d = ((b[0]-a[0])*(p[1]-a[1]) - (b[1]-a[1])*(p[0]-a[0])) / L
        if d < m:
            return False
    return True

MG = max(HS, 1.618) * 2.5
HCORE = {k for k, p in HEX.items() if margin(p, MG)}
PCORE = {q for q, p in PXY.items() if margin(p, MG)}
print(f"縁の帯 {MG:.2f} を除外：六角 {len(HEX)-len(HCORE)} / ペンローズ {len(PXY)-len(PCORE)} セル")

hex_match = len((EH1b & EH0) & HCORE) / max(1, len(EH0 & HCORE))

def shift_back(cells, ex, dx, dy):
    """ずらした励起を −(dx,dy) 平行移動し、最も近いセルへ着地させ直す"""
    out = set()
    for k in ex:
        x, y = cells[k][0] - dx, cells[k][1] - dy
        best, bd = None, None
        for k2, p2 in cells.items():
            d = (p2[0] - x) ** 2 + (p2[1] - y) ** 2
            if bd is None or d < bd:
                best, bd = k2, d
        out.add(best)
    return out

EP1b = shift_back(PXY, EP1, DX, DY)
pen_match = len((EP1b & EP0) & PCORE) / max(1, len(EP0 & PCORE))

print(f"検定F2 六角   一致率 {hex_match*100:.1f}%")
print(f"       ペンローズ 一致率 {pen_match*100:.1f}%")
f2 = hex_match > 0.99 and pen_match < 0.99
print(f"       {'OK' if f2 else 'NG'}（六角は台地・ペンローズは台地でない）")

# ── 検定F3 五芒星の符牒 ────────────────────────────────
rows2, place2, offs2 = U.build_stack()
faces = U.gaps(U.fits(sum(place2, [])))
SC = []
for area, cyc in faces:
    if abs(area - 2.9389) < 0.01:
        P = [U.xy(p) for p in cyc]
        SC.append((sum(a for a, _ in P) / len(P) - cx,
                   sum(b for _, b in P) / len(P) - cy))
STARS = SC
print(f"五芒星 {len(SC)} 個")

def around(c, k=5):
    d = sorted(((math.hypot(p[0] - c[0], p[1] - c[1]), q) for q, p in PXY.items()))
    return [q for _, q in d[:k]]

RING = [around(c) for c in SC]

def code(ex):
    return [tuple(1 if q in ex else 0 for q in r) for r in RING]

c0, c1 = code(EP0), code(EP1b)
changed = sum(1 for a, b in zip(c0, c1) if a != b)
same = sum(1 for a, b in zip(c0, code(EP0)) if a != b)
print(f"検定F3 ずれ1歩で符牒が変わった五芒星 {changed} / {len(STARS)}")
print(f"       ずれ0 で変わった五芒星 {same} 個   {'OK' if same == 0 else 'NG'}（必ず落ちる設定）")

# ── 作図 ───────────────────────────────────────────────
W, H = 980, 520
PANW = 460
CXs = [250, 730]
CYp = 285
SCL = 195 / RAD

def svg_open(w, h):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
            f'viewBox="0 0 {w} {h}" font-family="Noto Sans CJK JP, sans-serif">',
            f'<rect width="{w}" height="{h}" fill="#ffffff"/>']

def hexpath(x, y, r):
    pts = [(x + r * math.cos(math.pi / 6 + k * math.pi / 3),
            y + r * math.sin(math.pi / 6 + k * math.pi / 3)) for k in range(6)]
    return "M" + "L".join(f"{a:.2f},{b:.2f}" for a, b in pts) + "Z"

def pentpath(q, x, y, r):
    pts = [(x + r * math.cos(math.radians(-18 + 72 * k + 36 * (sum(q) % 2))),
            y + r * math.sin(math.radians(-18 + 72 * k + 36 * (sum(q) % 2)))) for k in range(5)]
    return "M" + "L".join(f"{a:.2f},{b:.2f}" for a, b in pts) + "Z"

ON, OFF, EDGE = "#1a3f6b", "#eef1f4", "#c9d2da"

def draw_hex(out, ox, oy, ex):
    r = HS * SCL * 0.56
    for k, (x, y) in HEX.items():
        X, Y = ox + x * SCL, oy - y * SCL
        out.append(f'<path d="{hexpath(X, Y, r)}" fill="{ON if k in ex else OFF}" '
                   f'stroke="{EDGE}" stroke-width="0.4"/>')

def draw_pen(out, ox, oy, ex):
    r = 1.0 * SCL * 0.94
    for q, (x, y) in PXY.items():
        X, Y = ox + x * SCL, oy - y * SCL
        out.append(f'<path d="{pentpath(q, X, Y, r)}" fill="{ON if q in ex else OFF}" '
                   f'stroke="{EDGE}" stroke-width="0.4"/>')

# 検定F4 版面
f4 = all(0 < c - 205 and c + 205 < W for c in CXs) and (CXs[1] - CXs[0]) > 410
print(f"検定F4 版面の重なり   {'OK' if f4 else 'NG'}")
if not f4:
    raise SystemExit("検定F4 NG。図を出さない")

# 図1 同じ像を二つの網膜で見る
o = svg_open(W, H)
o.append(f'<text x="{W/2}" y="34" text-anchor="middle" font-size="19" fill="#111">'
         f'同じ像を、同じ広さ・同じ程度のセル数で見る</text>')
o.append(f'<text x="{CXs[0]}" y="66" text-anchor="middle" font-size="15" fill="#444">'
         f'六角格子　{NH} セル</text>')
o.append(f'<text x="{CXs[1]}" y="66" text-anchor="middle" font-size="15" fill="#444">'
         f'ペンローズ　{NP} セル</text>')
draw_hex(o, CXs[0], CYp, EH0)
draw_pen(o, CXs[1], CYp, EP0)
o.append(f'<text x="{CXs[0]}" y="500" text-anchor="middle" font-size="13" fill="#666">'
         f'励起 {len(EH0)} セル</text>')
o.append(f'<text x="{CXs[1]}" y="500" text-anchor="middle" font-size="13" fill="#666">'
         f'励起 {len(EP0)} セル</text>')
o.append('</svg>')
open('fig1_two_retinas.svg', 'w').write("\n".join(o))

# 図2 一歩ずらす
o = svg_open(W, H)
o.append(f'<text x="{W/2}" y="34" text-anchor="middle" font-size="19" fill="#111">'
         f'像を六角格子の一歩ぶん動かして、模様を戻して重ねる</text>')
o.append(f'<text x="{CXs[0]}" y="66" text-anchor="middle" font-size="15" fill="#444">'
         f'六角格子　一致 {hex_match*100:.0f}%</text>')
o.append(f'<text x="{CXs[1]}" y="66" text-anchor="middle" font-size="15" fill="#444">'
         f'ペンローズ　一致 {pen_match*100:.0f}%</text>')

DIFF = "#d4553f"
def draw_hex_diff(out, ox, oy, a, b):
    r = HS * SCL * 0.56
    for k, (x, y) in HEX.items():
        X, Y = ox + x * SCL, oy - y * SCL
        c = OFF
        if k in a and k in b: c = ON
        elif k in a or k in b: c = DIFF
        out.append(f'<path d="{hexpath(X, Y, r)}" fill="{c}" stroke="{EDGE}" stroke-width="0.4"/>')

def draw_pen_diff(out, ox, oy, a, b):
    r = 1.0 * SCL * 0.94
    for q, (x, y) in PXY.items():
        X, Y = ox + x * SCL, oy - y * SCL
        c = OFF
        if q in a and q in b: c = ON
        elif q in a or q in b: c = DIFF
        out.append(f'<path d="{pentpath(q, X, Y, r)}" fill="{c}" stroke="{EDGE}" stroke-width="0.4"/>')

draw_hex_diff(o, CXs[0], CYp, EH0, EH1b)
draw_pen_diff(o, CXs[1], CYp, EP0, EP1b)
o.append(f'<rect x="{CXs[0]-110}" y="486" width="16" height="12" fill="{ON}"/>')
o.append(f'<text x="{CXs[0]-88}" y="497" font-size="13" fill="#666">両方で励起</text>')
o.append(f'<rect x="{CXs[0]+30}" y="486" width="16" height="12" fill="{DIFF}"/>')
o.append(f'<text x="{CXs[0]+52}" y="497" font-size="13" fill="#666">片方だけ</text>')
o.append(f'<text x="{CXs[1]}" y="497" text-anchor="middle" font-size="13" fill="#666">'
         f'符牒が変わった五芒星 {changed} / {len(STARS)}</text>')
o.append('</svg>')
open('fig2_shift.svg', 'w').write("\n".join(o))


# ── 検定F5 模様だけで、セルの間隔より細かい位置に分かれるか ──
#      像を x 方向へ 0.1 セルずつ 20 段（＝2 セル分）動かす。
#      六角では、模様を格子ベクトルで重ねてよいので、重なるものは同じ模様とみなす。
#      測る量：現れた模様の種類
#      予想：六角は 1 セル進むごとに同じ列を繰り返すので約 10 通りで頭打ち
#            ペンローズは平行移動で重ねられないので 21 通りとも異なる
#      NG なら：ペンローズ側も頭打ち → 模様が位置を持たない

def canon_hex(e):
    if not e: return frozenset()
    mi = min(i for i, _ in e); mj = min(j for _, j in e)
    return frozenset((i - mi, j - mj) for i, j in e)

def count_patterns(cells, canon, core, step, n):
    seen = []
    for k in range(n + 1):
        sg = seglist(digit3(SCALE, k * step, OY))
        seen.append(canon(frozenset(excite(cells, sg))))   # 縁は像が触れないので隠さない
    return len(set(seen)), len(seen)

nh, th = count_patterns(HEX, canon_hex, HCORE, HS / 10, 20)
npn, tp = count_patterns(PXY, lambda e: e, PCORE, HS / 10, 20)
print("検定F5 0.1セルずつ20段（2セル分）動かしたときに現れる模様の種類")
print(f"       六角（平行移動で重なるものは同じ）  {nh} / {th} 段")
print(f"       ペンローズ                        {npn} / {tp} 段")
f5 = nh < npn
print(f"       {'OK' if f5 else 'NG'}（六角は頭打ち・ペンローズは全段で異なる）")

print("→ fig1_two_retinas.svg / fig2_shift.svg")
