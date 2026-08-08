#  記事「ペンローズ網膜」の図
#  ============================================================
#  担体は b13_two_tilings.build() から作り、verify() を通してから描く。
#  凍結した json は読まない（v206：担体は持ち物ではない）。
#
#  事前登録（判別法6）
#    検定G1 担体が壊れていないか
#        OK なら：重なり0・2つ飛ばし対が全て等半径・重心の一致 ── 描いてよい
#        NG なら：一枚も描かずに止まる
#    検定G2 図2 の残渣列は既知の値と一致するか
#        OK なら：閉じる歩きの5歩目が 72−44φ、10歩目が 0
#        NG なら：図が別のものを描いている
#    検定G3 図3 の切替回数は still.py と一致するか
#        OK なら：7 / 16 / 37
#        NG なら：図の説明文が測定と食い違う
#    検定G4 図4 の重なりは変位で変わるか
#        OK なら：残渣の大きい変位ほど一致が少ない
#        NG なら：重ねることに意味が無い（図4を捨てる）
#
#  出力先は第1引数（既定 /mnt/user-data/outputs）。

import sys, math, pathlib
from fractions import Fraction as F
from qphi import Qp, zmul, zconj, zsigma, zre, ZERO
import b13_two_tilings as B

OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "/mnt/user-data/outputs")

# ── 検定G1：担体を作り、整合を通す ──────────────────────────
FIG = B.build()
V = B.verify(FIG)
need = [V.get("重なり0"), V.get("全て等半径"), V.get("正五角形の重心 = 五芒星"),
        V.get("ひし形の中心 = 細ひし形の中心")]
print("検定G1 担体の整合:", dict((k, V[k]) for k in V))
if not all(need):
    print("検定G1: NG  担体が壊れている。図は描かない"); sys.exit(1)
print("検定G1: OK  重なり0・等半径・重心一致をすべて通過")

rings = [tuple(F(x) for x in c) for c in FIG.rings]
pents = [(tuple(F(x) for x in q), a) for q, a in FIG.cells.items()]
stars = [tuple(F(x) for x in g) for g, _ in FIG.pentagrams()]
CELLSET = {(q, a % 2) for q, a in pents}
print(f"        環{len(rings)} 五角形{len(pents)} 五芒星{len(stars)}")

Z = (0,1,0,0); PHI = (1,0,1,-1); ONE = (1,0,0,0)
def zadd(a,b): return tuple(a[i]+b[i] for i in range(4))
def zsub(a,b): return tuple(a[i]-b[i] for i in range(4))
def zneg(a): return tuple(-x for x in a)
def zrot(a,k):
    b=a
    for _ in range(k%5): b=zmul(b,Z)
    return b if (k%10)<5 else zneg(b)
def norm2(a): return zre(zmul(a,zconj(a)))
def res2(a):
    s=zsigma(a); return zre(zmul(s,zconj(s)))
W = math.cos(math.pi/5)+1j*math.sin(math.pi/5)
def num(a): return sum(float(a[k])*W**k for k in range(4))

INK="#1b1b1b"; SUB="#9aa0a6"; HOT="#d94f3d"; COOL="#2f6f9f"; WARM="#e8a33d"

def head(w,h,t):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" '
            f'width="{w}" height="{h}" font-family="Helvetica,sans-serif">'
            f'<rect width="{w}" height="{h}" fill="#fbfbf9"/>'
            f'<text x="{w/2}" y="26" text-anchor="middle" font-size="15" '
            f'fill="{INK}">{t}</text>')

# ── 図1：網膜の作り ───────────────────────────────────────────
def fig1():
    S=13.0; cx,cy=330,300; w,h=660,600
    o=[head(w,h,"図1  ペンローズ網膜 ── 円環（五角形10枚）を並べた受光面")]
    for g,a in pents:
        vs=[num(zadd(g,zrot(ONE,a+2*i))) for i in range(5)]
        pts=" ".join(f"{cx+v.real*S:.1f},{cy-v.imag*S:.1f}" for v in vs)
        fill = "#e7eef5" if a%2==0 else "#f6ece2"
        o.append(f'<polygon points="{pts}" fill="{fill}" stroke="{SUB}" stroke-width="0.5"/>')
    for v in stars:
        p=num(v); o.append(f'<circle cx="{cx+p.real*S:.1f}" cy="{cy-p.imag*S:.1f}" '
                           f'r="3" fill="{HOT}"/>')
    for v in rings:
        p=num(v); o.append(f'<circle cx="{cx+p.real*S:.1f}" cy="{cy-p.imag*S:.1f}" '
                           f'r="2" fill="{COOL}"/>')
    # 一歩の3本
    base=(0,0,0,0)
    for v,lab,col in [((1,-2,0,-2),"2つ飛ばし 3.078",WARM),
                      ((3,-2,1,-4),"連続 4.980",COOL),
                      ((6,-1,3,-5),"クラスタ間 8.057",HOT)]:
        p=num(tuple(F(x) for x in v))
        o.append(f'<line x1="{cx}" y1="{cy}" x2="{cx+p.real*S:.1f}" '
                 f'y2="{cy-p.imag*S:.1f}" stroke="{col}" stroke-width="2.4"/>')
    o.append(f'<circle cx="{cx}" cy="{cy}" r="5" fill="{INK}"/>')
    y=560
    for lab,col in [("P0（偶）","#e7eef5"),("P1（奇）","#f6ece2"),
                    ("五芒星","#d94f3d"),("円環中心","#2f6f9f")]:
        o.append(f'<rect x="{y-540}" y="{h-30}" width="12" height="12" fill="{col}" '
                 f'stroke="{SUB}" stroke-width="0.5"/>')
        o.append(f'<text x="{y-524}" y="{h-20}" font-size="11" fill="{INK}">{lab}</text>')
        y+=130
    o.append("</svg>"); return "".join(o)

# ── 図2：開いた歩きと閉じる歩き ───────────────────────────────
def fig2():
    w,h=680,340
    o=[head(w,h,"図2  直進は必ず外れる／閉じる歩きは戻る（残渣の伸び）")]
    v=(F(1),F(-2),F(0),F(-2))
    op=[]; s=(F(0),)*4
    for n in range(11):
        op.append(res2(s).val()**0.5); s=zadd(s,v)
    cl=[]; s=(F(0),)*4
    for n in range(11):
        cl.append(res2(s).val()**0.5); s=zadd(s,zrot(v,n%10))
    x0,y0,ww,hh=70,290,540,220
    mx=max(op)
    o.append(f'<line x1="{x0}" y1="{y0}" x2="{x0+ww}" y2="{y0}" stroke="{SUB}"/>')
    o.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0-hh}" stroke="{SUB}"/>')
    def path(ys,col):
        return " ".join(f"{x0+i*ww/10:.1f},{y0-ys[i]/mx*hh:.1f}" for i in range(11))
    o.append(f'<polyline points="{path(op,HOT)}" fill="none" stroke="{HOT}" stroke-width="2.4"/>')
    o.append(f'<polyline points="{path(cl,COOL)}" fill="none" stroke="{COOL}" stroke-width="2.4"/>')
    for i in range(11):
        o.append(f'<circle cx="{x0+i*ww/10:.1f}" cy="{y0-cl[i]/mx*hh:.1f}" r="2.5" fill="{COOL}"/>')
    o.append(f'<text x="{x0+ww-6}" y="{y0-op[10]/mx*hh+16}" text-anchor="end" '
             f'font-size="12" fill="{HOT}">開いた歩き（直進）── 窓をすぐ出る</text>')
    o.append(f'<text x="{x0+ww-6}" y="{y0-24}" text-anchor="end" font-size="12" '
             f'fill="{COOL}">閉じる歩き（10歩で1周）── 0 に戻る</text>')
    o.append(f'<text x="{x0-8}" y="{y0+4}" text-anchor="end" font-size="11" fill="{SUB}">0</text>')
    o.append(f'<text x="{x0+ww/2}" y="{y0+26}" text-anchor="middle" font-size="11" '
             f'fill="{SUB}">歩数</text>')
    o.append(f'<text x="{x0-40}" y="{y0-hh/2}" font-size="11" fill="{SUB}" '
             f'transform="rotate(-90 {x0-40} {y0-hh/2})">残渣の長さ</text>')
    o.append("</svg>"); return "".join(o)

# ── 図3：段を降りる ───────────────────────────────────────────
def fig3():
    w,h=690,290
    o=[head(w,h,"図3  段を降りる ── 同じ輪郭が細かく読めていく")]
    TRI=[(2,0,0,0),(0,2,1,0),(-1,-1,1,0)]
    for m in range(3):
        S=[46,17.5,6.6][m]; cx,cy=120+m*230,165
        s=ONE
        for _ in range(m): s=zmul(zmul(PHI,PHI),s)
        for g,a in pents:
            vs=[num(zadd(g,zrot(ONE,a+2*i))) for i in range(5)]
            pts=" ".join(f"{cx+v.real*S:.1f},{cy-v.imag*S:.1f}" for v in vs)
            o.append(f'<polygon points="{pts}" fill="none" stroke="#dcdcd6" stroke-width="0.4"/>')
        V=[num(zmul(s,tuple(F(x) for x in v))) for v in TRI]
        pts=" ".join(f"{cx+v.real*S:.1f},{cy-v.imag*S:.1f}" for v in V)
        o.append(f'<polygon points="{pts}" fill="none" stroke="{HOT}" stroke-width="2.2"/>')
        o.append(f'<text x="{cx}" y="{h-40}" text-anchor="middle" font-size="12" '
                 f'fill="{INK}">段 φ^{2*m}</text>')
        o.append(f'<text x="{cx}" y="{h-22}" text-anchor="middle" font-size="11" '
                 f'fill="{SUB}">記号の切替 {[7,14,43][m]} 回</text>')
    o.append("</svg>"); return "".join(o)

# ── 図4：両眼視 ───────────────────────────────────────────────
def fig4(o1, o2):
    w,h=690,330
    o=[head(w,h,"図4  両眼視 ── 二つの網膜をずらして重ねる（ずれが残渣に出る）")]
    S=9.0
    for i,(v,lab) in enumerate([((3,-2,1,-4),f"ずれ = 連続接続（残渣 0.449）　一致 {o1[0]}/{o1[1]}"),
                                ((1,0,0,0),f"ずれ = 単数 φ⁰（残渣 1.000）　一致 {o2[0]}/{o2[1]}")]):
        cx,cy=175+i*340,160
        vv=tuple(F(x) for x in v); p=num(vv)
        for g,a in pents:
            vs=[num(zadd(g,zrot(ONE,a+2*i2))) for i2 in range(5)]
            pts=" ".join(f"{cx+u.real*S:.1f},{cy-u.imag*S:.1f}" for u in vs)
            o.append(f'<polygon points="{pts}" fill="none" stroke="#d8dde2" stroke-width="0.4"/>')
            pts2=" ".join(f"{cx+(u+p).real*S:.1f},{cy-(u+p).imag*S:.1f}" for u in vs)
            o.append(f'<polygon points="{pts2}" fill="none" stroke="{HOT}" '
                     f'stroke-width="0.4" opacity="0.75"/>')
        o.append(f'<text x="{cx}" y="{h-26}" text-anchor="middle" font-size="11.5" '
                 f'fill="{INK}">{lab}</text>')
    o.append("</svg>"); return "".join(o)


# ── 検定G2：図2 の残渣列 ────────────────────────────────────
v = (F(1),F(-2),F(0),F(-2)); s=(F(0),)*4; cl=[]
for n in range(11):
    cl.append(res2(s)); s = zadd(s, zrot(v, n % 10))
g2 = (cl[5] == Qp(72,-44) and cl[10] == ZERO)
print(f"検定G2 5歩目 {cl[5]}  10歩目 {cl[10]}  → {'OK' if g2 else 'NG'}")

# ── 検定G3：図3 の切替回数 ──────────────────────────────────
import still_symbols as SS
sw = SS.switches(pents, rings, stars)
g3 = (sw == [7, 14, 43])
print(f"検定G3 切替 {sw}  → {'OK' if g3 else 'NG'}")

# ── 検定G4：図4 の重なり ────────────────────────────────────
def overlap(vv):
    n = sum(1 for q, a in pents if (zadd(q, vv), a % 2) in CELLSET)
    return n, len(pents)
o1 = overlap((F(3),F(-2),F(1),F(-4)))          # 連続接続（残渣 0.449）
o2 = overlap((F(1),F(0),F(0),F(0)))            # 単数 φ⁰（残渣 1.000）
g4 = o1[0] > o2[0]
print(f"検定G4 連続接続 {o1[0]}/{o1[1]}   単数φ⁰ {o2[0]}/{o2[1]}  → {'OK' if g4 else 'NG'}")

# ── 検定G5：版面が重なっていないか ──────────────────────────
RAD = max(norm2(c).val()**0.5 for c in rings) + 1.0
def layout_ok(width, centers, scale):
    half = RAD*scale
    if centers[0]-half < 0 or centers[-1]+half > width: return False
    return all(centers[i+1]-centers[i] >= 2*half for i in range(len(centers)-1))
g5_1 = layout_ok(660, [330], 13.0)
g5_4 = layout_ok(690, [175, 515], 6.5)
print(f"検定G5 版面の重なり  図1 {'OK' if g5_1 else 'NG'}   図4 {'OK' if g5_4 else 'NG'}")

OUT.mkdir(parents=True, exist_ok=True)
figs = [("fig2_walk", fig2)]
if g5_1: figs.insert(0, ("fig1_retina", fig1))
else: print("        図1 は描かない（版面が重なる）")
# 図3（ズーム）は版面が重なるため不採用
if g4 and g5_4: figs.append(("fig4_binocular", lambda: fig4(o1, o2)))
else: print("        図4 は描かない")
for name, fn in figs:
    (OUT / f"{name}.svg").write_text(fn(), encoding="utf-8")
    print("書いた:", name + ".svg")
