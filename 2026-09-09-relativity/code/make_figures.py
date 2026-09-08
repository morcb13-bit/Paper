# 動く図で見るB13講座（その6）の静止図を作る。
#   図1 fig1_two_circle_families.svg  二族の同心円と交点、和が一定の輪／差が一定の枝
#   図2 fig2_retarded_direction.svg   外へ行くほど古い向きが届くこと
# 標準ライブラリのみ。出力は SVG。
#   図1 の交点は meet() の実計算。曲線は焦点と和／差から引いた解析形で、両者が重なることが確認になる。
#   図2 の巻き角 delay は説明のために選んだ値であって、測定値ではない。

import math

COL_A = "#5ad6ff"; COL_B = "#ffab4a"; WHITE="#ffffff"
BG="#04070c"; INK="#dfe9f0"; DIM="#7e93a3"

def meet(cA,rA,cB,rB):
    dx,dy=cB[0]-cA[0],cB[1]-cA[1]; d=math.hypot(dx,dy)
    if d<1e-9 or d>rA+rB or d<abs(rA-rB): return None
    a=(rA*rA-rB*rB+d*d)/(2*d); h2=rA*rA-a*a
    if h2<0: return None
    h=math.sqrt(h2); mx=cA[0]+a*dx/d; my=cA[1]+a*dy/d
    ux,uy=-dy/d,dx/d
    return [(mx+h*ux,my+h*uy),(mx-h*ux,my-h*uy)]

# ---------------- 図1 : 二つの同心円族と交点 ----------------
def fig1(path):
    W,H=900,520; S=150; cx,cy=W/2,H/2
    d=1.0; step=0.36; NMAX=11
    A=(-d,0.0); B=(d,0.0)
    def px(p): return (cx+p[0]*S, cy-p[1]*S)
    o=[]
    o.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')
    o.append(f'<rect width="{W}" height="{H}" fill="{BG}"/>')
    o.append(f'<clipPath id="c"><rect x="0" y="0" width="{W}" height="{H}"/></clipPath><g clip-path="url(#c)">')
    # 同心円
    for n in range(1,NMAX+1):
        r=n*step
        for c,col in ((A,COL_A),(B,COL_B)):
            p=px(c)
            o.append(f'<circle cx="{p[0]:.1f}" cy="{p[1]:.1f}" r="{r*S:.1f}" fill="none" stroke="{col}" stroke-opacity="0.42" stroke-width="1"/>')
    # 和が一定（閉じた輪）
    for m in range(6,15):
        s=m*step
        if s<=2*d+1e-9: continue
        a=s/2; b=math.sqrt(a*a-d*d)
        pts=[px((a*math.cos(th),b*math.sin(th))) for th in [i*2*math.pi/240 for i in range(241)]]
        dd=" ".join(f"{'M' if i==0 else 'L'}{x:.1f},{y:.1f}" for i,(x,y) in enumerate(pts))
        o.append(f'<path d="{dd}" fill="none" stroke="#8de6b8" stroke-opacity="0.55" stroke-width="1.4"/>')
    # 差が一定（開いた枝）
    for k in range(0,6):
        a=k*step/2
        if a>=d: continue
        col = WHITE if k==0 else COL_A
        w = 1.8 if k==0 else 1.4
        if k==0:
            o.append(f'<line x1="{cx}" y1="0" x2="{cx}" y2="{H}" stroke="{WHITE}" stroke-opacity="0.75" stroke-width="1.8"/>')
            continue
        b=math.sqrt(d*d-a*a)
        for sgn,col in ((1,COL_A),(-1,COL_B)):
            pts=[]
            for i in range(121):
                u=-2.3+4.6*i/120
                pts.append(px((sgn*a*math.cosh(u), b*math.sinh(u))))
            dd=" ".join(f"{'M' if i==0 else 'L'}{x:.1f},{y:.1f}" for i,(x,y) in enumerate(pts))
            o.append(f'<path d="{dd}" fill="none" stroke="{col}" stroke-opacity="0.8" stroke-width="{w}"/>')
    # 交点
    for i in range(1,NMAX+1):
        for j in range(1,NMAX+1):
            r=meet(A,i*step,B,j*step)
            if not r: continue
            for q in r:
                x,y=px(q)
                if -20<x<W+20 and -20<y<H+20:
                    o.append(f'<rect x="{x-2:.1f}" y="{y-2:.1f}" width="4" height="4" fill="{INK}" fill-opacity="0.95"/>')
    for c,col in ((A,COL_A),(B,COL_B)):
        p=px(c)
        o.append(f'<circle cx="{p[0]:.1f}" cy="{p[1]:.1f}" r="6" fill="{col}"/>')
    o.append('</g>')
    f='font-family="sans-serif"'
    o.append(f'<g {f} font-size="15" fill="{DIM}">')
    o.append(f'<rect x="14" y="14" width="250" height="92" fill="#04070c" fill-opacity="0.75"/>')
    o.append(f'<line x1="26" y1="38" x2="52" y2="38" stroke="#8de6b8" stroke-width="2"/><text x="62" y="43">距離の和が一定</text>')
    o.append(f'<line x1="26" y1="64" x2="52" y2="64" stroke="{COL_A}" stroke-width="2"/><text x="62" y="69">距離の差が一定</text>')
    o.append(f'<line x1="26" y1="90" x2="52" y2="90" stroke="{WHITE}" stroke-width="2"/><text x="62" y="95">差がゼロ</text>')
    o.append('</g></svg>')
    open(path,"w",encoding="utf-8").write("\n".join(o))

# ---------------- 図2 : 外へ行くほど古い向きが届く ----------------
def fig2(path):
    W,H=900,560; cx,cy=W/2,H/2+10
    o=[]
    o.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')
    o.append(f'<rect width="{W}" height="{H}" fill="{BG}"/>')
    Rorb=44.0
    rings=[92,150,208,266,324,382]
    delay=0.62   # 一段外へ出るごとに戻る角（rad）
    tips=[[],[]]
    for idx,r in enumerate(rings):
        th=-delay*(idx+1)
        o.append(f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="none" stroke="{DIM}" stroke-opacity="0.28" stroke-width="1"/>')
        # その半径に届いている「向き」
        for s,(col) in enumerate((COL_A,COL_B)):
            sg = 1 if s==0 else -1
            x=cx+sg*r*math.cos(th); y=cy-sg*r*math.sin(th)
            o.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="{col}" fill-opacity="0.9"/>')
            tips[s].append((x,y))
        x1=cx+r*math.cos(th); y1=cy-r*math.sin(th)
        x2=cx-r*math.cos(th); y2=cy+r*math.sin(th)
        o.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{INK}" stroke-opacity="0.35" stroke-width="1.4"/>')
    # 巻いていく様子
    for s,col in enumerate((COL_A,COL_B)):
        pts=[]
        sg = 1 if s==0 else -1
        for i in range(0,121):
            u=i/120
            r=Rorb+(rings[-1]-Rorb)*u
            th=-delay*((r-Rorb)/(rings[1]-rings[0]))
            pts.append((cx+sg*r*math.cos(th), cy-sg*r*math.sin(th)))
        dd=" ".join(f"{'M' if i==0 else 'L'}{x:.1f},{y:.1f}" for i,(x,y) in enumerate(pts))
        o.append(f'<path d="{dd}" fill="none" stroke="{col}" stroke-opacity="0.55" stroke-width="1.6" stroke-dasharray="5 5"/>')
    # 中心（いまの向き）
    o.append(f'<circle cx="{cx}" cy="{cy}" r="{Rorb}" fill="none" stroke="{INK}" stroke-opacity="0.35" stroke-width="1"/>')
    o.append(f'<line x1="{cx-Rorb}" y1="{cy}" x2="{cx+Rorb}" y2="{cy}" stroke="{INK}" stroke-opacity="0.5" stroke-width="1.6"/>')
    o.append(f'<circle cx="{cx+Rorb}" cy="{cy}" r="11" fill="{COL_A}"/>')
    o.append(f'<circle cx="{cx-Rorb}" cy="{cy}" r="11" fill="{COL_B}"/>')
    f='font-family="sans-serif"'
    o.append(f'<g {f} font-size="16" fill="{INK}">')
    o.append(f'<text x="{cx+Rorb+18}" y="{cy-14}">いまの向き</text>')
    o.append(f'<text x="26" y="40" fill="{DIM}">外へ行くほど、届いているのは古い向き</text>')
    o.append(f'<text x="26" y="66" fill="{DIM}">向きは同じでも、場所が捻れる</text>')
    o.append(f'<text x="{cx-rings[-1]-8}" y="{cy+rings[-1]-6}" fill="{DIM}" font-size="14">一段外＝一つ前の時刻に出た段</text>')
    o.append('</g></svg>')
    open(path,"w",encoding="utf-8").write("\n".join(o))

fig1("/mnt/user-data/outputs/fig1_two_circle_families.svg")
fig2("/mnt/user-data/outputs/fig2_retarded_direction.svg")
print("ok")
