#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
b13_chain_units.py ── 数の鎖を作り、積み、円環の内側を篩う
============================================================
標準ライブラリのみ。他のファイルを一切 import しない（判別法4）。
Z[ζ₁₀] の整数4成分で持ち、実数に落とすのは描画と面積表示だけ。

  python3 b13_chain_units.py            # 検定を走らせる
  python3 b13_chain_units.py out.svg    # 篩の図も書き出す

事前登録（判別法6）
-------------------
  検定1 単位の枚数
      OK なら：単位 n の五角形は 8n+2 枚で、奇偶によらない
      NG なら：鎖の作り方（奇数=連続のみ / 偶数=2つ飛ばしの反復）が違う
  検定2 円環中心の折れ角
      OK なら：奇数は全箇所108°、偶数は全箇所180°
      NG なら：接続ベクトルの選び方が指定と違う
  検定3 積んだ図の整合
      OK なら：91環・五角形628枚・重なり0で三角形が閉じる
      NG なら：落差か横ずらしの取り方が違う
  検定4 外形の共線性
      OK なら：左端・右端がそれぞれ厳密に一直線（Z[φ] 上の等式で判定）
      NG なら：三角形の辺が直線でない＝積み方が違う
  検定5 隙間の分類
      OK なら：隙間は4種（細ひし形・五芒星・舟・正十角形）だけで142個
      NG なら：五角形の重なりか、隙間の切り出しに誤りがある
  検定6 円環の内側の個数
      OK なら：奇数行 n 個・偶数行 3n/2 個、合計112個
      NG なら：内側の判定（重心から最近接円環中心まで2未満）が効いていない
  検定7 五芒星の位置
      OK なら：五芒星30個はすべて円環の外側（行と行のあいだ）
      NG なら：内側の判定が五芒星を拾っている

  NG が一件も出ない検査は反証不能である可能性を疑うこと。
  緩めて落ちるかを見る例：ROW_ORDER を逆順にすると検定3が落ちる。
"""
import sys, math, cmath, json
from collections import Counter, defaultdict

# ── Z[ζ₁₀]  ζ = e^{iπ/5}、ζ⁴ = -1+ζ-ζ²+ζ³ ────────────────────
POW = [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(-1,1,-1,1),
       (-1,0,0,0),(0,-1,0,0),(0,0,-1,0),(0,0,0,-1),(1,-1,1,-1)]
ZERO = (0,0,0,0); ONE = (1,0,0,0); PHI = (1,0,1,-1)

def zt(k): return POW[k % 10]
def zadd(a,b): return (a[0]+b[0],a[1]+b[1],a[2]+b[2],a[3]+b[3])
def zsub(a,b): return (a[0]-b[0],a[1]-b[1],a[2]-b[2],a[3]-b[3])
def zmul(a,b):
    c=[0]*8
    for i in range(4):
        if a[i]:
            for j in range(4): c[i+j]+=a[i]*b[j]
    r=[c[0],c[1],c[2],c[3]]
    r[3]+=c[4]; r[2]-=c[4]; r[1]+=c[4]; r[0]-=c[4]
    r[0]-=c[5]; r[1]-=c[6]; r[2]-=c[7]
    return tuple(r)
def zrot(a,k): return zmul(a,zt(k))
def zconj(a):
    r=ZERO
    for k in range(4):
        if a[k]: r=zadd(r,zmul((a[k],0,0,0),zt(-k)))
    return r
def norm2(a):
    z=zmul(a,zconj(a))
    assert z[1]==0 and z[3]==-z[2], "ノルムが実数でない"
    return (z[0]-z[2], z[2])                      # p + qφ

def phi_sign(x):
    """p+qφ の符号。φ=(1+√5)/2 なので 2p+q と q で判定する（浮動小数を使わない）"""
    A=2*x[0]+x[1]; B=x[1]
    if A==0 and B==0: return 0
    if A>=0 and B>=0: return 1
    if A<=0 and B<=0: return -1
    return (1 if A*A>5*B*B else -1) if A>0 else (-1 if A*A>5*B*B else 1)
def phi_lt(a,b): return phi_sign((a[0]-b[0], a[1]-b[1])) < 0

PHI2  = zmul(PHI,PHI)                              # 円環の半径
NCELL = norm2(PHI)                                 # 五角形の中心間の最小距離

# ── 接続は2種だけ ────────────────────────────────────────────
CONT = [zrot(zmul(PHI2, zsub(zt(4),ONE)), k) for k in range(10)]   # 連続 4.9798
SKIP = [zrot(zmul(PHI2, zsub(zt(2),ONE)), k) for k in range(10)]   # 2つ飛ばし 3.0777
N_CONT = norm2(CONT[0]); N_SKIP = norm2(SKIP[0])

# ── 描画座標（-18° 回して単位の軸を水平にする。y は画面下向き） ──
ROT = cmath.exp(-1j*math.radians(18))
def xy(a):
    z=sum(a[k]*cmath.exp(1j*math.pi*k/5) for k in range(4))*ROT
    return (z.real, z.imag)

# ── 単位（数の鎖） ───────────────────────────────────────────
# 奇数：連続接続だけ。折れ角は全箇所108°（=±72°ずつ交互に折る）
# 偶数：2つ飛ばしの「2」を連続接続でつなぐ。飛・連・飛・連…の直線
C_A, C_B = CONT[5], CONT[7]        # 342° と 54°。軸は 18°
S_AX, C_AX = SKIP[7], CONT[6]      # ともに 18°

def unit(n):
    cs=[ZERO]
    v = [C_A,C_B] if n%2 else [S_AX,C_AX]
    for i in range(n-1): cs.append(zadd(cs[-1], v[i%2]))
    return cs

def ring_cells(c): return [(zadd(c, zrot(PHI2,k)), k) for k in range(10)]

def fits(centers):
    """五角形を置く。同位置は番地の偶奇が合えば可、それ以外は距離 φ 未満を禁じる"""
    cells={}
    for c in centers:
        for q,k in ring_cells(c):
            if q in cells:
                if cells[q]%2 != k%2: return None
            else: cells[q]=k
    qs=list(cells)
    grid=defaultdict(list)
    for q in qs:
        p=xy(q); grid[(int(p[0]//2), int(p[1]//2))].append(q)
    for q in qs:
        p=xy(q); gx,gy=int(p[0]//2), int(p[1]//2)
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for q2 in grid.get((gx+dx,gy+dy),()):
                    if q2!=q and phi_lt(norm2(zsub(q,q2)), NCELL): return None
    return cells

# ── 積む：1 を一番上に、中央を揃えて、落差が最小になる位置へ ──
ROW_ORDER = list(range(1,14))

def build_stack():
    V1=(1,-1,0,-3); V2=(3,-3,0,-4)                 # 画面下向きの格子歩 3.618 / 5.854
    def smul(k,z): return (k*z[0],k*z[1],k*z[2],k*z[3])
    base=set()
    for a in CONT+SKIP:
        base.add(a)
        for b in CONT+SKIP:
            base.add(zadd(a,b))
            for c in CONT+SKIP: base.add(zadd(zadd(a,b),c))
    ver=sorted({(round(xy(zadd(smul(i,V1),smul(j,V2)))[1],9), zadd(smul(i,V1),smul(j,V2)))
                for i in range(10) for j in range(10)})
    rows=[unit(n) for n in ROW_ORDER]
    def mid(cs):
        x=[xy(c)[0] for c in cs]; return (min(x)+max(x))/2
    place=[rows[0]]; offs=[ZERO]; cur=ZERO
    for i in range(1,len(rows)):
        want = mid(place[-1]) - mid([zadd(c,cur) for c in rows[i]])
        cands={(round(xy(zadd(t,v))[1],9), zadd(t,v))
               for t in base if abs(xy(t)[0]-want)<1e-6 for _,v in ver}
        cands |= {(round(xy(t)[1],9), t) for t in base if abs(xy(t)[0]-want)<1e-6}
        ok=False
        for dy,t in sorted(x for x in cands if x[0]>1e-9):
            trial=[zadd(c, zadd(cur,t)) for c in rows[i]]
            if fits(sum(place,[])+trial) is None: continue
            cur=zadd(cur,t); place.append(trial); offs.append(cur); ok=True; break
        if not ok: raise RuntimeError("行 %d が置けない" % ROW_ORDER[i])
    return rows, place, offs

# ── 隙間：五角形の辺の使用回数だけから切り出す ────────────────
def gaps(cells):
    poly={q:[zadd(q,zt(a+2*i)) for i in range(5)] for q,a in cells.items()}
    use=Counter(); nbr=defaultdict(set)
    for q,vs in poly.items():
        for i in range(5):
            u,v=vs[i],vs[(i+1)%5]
            use[frozenset((u,v))]+=1
    for e,n in use.items():
        if n==1:
            u,v=tuple(e); nbr[u].add(v); nbr[v].add(u)
    ang={}
    for u in nbr:
        pu=xy(u)
        ang[u]=sorted(nbr[u], key=lambda w: math.atan2(xy(w)[1]-pu[1], xy(w)[0]-pu[0]))
    faces=[]; seen=set()
    for u in nbr:
        for v in nbr[u]:
            if (u,v) in seen: continue
            cyc=[]; a,b=u,v
            while True:
                seen.add((a,b)); cyc.append(a)
                ring=ang[b]; i=ring.index(a)
                a,b = b, ring[(i-1) % len(ring)]
                if (a,b)==(u,v): break
                if len(cyc)>2000: raise RuntimeError("面の追跡が閉じない")
            P=[xy(p) for p in cyc]
            s=sum(P[i][0]*P[(i+1)%len(P)][1]-P[(i+1)%len(P)][0]*P[i][1] for i in range(len(P)))/2
            if s<0 and abs(s)<100: faces.append((abs(s), cyc))   # 外周（面積2332）を除く
    return faces

GAP_NAME={0.8123:"細ひし形", 2.9389:"五芒星", 5.0656:"舟", 10.6331:"正十角形"}

# ── 2,3 の篩 ────────────────────────────────────────────────
def sieve_class(n):
    if n==1: return "unit"
    if n in (2,3): return "sieve%d"%n
    if n%2==0: return "out2"
    if n%3==0: return "out3"
    return "prime"
COLOR={"unit":"#555b66","sieve2":"#4aa3e0","sieve3":"#5cc06e",
       "out2":"#1e2c3a","out3":"#1e3324","prime":"#f2b544"}

# ── 検定 ────────────────────────────────────────────────────
def main(svg_path=None):
    ok=lambda c: "OK" if c else "NG"
    res=[]

    # 検定1・2
    n_ok=True; a_ok=True
    for n in ROW_ORDER:
        cs=unit(n); cells=fits(cs)
        if cells is None or len(cells)!=8*n+2: n_ok=False
        for i in range(1,len(cs)-1):
            u=complex(*xy(zsub(cs[i-1],cs[i]))); w=complex(*xy(zsub(cs[i+1],cs[i])))
            d=round(math.degrees(abs(cmath.phase(w/u))),6)
            if abs(d-(108 if n%2 else 180))>1e-6: a_ok=False
    res.append(("検定1 単位の枚数 8n+2", n_ok))
    res.append(("検定2 折れ角 奇数108°/偶数180°", a_ok))

    rows, place, offs = build_stack()
    allc=sum(place,[]); cells=fits(allc)
    res.append(("検定3 91環・628枚・重なり0",
                cells is not None and len(allc)==91 and len(cells)==628))

    # 検定4：外形の共線性。左端・右端の円環中心が一直線か
    def collinear(pts):
        """外積 Im(conj(d0)·d) = 0 を整数で判定する（実数元は w[1]==0 かつ w[3]==-w[2]）"""
        if len(pts)<3: return True
        d0=zsub(pts[1],pts[0])
        for p in pts[2:]:
            w=zmul(zconj(d0), zsub(p,pts[0]))
            if not (w[1]==0 and w[3]==-w[2]): return False
        return True
    L=[min(r,key=lambda c: xy(c)[0]) for r in place]
    R=[max(r,key=lambda c: xy(c)[0]) for r in place]
    res.append(("検定4 左右の端が共線", collinear(L) and collinear(R)))

    # 検定5・6・7
    fs=gaps(cells)
    kinds=Counter(round(a,4) for a,_ in fs)
    res.append(("検定5 隙間4種142個", len(fs)==142 and set(kinds)==set(GAP_NAME)))
    rowof={}
    for ri,r in enumerate(place):
        for c in r: rowof[c]=ri
    inner=[]; outer=0
    for a,cyc in fs:
        P=[xy(p) for p in cyc]
        g=(sum(p[0] for p in P)/len(P), sum(p[1] for p in P)/len(P))
        c=min(rowof, key=lambda c: math.dist(g, xy(c)))
        if math.dist(g, xy(c))<2.0: inner.append((rowof[c], round(a,4), cyc))
        else: outer+=1
    per=Counter(r for r,_,_ in inner)
    want={i: (ROW_ORDER[i] if ROW_ORDER[i]%2 else 3*ROW_ORDER[i]//2) for i in range(13)}
    res.append(("検定6 内側 奇数n/偶数3n/2・計112",
                len(inner)==112 and all(per[i]==want[i] for i in range(13))))
    res.append(("検定7 五芒星30個は外側",
                kinds[2.9389]==30 and all(a!=2.9389 for _,a,_ in inner)))

    for name,c in res: print("%-32s %s" % (name, ok(c)))
    print("\n行ごとの内側:", {ROW_ORDER[i]: per[i] for i in range(13)})
    print("隙間の内訳:", {GAP_NAME.get(k,str(k)): v for k,v in sorted(kinds.items())})
    print("篩の結果:", {n: sieve_class(n) for n in ROW_ORDER})
    print("\nNG %d / %d" % (sum(1 for _,c in res if not c), len(res)))

    if svg_path: write_svg(svg_path, cells, inner, place)

def write_svg(path, cells, inner, place):
    pent={q:[xy(zadd(q,zt(a+2*i))) for i in range(5)] for q,a in cells.items()}
    pts=[p for v in pent.values() for p in v]
    x0=min(p[0] for p in pts)-2; x1=max(p[0] for p in pts)+2
    y0=min(p[1] for p in pts)-2; y1=max(p[1] for p in pts)+2
    o=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0:.2f} {y0:.2f} {x1-x0:.2f} {y1-y0:.2f}">',
       f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{x1-x0:.2f}" height="{y1-y0:.2f}" fill="#0e1014"/>',
       '<g stroke="#3a414c" stroke-width="0.05" fill="none">']
    for v in pent.values():
        o.append('<polygon points="'+" ".join(f"{x:.4f},{y:.4f}" for x,y in v)+'"/>')
    o.append('</g><g stroke="none">')
    for ri,_,cyc in inner:
        d="M"+"L".join(f"{x:.4f} {y:.4f}" for x,y in (xy(p) for p in cyc))+"Z"
        o.append(f'<path d="{d}" fill="{COLOR[sieve_class(ROW_ORDER[ri])]}"/>')
    o.append('</g><g fill="#7b8492" font-family="system-ui" font-size="2.2">')
    for ri,r in enumerate(place):
        o.append(f'<text x="{x0+0.6:.2f}" y="{xy(r[0])[1]+0.8:.2f}">{ROW_ORDER[ri]}</text>')
    o.append('</g></svg>')
    open(path,"w").write("\n".join(o))
    print("→", path)

if __name__=="__main__":
    main(sys.argv[1] if len(sys.argv)>1 else None)
