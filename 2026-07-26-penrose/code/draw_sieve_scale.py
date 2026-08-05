#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
draw_sieve_scale.py ── 篩の階段と、二つの対数
  上：横軸 √n。段は p に立ち、次の段までの距離が ln p に伸びる
  下：横軸 √n。本数 π(√n) と 帯の周期の対数 θ(√n)。ln n は下に沈む
事前登録（判別法6）
  検定1 段の間隔の平均が ln p に一致するか（比が 1 前後に収まる）
  検定2 θ(√n)/√n が単調に増えて 1 に近づくか（下から寄る）
  検定3 対照 ── 本数 π(√n) は ln n では説明できないか（比が発散する）
"""
import math

def isp(n):
    if n < 2: return False
    d = 2
    while d*d <= n:
        if n % d == 0: return False
        d += 1
    return True

def sieve(N):
    b = bytearray([1])*(N+1); b[0]=b[1]=0
    for i in range(2, int(N**0.5)+1):
        if b[i]: b[i*i::i] = bytearray(len(b[i*i::i]))
    return [i for i in range(N+1) if b[i]]
BIG = sieve(10**6)
PR = [p for p in BIG if p <= 1000]
def pi(x):  return sum(1 for p in PR if p <= x)
def th(x):  return sum(math.log(p) for p in PR if p <= x)

# ── 検定 ─────────────────────────────────────────────
small = [p for p in PR if p <= 100]
gaps  = [small[i+1]-small[i] for i in range(len(small)-1)]
r1 = sum(g/math.log(p) for g, p in zip(gaps, small)) / len(gaps)
print("検定1 段の間隔 / ln p の平均 = %.3f          %s" % (r1, "OK" if 0.7 < r1 < 1.4 else "NG"))
sl = [th(x)/x for x in (100, 300, 600, 1000)]
mono = all(sl[i] < sl[i+1] for i in range(len(sl)-1)) and 0.8 < sl[0] and sl[-1] < 1.0
print("検定2 θ(√n)/√n = %s   %s"
      % (["%.3f" % v for v in sl], "OK 下から1へ" if mono else "NG"))
import bisect
def PI(x): return bisect.bisect_right(BIG, x)
rat = [PI(math.sqrt(N))/math.log(N) for N in (10**2, 10**4, 10**6, 10**8, 10**10, 10**12)]
print("検定3 対照 π(√N)/ln N = %s  %s"
      % (["%.1f" % v for v in rat], "OK 発散" if rat[-1] > 5*rat[0] else "NG"))

# ── 作図 ─────────────────────────────────────────────
GOLD, BLUE, GREEN, RED, GRAY = "#c8963c", "#2f5fa8", "#3f8f56", "#b4523c", "#8a8a8a"
W, H = 900, 780
L, R = 78, 46
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d">' % (W, H, W),
     '<rect width="%d" height="%d" fill="#faf8f4"/>' % (W, H),
     '<style>text{font-family:sans-serif}</style>']

# ---- 上：√n 0..100、段の間隔が ln p ----
X0, Y0, PW, PH = L, 60, W-L-R, 290
def ax(x): return X0 + x/100*PW
def ay(k): return Y0 + PH - k/26*PH
o.append('<text x="%d" y="34" font-size="15" fill="#333">段の間隔は ln p ── 横軸 √n（n は 10⁴ まで）</text>' % L)
o.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#bbb"/>' % (X0, Y0+PH, X0+PW, Y0+PH))
o.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#bbb"/>' % (X0, Y0, X0, Y0+PH))
pts, k = [], 0
for p in small:
    pts.append((ax(p), ay(k))); k += 1; pts.append((ax(p), ay(k)))
pts.append((ax(100), ay(k)))
o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2"/>'
         % (" ".join("%.1f,%.1f" % q for q in pts), BLUE))
for i, p in enumerate(small[:-1]):
    g, lg = small[i+1]-p, math.log(p)
    o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" opacity="0.30"/>'
             % (ax(p), Y0+PH+6, max(1.5, ax(small[i+1])-ax(p)), g*4.4, GREEN))
    o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.8"/>'
             % (ax(p), Y0+PH+6+lg*4.4, ax(small[i+1]), Y0+PH+6+lg*4.4, RED))
for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 61, 83):
    o.append('<text x="%.1f" y="%d" font-size="8.5" fill="%s" text-anchor="middle">%d</text>'
             % (ax(p), Y0+PH+70, GRAY, p*p))
o.append('<text x="%d" y="%d" font-size="9.5" fill="%s">緑＝実際の段の間隔　赤線＝ln p　'
         '下の数字＝段の立つ行 p²</text>' % (X0, Y0+PH+84, GRAY))
o.append('<text x="%d" y="%d" font-size="10" fill="%s" text-anchor="end">本数</text>' % (X0-8, Y0+10, GRAY))
for k in (0, 5, 10, 15, 20, 25):
    o.append('<text x="%d" y="%.1f" font-size="9" fill="%s" text-anchor="end">%d</text>' % (X0-8, ay(k)+3, GRAY, k))
x1, y1 = ax(math.sqrt(1103)), ay(11)
o.append('<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/>' % (x1, y1, GOLD))
o.append('<text x="%.1f" y="%.1f" font-size="10.5" fill="%s">1103 → 11本</text>' % (x1+8, y1-6, GOLD))

# ---- 下：√n 0..1000、二つの対数 ----
Y1, PH2 = 470, 230
def bx(x): return X0 + x/1000*PW
def byL(k): return Y1 + PH2 - k/175*PH2
def byR(v): return Y1 + PH2 - v/1050*PH2
o.append('<text x="%d" y="%d" font-size="15" fill="#333">帯の対数は √n の直線 ── 横軸 √n（n は 10⁶ まで）</text>' % (L, Y1-24))
o.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#bbb"/>' % (X0, Y1+PH2, X0+PW, Y1+PH2))
o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2.2"/>'
         % (" ".join("%.1f,%.1f" % (bx(x), byR(th(x))) for x in range(1, 1001, 5)), GREEN))
o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.1" stroke-dasharray="5 4"/>'
         % (" ".join("%.1f,%.1f" % (bx(x), byR(x)) for x in range(1, 1001, 20)), "#c9c9c9"))
o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2"/>'
         % (" ".join("%.1f,%.1f" % (bx(x), byL(pi(x))) for x in range(1, 1001)), BLUE))
o.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.6" stroke-dasharray="3 3"/>'
         % (" ".join("%.1f,%.1f" % (bx(x), byL(2*math.log(x))) for x in range(2, 1001, 4)), RED))
o.append('<text x="%.1f" y="%.1f" font-size="11" fill="%s">θ(√n) ＝ 帯の周期の対数（右）</text>' % (bx(430), byR(560), GREEN))
o.append('<text x="%.1f" y="%.1f" font-size="11" fill="%s">π(√n) ＝ 必要な本数（左）</text>' % (bx(600), byL(125), BLUE))
o.append('<text x="%.1f" y="%.1f" font-size="11" fill="%s">ln n ── 本数はこれでは説明できない</text>' % (bx(430), byL(28), RED))
for k in (0, 50, 100, 150):
    o.append('<text x="%d" y="%.1f" font-size="9" fill="%s" text-anchor="end">%d</text>' % (X0-8, byL(k)+3, BLUE, k))
for v in (0, 250, 500, 750, 1000):
    o.append('<text x="%d" y="%.1f" font-size="9" fill="%s">%d</text>' % (X0+PW+7, byR(v)+3, GREEN, v))
for x in (0, 200, 400, 600, 800, 1000):
    o.append('<text x="%.1f" y="%d" font-size="9" fill="%s" text-anchor="middle">√n=%d</text>'
             % (bx(x), Y1+PH2+16, GRAY, x))
    o.append('<text x="%.1f" y="%d" font-size="8" fill="%s" text-anchor="middle">n=%s</text>'
             % (bx(x), Y1+PH2+28, GRAY, "{:,}".format(x*x)))
o.append('<text x="%d" y="%d" font-size="10.5" fill="%s">'
         'θ(100)=83.7 は「1万まで篩う帯の周期 2.3×10³⁶」の対数そのもの</text>' % (L, H-16, GRAY))
o.append('</svg>')
open("/mnt/user-data/outputs/fig_sieve_scale.svg", "w").write("\n".join(o))
print("書き出しました %d×%d" % (W, H))
