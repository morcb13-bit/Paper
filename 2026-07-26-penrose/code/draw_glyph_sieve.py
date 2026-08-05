#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
draw_glyph_sieve.py ── 三角形を字形だけに簡素化した篩の図
  1 = ・   2 = ー   3 = v とその逆さ ∧
事前登録（判別法6）
  検定1 字形の個数が n を復元するか（3×V + 桁 = n、偶数は 2×ー = n）
  検定2 v と ∧ の交互が破れる箇所
      OK なら：桁0 の行では0箇所、桁±1 の行ではちょうど1箇所で、
               破れる位置が中心の印の位置と一致する
      NG なら：一致しない n を出す
  検定3 字形だけで落とせる行と、落ちずに残る合成数
      OK なら：ー を含む行＝2の倍数、中心に印が無い行＝3の倍数。
               それ以外に残る最小の合成数は 25
      NG なら：合わない n を出す
"""
import sys, math
sys.path.insert(0, "/home/claude/Paper/2026-07-26-penrose/code")
sys.path.insert(0, "/mnt/skills/user/b13-verify/scripts")
import chain_balanced as C

N   = int(sys.argv[1]) if len(sys.argv) > 1 else 120
U   = 7.0        # 環1個ぶんの横幅
AMP = 2.8        # 山と谷の高さ差
RH  = 7.0        # 行の高さ（横幅と揃えて三角形を正方形の枠に収める）
GAP = 0.5        # 字形どうしの隙間
BAND = 12        # 帯の周期

def isp(n):
    if n < 2: return False
    d = 2
    while d*d <= n:
        if n % d == 0: return False
        d += 1
    return True

def row_glyphs(n):
    """(種類, 開始環, 長さ) の列。種類は 'v' / '^' / '-' / '.' / 'o'"""
    if n == 1: return [(".", 1, 1)]
    if n % 2 == 0:
        return [("-", 2*i+1, 2) for i in range(n//2)]
    d, cov, cent = C.balanced(n)
    g = []
    for a, b in cov:
        g.append(("v" if a % 2 == 1 else "^", a, 3))
    if d == +1: g.append((".", cent, 1))
    if d == -1: g.append(("o", cent, 1))
    return sorted(g, key=lambda t: (t[1], t[0]))

# ── 検定 ────────────────────────────────────────────────
ok1 = ok2 = True
for n in range(1, N+1):
    g = row_glyphs(n)
    if n % 2 == 0:
        if 2*sum(1 for t in g if t[0] == "-") != n: ok1 = False
    elif n > 1:
        d, cov, cent = C.balanced(n)
        if 3*sum(1 for t in g if t[0] in "v^") + d != n: ok1 = False
        vs = [(t[1], t[0]) for t in g if t[0] in "v^"]
        br = [i for i in range(len(vs)-1) if vs[i][1] == vs[i+1][1]]
        if len(br) != (0 if d == 0 else 1): ok2 = False
        if d != 0 and br:
            a1, a2 = vs[br[0]][0], vs[br[0]+1][0]
            if not (a1 + 2 <= cent <= a2 + 2 and a2 - 2 <= cent <= a1 + 4): ok2 = False
print("検定1 字形の個数が n を復元する        %s" % ("OK" if ok1 else "NG"))
print("検定2 交互の破れは印の位置に1箇所だけ  %s" % ("OK" if ok2 else "NG"))

def falls(n):
    """字形だけで落とせるか"""
    if n == 1: return "one"
    g = row_glyphs(n)
    if any(t[0] == "-" for t in g): return "two"          # ー を含む＝偶数
    if not any(t[0] in ".o" for t in g): return "three"   # 中心に印が無い＝3の倍数
    return "kept"
kept_comp = [n for n in range(5, N+1) if falls(n) == "kept" and not isp(n)]
print("検定3 字形で落ちずに残る合成数         %s" % kept_comp[:10])

# ── 作図 ────────────────────────────────────────────────
COL = {"prime": "#c8963c", "two": "#2f5fa8", "three": "#3f8f56",
       "one": "#999999", "comp": "#b4523c"}
RED = "#d02020"
def color(n):
    if n == 1: return COL["one"]
    if isp(n):  return COL["prime"]
    f = falls(n)
    return COL["two"] if f == "two" else (COL["three"] if f == "three" else COL["comp"])

W = N*U + 90; H = N*RH + 70
out = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="%d %d %d %d" width="%d">'
       % (-W//2, -30, W, H, int(W)),
       '<rect x="%d" y="-30" width="%d" height="%d" fill="#faf8f4"/>' % (-W//2, W, H)]
for b in range(0, N//BAND + 1):
    if b % 2: continue
    out.append('<rect x="%d" y="%.1f" width="%d" height="%.1f" fill="#000" opacity="0.028"/>'
               % (-W//2, b*BAND*RH + RH/2, W, BAND*RH))
for n in range(1, N+1):
    y = n*RH; c = color(n); x0 = -(n-1)*U/2
    def px(i): return x0 + (i-1)*U
    def py(i): return y + (AMP if i % 2 == 0 else 0)
    for kind, a, ln in row_glyphs(n):
        if kind == "-":
            out.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" '
                       'stroke-width="1.5" stroke-linecap="round"/>'
                       % (px(a)+GAP, y+AMP/2, px(a+1)-GAP, y+AMP/2, c))
        elif kind in "v^":
            out.append('<polyline points="%.2f,%.2f %.2f,%.2f %.2f,%.2f" fill="none" '
                       'stroke="%s" stroke-width="1.25" stroke-linejoin="round"/>'
                       % (px(a)+GAP, py(a), px(a+1), py(a+1), px(a+2)-GAP, py(a+2), c))
        elif kind == ".":
            out.append('<circle cx="%.2f" cy="%.2f" r="1.5" fill="%s"/>' % (px(a), py(a), c))
        else:
            out.append('<circle cx="%.2f" cy="%.2f" r="1.9" fill="none" stroke="%s" '
                       'stroke-width="1.1"/>' % (px(a), py(a), c))
    p = isp(n)
    out.append('<text x="%.2f" y="%.2f" font-family="sans-serif" font-size="%s" '
               'fill="%s" %stext-anchor="end">%d</text>'
               % (x0-7, y+AMP/2+2.2, "6.5" if p else "5.5", RED if p else "#9a9a9a",
                  'font-weight="bold" ' if p else "", n))
out.append('<text x="%d" y="%.1f" font-family="sans-serif" font-size="11" fill="#555">'
           '1=・  2=ー  3=v/∧  ｜ 中心の印 ・=桁+1  ○=桁−1  ｜ 帯は12行'
           '</text>' % (-W//2+18, N*RH+40))
out.append('<text x="%d" y="%.1f" font-family="sans-serif" font-size="11" fill="#555">'
           '赤い数字=素数（初出）  金=素数  青=2で落ちる  緑=3で落ちる  朱=字形では落ちない合成数'
           '</text>' % (-W//2+18, N*RH+56))
out.append('</svg>')
open("/mnt/user-data/outputs/fig_glyph_sieve.svg", "w").write("\n".join(out))
print("行 %d  図の大きさ %d×%d  字形 %d 個（環なら %d 個ぶん）"
      % (N, W, H, sum(len(row_glyphs(n)) for n in range(1, N+1)), N*(N+1)//2))
