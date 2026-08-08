#  対象を一行で（判別法9）：
#      中心のペンタゴン（網膜の中心）を意図的に格子ベクトル d だけずらしたとき、
#      ずらす前と後の両方で成立している場所の集合が何になるか。
#
#  機構：場所 v が在る ⟺ 残渣 σ(v) が窓 W の中。ずらした側に在る ⟺ σ(v)−σ(d) が W の中。
#        両方成立 ⟺ σ(v) ∈ W ∩ (W+σ(d))。ずらしたぶんだけ窓が細くなる。
#
#  座標：X = Re σ(v)、Y = Im σ(v)/sin36°。どちらも Q(φ) の元。
#        これは平面の線形写像なので、凸性・内外・面積比はそのまま保たれる。
#        すべて厳密。浮動小数は表示のときだけ（判別法16）。
#
#  検定SW1 一致集合は新しい窓で書けるか
#      OK なら：担体を直に突き合わせた集合と、窓の判定だけで作った集合が一致する。
#               ずらしは窓を作る操作である
#      NG なら：窓では書けない
#      反証の用意：格子に乗らないずらしでは一致0になるはず（対照を最後に置く）
#
#  検定SW4 窓は五角形か
#      OK なら：担体の残渣の外郭がちょうど5頂点
#      NG なら：5頂点にならない
#
#  検定SW8 ずらした窓に新しい辺の向きが出るか
#      OK なら：出ない。辺の向きは元の窓の5通りだけ
#      NG なら：出る（窓の向きがずらしで変わる）
#
#  検定SW5 閉じる条件は W+(−W) で書けるか
#      OK なら：W∩(W+σ(d)) が空 ⟺ σ(d) が W+(−W) の外。例外なし
#      NG なら：例外がある
#
#  検定SW9 窓が開く実長の下限
#      OK なら：実長×残渣=√N と閉じる境界から出る下限 √N/R₁₀ を、
#               担体の最短の一歩がすべて上回る
#      NG なら：下回るものがある
#
#  検定SW6 窓の広さは残渣の長さだけで決まるか
#      OK なら：残渣の長さが同じなら面積比も同じ（等方）
#      NG なら：向きが効く
#
#  検定SW7 面積比は一致の割合を予言するか
#      OK なら：担体で数えた割合と面積比の差が 0.02 未満
#      NG なら：合わない（担体が有限で画角の効果が混ざる場合を含む）
#
#  使い方: b13_two_tilings.py と qphi.py と同じ場所で python3 shift_window.py

from fractions import Fraction as F
import b13_two_tilings as T
from qphi import Qp, zmul, zconj, zsigma, zre

# ── 担体は道具から作り、整合を通してから使う（v206） ──────────────────
FIG = T.build()
V = T.verify(FIG)
assert V["重なり0"] and V["全て等半径"], "担体が整合していない"
rings = [tuple(F(x) for x in v) for v in FIG.rings]
S = set(rings)
print(f"担体：{FIG.stats()}  整合 OK")

def zim(a):     return Qp(a[1], a[2] + a[3])
def zsub(a, b): return tuple(x - y for x, y in zip(a, b))
def norm2(a):   return zre(zmul(a, zconj(a)))
def res2(a):
    s = zsigma(a); return zre(zmul(s, zconj(s)))
def absN(a):
    x = norm2(a); return int(x.p*x.p + x.p*x.q - x.q*x.q)
def XY(v):
    s = zsigma(v); return (zre(s), zim(s))

# ── 凸多角形の道具（すべて Q(φ)） ─────────────────────────────────
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
            out.append((p[0] + (q[0]-p[0])*t, p[1] + (q[1]-p[1])*t))
    return out
def inter(A, B):
    poly = A
    for i in range(len(B)):
        if not poly: return []
        poly = clip(poly, B[i], B[(i+1) % len(B)])
    return poly
def shift(poly, s): return [(p[0]+s[0], p[1]+s[1]) for p in poly]
def edge_dirs(poly):
    out = set()
    for i in range(len(poly)):
        a, b = poly[i], poly[(i+1) % len(poly)]
        dx, dy = b[0]-a[0], b[1]-a[1]
        if dx.sign() == 0 and dy.sign() == 0: continue
        k = dx if dx.sign() != 0 else dy
        if k.sign() < 0: dx, dy, k = -dx, -dy, -k
        out.add((str(dx/k), str(dy/k)))
    return out
def inside(poly, p):
    n = len(poly)
    return all(cross(poly[i], poly[(i+1) % n], p).sign() >= 0 for i in range(n))

RES = {v: XY(v) for v in rings}
W = hull(list(RES.values()))
AW = area2(W)
RSET = set(RES.values())

# ── 検定SW4 ─────────────────────────────────────────────────
print(f"\n検定SW4 窓の外郭の頂点数 = {len(W)}")
for p in W:
    print(f"   ({p[0].val():+.6f}, {p[1].val():+.6f})   厳密 X={p[0]}  Y={p[1]}")
print(f"検定SW4: {'OK' if len(W)==5 else 'NG'}   窓の面積(2×) = {AW} = {AW.val():.6f}\n")

# ── ずらしの候補（担体の差ベクトル＝格子に乗るずらし） ─────────────────
cand = {}
for a in rings:
    for b in rings:
        if a == b: continue
        d = zsub(a, b)
        cand.setdefault(norm2(d), []).append(d)
D = [(n, vs) for n, vs in sorted(cand.items(), key=lambda kv: kv[0].val())]
Rmax = max(norm2(r).val() for r in rings)

# ── 検定SW1・SW7 ─────────────────────────────────────────────
print("検定SW1／SW7")
print(f"{'実長':>9}{'残渣':>9}{'ノルム':>7}{'辺':>4}{'面積比':>9}"
      f"{'一致':>6}{'窓から':>7}{'母数':>6}{'割合':>8}")
bad1 = 0; dev = []
for n, vs in D[:16]:
    d = vs[0]; s = XY(d)
    I = inter(W, shift(W, s))
    ratio = (area2(I)/AW).val() if I else 0.0
    A = {v for v in rings if zsub(v, d) in S}
    B = {v for v in rings if (RES[v][0]-s[0], RES[v][1]-s[1]) in RSET}
    if A != B: bad1 += 1
    denom = sum(1 for v in rings if norm2(zsub(v, d)).val() <= Rmax)
    frac = len(A)/denom if denom else 0.0
    if ratio > 0: dev.append(abs(frac-ratio))
    print(f"{n.val()**0.5:9.4f}{res2(d).val()**0.5:9.4f}{absN(d):7d}{len(I):4d}"
          f"{ratio:9.4f}{len(A):6d}{len(B):7d}{denom:6d}{frac:8.4f}")
print(f"検定SW1: {'OK' if bad1==0 else 'NG'}  食い違い {bad1} 件")
print(f"検定SW7: {'OK' if max(dev) < 0.02 else 'NG'}  面積比と割合の差の最大 {max(dev):.4f}\n")

# ── 検定SW8・SW5・SW6 ────────────────────────────────────────
DW = edge_dirs(W)
mW = [(-p[0], -p[1]) for p in W]
DB = hull([(a[0]+b[0], a[1]+b[1]) for a in W for b in mW])
extra = 0; tested = 0; sides = set(); bad5 = 0; bykey = {}
for n, vs in D[:80]:
    for d in vs[:4]:
        s = XY(d)
        I = inter(W, shift(W, s))
        nonempty = bool(I) and area2(I).sign() > 0
        if nonempty and not inside(DB, (-s[0], -s[1])): bad5 += 1
        if not nonempty: continue
        tested += 1; sides.add(len(I))
        if not edge_dirs(I) <= DW: extra += 1
        bykey.setdefault(round(res2(d).val(), 9), set()).add(round((area2(I)/AW).val(), 6))
print(f"検定SW8: {'OK' if extra==0 else 'NG'}  {tested} 本中 新しい辺の向き {extra} 件"
      f"（ずらした窓の辺の数は {sorted(sides)}）")
print(f"検定SW5: {'OK' if bad5==0 else 'NG'}  W+(−W) の頂点数 = {len(DB)}  例外 {bad5} 件")
aniso = {k: v for k, v in bykey.items() if len(v) > 1}
print(f"検定SW6: {'OK（等方）' if not aniso else 'NG（向きが効く）'}  "
      f"長さが同じで開き方が違う残渣値 {len(aniso)} 通り")

R10 = max((p[0].val()**2 + (p[1].val()*0.5877852522924731)**2)**0.5 for p in DB)
print(f"   W+(−W) の外接半径 R₁₀ = {R10:.6f}\n")

# ── 検定SW9 ─────────────────────────────────────────────────
print("検定SW9 窓が開く実長の下限")
seen = {}
for n, vs in D:
    for d in vs:
        N = absN(d); L = n.val()**0.5
        if N not in seen or L < seen[N]: seen[N] = L
bad9 = 0
for N in sorted(seen)[:8]:
    lim = (N**0.5)/R10
    if seen[N] <= lim: bad9 += 1
    print(f"   ノルム {N:>4}  下限 √N/R₁₀ = {lim:8.4f}   担体の最短 {seen[N]:8.4f}")
print(f"検定SW9: {'OK' if bad9==0 else 'NG'}  下回るもの {bad9} 件\n")

# ── 対照：格子に乗らないずらし（反証の用意） ────────────────────────
print("対照 格子に乗らないずらし")
for d in [(1,0,0,0), (0,1,0,0), (1,-1,0,0), (0,0,1,0)]:
    d = tuple(F(x) for x in d)
    m = sum(1 for v in rings if zsub(v, d) in S)
    print(f"   {tuple(int(x) for x in d)}  実長 {norm2(d).val()**0.5:.4f}"
          f"  残渣 {res2(d).val()**0.5:.4f}  一致 {m}")
