#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
draw_conics.py ── 2つ飛ばしの円環対を焦点に、楕円と双曲線を描く
============================================================
2026-07-03 の「半径10・中心距離10の二中心」モデルを、ペンローズの中に置き直す。
和が一定 → 楕円、差が一定 → 双曲線。隣の2つ飛ばし組の族を重ねる。

事前登録（判別法6）
  検定1 2つ飛ばし接続は環あたり何本か
      OK なら：ちょうど1本＝210組が環を共有しない完全マッチング
      NG なら：2本以上持つ環がある＝族が焦点を共有する
  検定2 五角形の中心は描いた円錐曲線の上に乗るか（多重度で見る）
      OK なら：18枚が10本に分かれ、和側・差側とも多重度が {4×4本, 2×1本}
               ＝1本あたり4枚が乗る（共有五角形の2本だけ2枚）。曲線は当てはめではない
      NG なら：多重度が1に散る＝点ごとに1本ずつ曲線ができているだけ
      対照（判別法2）：同じ焦点対に格子から外した18点を置くと NG に落ちること
      注：残差では見ない。和と差をその点集合自身から集めている以上、
          残差は round(...,9) の丸めまでしか上がらず、どんな入力でも OK になる。
          旧版の「残差 5.00e-10」はでたらめな20点でも同じ桁が出る（誤り55 の具体形）
  検定3 隣の組の族と交わるか
      OK なら：描画窓の中で組1と組2の曲線が交点を持つ＝相互作用が図に出る
      NG なら：族が分離していて重ならない
      対照（判別法2）：中点が最も遠い組を相手にすると、同じ窓で交点0＝NG になること
"""
import math, sys
from collections import Counter, defaultdict
import b13_chain_units as B

PHI = (1 + 5 ** 0.5) / 2
T = (-2, 2, 0, 3); S = (-5, 5, 0, 8)
rows, place, offs = B.build_stack()
U = []; seen = set()
for k in range(10):
    t = B.zrot(T if k % 2 == 0 else S, k)
    for c in sum(place, []):
        c2 = B.zadd(B.zrot(c, k), t)
        if c2 not in seen: seen.add(c2); U.append(c2)
R = set(U); cells = B.fits(U)

pairs = []
for c in U:
    for v in B.SKIP:
        o = B.zadd(c, v)
        if o in R and (o, c) not in pairs: pairs.append((c, o))
deg = Counter()
for a, b in pairs: deg[a] += 1; deg[b] += 1
print("検定1 2つ飛ばし %d 組／環あたりの本数 %s → %s"
      % (len(pairs), dict(Counter(deg.values())), "OK" if set(deg.values()) == {1} else "NG"))

# 原点に近い組と、その隣の組
pairs.sort(key=lambda p: math.hypot(*B.xy(p[0])) + math.hypot(*B.xy(p[1])))
P1 = pairs[0]
mid = lambda p: ((B.xy(p[0])[0] + B.xy(p[1])[0]) / 2, (B.xy(p[0])[1] + B.xy(p[1])[1]) / 2)
P2 = min((p for p in pairs[1:]), key=lambda p: math.dist(mid(p), mid(P1)))
print("組1 中点 (%.3f, %.3f) ／ 組2 中点 (%.3f, %.3f)  間隔 %.3f"
      % (*mid(P1), *mid(P2), math.dist(mid(P1), mid(P2))))

def family(pair):
    A, Bc = B.xy(pair[0]), B.xy(pair[1])
    d = math.dist(A, Bc)
    pts = set()
    for c in pair:
        for q, k in B.ring_cells(c): pts.add(q)
    sums, difs = set(), set()
    for q in pts:
        p = B.xy(q)
        sums.add(round(math.dist(p, A) + math.dist(p, Bc), 9))
        difs.add(round(abs(math.dist(p, A) - math.dist(p, Bc)), 9))
    return A, Bc, d, sorted(sums), sorted(x for x in difs if x > 1e-9), pts

# 検定2：五角形の中心が曲線に乗ることを、1本あたりの枚数で見る
def multiplicity(A, Bc, pts_xy):
    """焦点 A,B から見た和・差の多重度を返す"""
    cs, cd = Counter(), Counter()
    for p in pts_xy:
        ra, rb = math.dist(p, A), math.dist(p, Bc)
        cs[round(ra + rb, 9)] += 1
        cd[round(abs(ra - rb), 9)] += 1
    return cs, cd

def shape(cnt):
    """多重度の分布を {枚数: 本数} で返す"""
    return dict(sorted(Counter(cnt.values()).items()))

WANT = {2: 1, 4: 4}          # 4枚の曲線が4本＋2枚の曲線が1本
m_ok = True; m_rep = []
for nm, pr in (("組1", P1), ("組2", P2)):
    A, Bc, d, sums, difs, pts = family(pr)
    cs, cd = multiplicity(A, Bc, [B.xy(q) for q in pts])
    good = (shape(cs) == WANT and shape(cd) == WANT)
    m_ok &= good
    m_rep.append("%s 和%s 差%s" % (nm, shape(cs), shape(cd)))
print("検定2 曲線1本あたりの五角形 %s → %s" % ("／".join(m_rep), "OK" if m_ok else "NG"))

# 対照：格子から外した18点を同じ焦点対に置くと落ちること（判別法2）
import random as _rnd
_rnd.seed(0)
A, Bc, d, *_ = family(P1)
_pts = [(_rnd.uniform(A[0] - 6, A[0] + 6), _rnd.uniform(A[1] - 6, A[1] + 6)) for _ in range(18)]
_cs, _cd = multiplicity(A, Bc, _pts)
print("      対照 格子外18点 和%s 差%s → %s"
      % (shape(_cs), shape(_cd), "落ちる（検査は有効）" if shape(_cs) != WANT else "落ちない（反証不能）"))

# ── 描画 ───────────────────────────────────────────────────
def conic_paths(A, Bc, d, sums, difs, N=240):
    mx, my = (A[0] + Bc[0]) / 2, (A[1] + Bc[1]) / 2
    th = math.atan2(Bc[1] - A[1], Bc[0] - A[0])
    ct, st = math.cos(th), math.sin(th)
    def to(x, y): return (mx + x * ct - y * st, my + x * st + y * ct)
    out = []
    c = d / 2
    for s in sums:
        a = s / 2
        if a <= c + 1e-12: continue
        b = math.sqrt(a * a - c * c)
        out.append(("e", [to(a * math.cos(2 * math.pi * i / N), b * math.sin(2 * math.pi * i / N))
                          for i in range(N + 1)]))
    for t in difs:
        a = t / 2
        if a >= c - 1e-12: continue
        b = math.sqrt(c * c - a * a)
        for sg in (1, -1):
            out.append(("h", [to(sg * a * math.cosh(u), b * math.sinh(u))
                              for u in [(-2.2 + 4.4 * i / N) for i in range(N + 1)]]))
    return out

F1 = family(P1); F2 = family(P2)
C1 = conic_paths(*F1[:5]); C2 = conic_paths(*F2[:5])

# 検定3：描画窓の中で2族が交わるか（本数を数えるのではなく交点を数える）
cx, cy = mid(P1); half = 11.0
BOX = (cx - half, cy - half, cx + half, cy + half)

def _inbox(p, box=BOX):
    return box[0] <= p[0] <= box[2] and box[1] <= p[1] <= box[3]

def _cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

def _hit(p1, p2, p3, p4):
    d1, d2 = _cross(p3, p4, p1), _cross(p3, p4, p2)
    d3, d4 = _cross(p1, p2, p3), _cross(p1, p2, p4)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0))

def _segs(C):
    """描画窓の中に入る線分だけを、1.0刻みの格子に振り分ける"""
    g = defaultdict(list)
    for _, pts in C:
        for i in range(len(pts) - 1):
            a, b = pts[i], pts[i + 1]
            if not (_inbox(a) or _inbox(b)):
                continue
            for gx in range(int(min(a[0], b[0]) // 1), int(max(a[0], b[0]) // 1) + 1):
                for gy in range(int(min(a[1], b[1]) // 1), int(max(a[1], b[1]) // 1) + 1):
                    g[(gx, gy)].append((a, b))
    return g

def n_cross(CA, CB):
    ga, gb = _segs(CA), _segs(CB)
    seen = set()
    for cell, sa in ga.items():
        for a, b in sa:
            for c, dpt in gb.get(cell, ()):
                if _hit(a, b, c, dpt):
                    k = (round(a[0], 6), round(a[1], 6), round(c[0], 6), round(c[1], 6))
                    seen.add(k)
    return len(seen)

n12 = n_cross(C1, C2)
# 対照：中点が最も遠い組を相手にすると同じ窓で交わらないこと（判別法2）
PF = max(pairs, key=lambda p: math.dist(mid(p), mid(P1)))
CF = conic_paths(*family(PF)[:5])
n1f = n_cross(C1, CF)
print("検定3 組1×組2 の交点 %d 個（窓 %.1f 角）→ %s" % (n12, 2 * half, "OK" if n12 > 0 else "NG"))
print("      対照 組1×最遠の組（中点間 %.2f）の交点 %d 個 → %s"
      % (math.dist(mid(PF), mid(P1)), n1f, "落ちる（検査は有効）" if n1f == 0 else "落ちない"))
x0, y0 = cx - half, cy - half
o = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="%.2f %.2f %.2f %.2f">' % (x0, y0, 2 * half, 2 * half),
     '<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="#0e1014"/>' % (x0, y0, 2 * half, 2 * half),
     '<defs><clipPath id="v"><rect x="%.2f" y="%.2f" width="%.2f" height="%.2f"/></clipPath></defs>'
     % (x0, y0, 2 * half, 2 * half),
     '<g stroke="#2f3742" stroke-width="0.035" fill="none">']
for q, a in cells.items():
    v = [B.xy(B.zadd(q, B.zt(a + 2 * i))) for i in range(5)]
    if min(p[0] for p in v) > x0 - 3 and max(p[0] for p in v) < x0 + 2 * half + 3 \
       and min(p[1] for p in v) > y0 - 3 and max(p[1] for p in v) < y0 + 2 * half + 3:
        o.append('<polygon points="%s"/>' % " ".join("%.4f,%.4f" % p for p in v))
o.append('</g><g clip-path="url(#v)" fill="none" stroke-width="0.06">')
for col, C in (("#f2b544", C1), ("#4aa3e0", C2)):
    for kind, pts in C:
        dsh = ' stroke-dasharray="0.35 0.3"' if kind == "h" else ""
        o.append('<polyline points="%s" stroke="%s" opacity="%s"%s/>'
                 % (" ".join("%.4f,%.4f" % p for p in pts), col, "0.95" if kind == "e" else "0.7", dsh))
o.append('</g><g>')
for col, pr in (("#f2b544", P1), ("#4aa3e0", P2)):
    for c in pr:
        p = B.xy(c); o.append('<circle cx="%.4f" cy="%.4f" r="0.16" fill="%s"/>' % (*p, col))
o.append('</g></svg>')
out = sys.argv[1] if len(sys.argv) > 1 else "b13_conics.svg"
open(out, "w").write("\n".join(o))
print("→", out, " 実線＝楕円（和一定）／破線＝双曲線（差一定）／金＝組1・青＝組2")
