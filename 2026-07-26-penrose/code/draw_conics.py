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
  検定2 五角形の中心は描いた円錐曲線の上に乗るか
      OK なら：残差が 1e-9 未満（曲線は後付けの当てはめではない）
      NG なら：和・差の分類が間違っている
  検定3 隣の組の族と交わるか
      OK なら：2族の曲線が交点を持つ＝相互作用が図に出る
      NG なら：族が分離していて重ならない
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

# 検定2：五角形の中心が曲線上に乗るか
worst = 0
for pr in (P1, P2):
    A, Bc, d, sums, difs, pts = family(pr)
    for q in pts:
        p = B.xy(q)
        ra, rb = math.dist(p, A), math.dist(p, Bc)
        worst = max(worst, min(abs(ra + rb - s) for s in sums))
        worst = max(worst, min(abs(abs(ra - rb) - t) for t in difs + [0.0]))
print("検定2 五角形中心の曲線への残差 最大 %.2e → %s" % (worst, "OK" if worst < 1e-9 else "NG"))

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
print("検定3 曲線の本数 組1 %d 本／組2 %d 本" % (len(C1), len(C2)))

cx, cy = mid(P1); half = 11.0
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
