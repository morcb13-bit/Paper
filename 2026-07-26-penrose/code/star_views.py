#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
star_views.py ── 910環の図で、五芒星ごとの「景色」を分類する
============================================================
窓（担体）を使わない。10枚を組んだ図そのものの上で数える。

事前登録（判別法6）
  検定1 五芒星のまわりの円環の数
      OK なら：外周で切れていない五芒星はすべて φ³ に円環をちょうど5個持つ
      NG なら：5以外が出る＝「常に5つで囲まれる」が偽
  検定2 72°回転で自分に写る五芒星
      OK なら：中心の1個だけ
      NG なら：中心以外にも対称点がある＝「中心からのみ回転対称」が偽
  検定3 景色の種類
      見る半径を広げたとき、類の数が増えて飽和するか。
      すべて別々に割れるなら「それぞれ違った形」が成立する。
      割れずに残る対があれば、それは合同な景色＝分類の単位になる
  検定4 外周で切れた五芒星を除く条件が効いているか
      条件を外すと類の数が増えるはず（切れ方の違いが署名に入るため）。
      増えなければ除外条件が効いていない
"""
import math, sys
from collections import Counter, defaultdict
import b13_chain_units as B

PHI = (1 + 5 ** 0.5) / 2
T = (-2, 2, 0, 3)     # |T| = φ³
S = (-5, 5, 0, 8)     # |S| = φ⁵

rows, place, offs = B.build_stack()
U = []; seen = set()
for k in range(10):
    t = B.zrot(T if k % 2 == 0 else S, k)
    for c in sum(place, []):
        c2 = B.zadd(B.zrot(c, k), t)
        if c2 not in seen: seen.add(c2); U.append(c2)
cells = B.fits(U)
RC = [B.xy(c) for c in U]
print("環 %d / 五角形 %d" % (len(U), len(cells)))

# ── 隙間と外周 ─────────────────────────────────────────────
poly = {q: [B.zadd(q, B.zt(a + 2 * i)) for i in range(5)] for q, a in cells.items()}
use = Counter(); nbr = defaultdict(set)
for q, vs in poly.items():
    for i in range(5):
        use[frozenset((vs[i], vs[(i + 1) % 5]))] += 1
for e, n in use.items():
    if n == 1:
        u, v = tuple(e); nbr[u].add(v); nbr[v].add(u)
ang = {u: sorted(nbr[u], key=lambda w: math.atan2(B.xy(w)[1] - B.xy(u)[1],
                                                  B.xy(w)[0] - B.xy(u)[0])) for u in nbr}
faces = []; visited = set()
for u in nbr:
    for v in nbr[u]:
        if (u, v) in visited: continue
        cyc = []; a, b = u, v
        while True:
            visited.add((a, b)); cyc.append(a)
            ring = ang[b]; i = ring.index(a)
            a, b = b, ring[(i - 1) % len(ring)]
            if (a, b) == (u, v): break
        P = [B.xy(p) for p in cyc]
        s = sum(P[i][0] * P[(i + 1) % len(P)][1] - P[(i + 1) % len(P)][0] * P[i][1]
                for i in range(len(P))) / 2
        faces.append((s, cyc))
outer = max(faces, key=lambda f: abs(f[0]))[1]
OB = [B.xy(p) for p in outer]
print("外周の頂点 %d 個" % len(OB))

stars = []
for s, cyc in faces:
    if s < 0 and abs(abs(s) - 2.9389) < 1e-3:
        P = [B.xy(p) for p in cyc]
        stars.append((sum(p[0] for p in P) / len(P), sum(p[1] for p in P) / len(P)))
print("五芒星 %d 個" % len(stars))

def dist_to_outer(p):
    m = 1e9
    for i in range(len(OB)):
        a, b = OB[i], OB[(i + 1) % len(OB)]
        vx, vy = b[0] - a[0], b[1] - a[1]
        L = vx * vx + vy * vy
        t = 0 if L == 0 else max(0, min(1, ((p[0]-a[0])*vx + (p[1]-a[1])*vy) / L))
        m = min(m, math.hypot(p[0]-a[0]-t*vx, p[1]-a[1]-t*vy))
    return m

safe = [(p, dist_to_outer(p)) for p in stars]

# ── 検定1 まわりの円環 ─────────────────────────────────────
cnt = Counter()
for p, d in safe:
    if d < 6: continue
    cnt[sum(1 for q in RC if abs(math.dist(p, q) - PHI ** 3) < 1e-6)] += 1
print("\n検定1 φ³ にある円環の数（外周から6以上離れた五芒星）:", dict(cnt),
      "→", "OK" if set(cnt) == {5} else "NG")

# ── 署名 ───────────────────────────────────────────────────
def signature(p, R):
    rel = [(q[0]-p[0], q[1]-p[1]) for q in RC if math.dist(p, q) <= R + 1e-9]
    base = [math.atan2(y, x) for x, y in rel if abs(math.hypot(x, y) - PHI**3) < 1e-6]
    best = None
    for b in base:
        for mir in (1, -1):
            sig = tuple(sorted((round(math.hypot(x, y), 4),
                                round(math.degrees((mir*(math.atan2(y, x) - b))) % 360, 3))
                               for x, y in rel))
            if best is None or sig < best: best = sig
    return best

def rot_sym(p, R):
    rel = {(round(q[0]-p[0], 6), round(q[1]-p[1], 6)) for q in RC if math.dist(p, q) <= R + 1e-9}
    c, s = math.cos(math.radians(72)), math.sin(math.radians(72))
    return all((round(x*c - y*s, 6), round(x*s + y*c, 6)) in rel for x, y in rel)

R = 22.0
pool = [p for p, d in safe if d >= R]
print("外周から %.0f 以上離れた五芒星: %d 個" % (R, len(pool)))
sym = [p for p in pool if rot_sym(p, R)]
print("検定2 72°回転で自分に写る:", len(sym), "個",
      [tuple(round(v, 4) for v in p) for p in sym],
      "→", "OK" if len(sym) == 1 and math.hypot(*sym[0]) < 1e-6 else "NG")

print("\n検定3 景色の種類（D₅ を法にした類）")
for R2 in (7.0, 10.0, 15.0, 22.0):
    q = [p for p, d in safe if d >= R2]
    cls = Counter(signature(p, R2) for p in q)
    print("  見る半径 %4.1f  対象 %3d 個 → %3d 類   最大の類 %d 個  単独の類 %d"
          % (R2, len(q), len(cls), max(cls.values()), sum(1 for v in cls.values() if v == 1)))

print("\n検定4 外周の切れを除かない場合（半径15）")
cls_all = Counter(signature(p, 15.0) for p, d in safe)
cls_in  = Counter(signature(p, 15.0) for p, d in safe if d >= 15.0)
print("  除外なし 対象 %d 個 → %d 類 ／ 除外あり 対象 %d 個 → %d 類  → %s"
      % (len(safe), len(cls_all), sum(cls_in.values()), len(cls_in),
         "OK" if len(cls_all) > len(cls_in) else "NG"))
