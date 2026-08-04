#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
similar_probe.py ── 型紙（五芒星＋最初の数殻）を当てて回る
============================================================
型紙 = 三角形の頂点にある五芒星を原点に置き、そこから行1..k の円環中心を並べたもの。
これを図の中の全ての五芒星に、10通りの向き × φ^m 倍で当てる。

事前登録（判別法6）
  検定1 型紙は自分自身に当たるか（m=0・中心の五芒星）
      OK なら：型紙の作り方が正しい
      NG なら：座標の取り方が違う
  検定2 他の五芒星に当たるか、当たるならどの倍率か
      当たる五芒星の個数を m ごとに数える。
      m=0 だけなら平行移動しただけで相似ではない。
      m≧1 が出れば相似形が実在する
  検定3 当たった向きは図の中心を向いているか
      型紙の軸と、その五芒星から見た中心の向きの差を数える
"""
import math, sys
from collections import Counter, defaultdict
import b13_chain_units as B

PHI = (1 + 5 ** 0.5) / 2
T = (-2, 2, 0, 3)      # 五芒星 → 行1の円環（|T| = φ³、画面方位90°）
S = (-5, 5, 0, 8)

rows, place, offs = B.build_stack()
U = []; seen = set()
for k in range(10):
    t = B.zrot(T if k % 2 == 0 else S, k)
    for c in sum(place, []):
        c2 = B.zadd(B.zrot(c, k), t)
        if c2 not in seen: seen.add(c2); U.append(c2)
RSET = set(U)
print("環 %d" % len(RSET))

# ── 五芒星の位置（格子点として厳密に） ─────────────────────
OFF = [B.zrot(T, j) for j in range(10)]        # 長さφ³の10方向
cand = set()
for c in RSET:
    for v in OFF: cand.add(B.zsub(c, v))
STARS = []
for g in cand:
    for j in range(10):
        ring = [B.zadd(g, B.zrot(OFF[j], 2 * i)) for i in range(5)]
        if all(r in RSET for r in ring):
            STARS.append(g); break
STARS = sorted(set(STARS))
print("五芒星 %d 個（格子点として検出）" % len(STARS))

# ── 型紙 ───────────────────────────────────────────────────
def template(k):
    return [B.zadd(c, T) for r in place[:k] for c in r]

PHIZ = B.PHI
def phipow(m):
    p = B.ONE
    for _ in range(m): p = B.zmul(p, PHIZ)
    return p

def hits(tmpl, m, g):
    """五芒星 g に型紙を φ^m 倍で当てられる向き（0..9）を返す"""
    s = phipow(m); out = []
    for r in range(10):
        ok = True
        for v in tmpl:
            if B.zadd(g, B.zrot(B.zmul(s, v), r)) not in RSET: ok = False; break
        if ok: out.append(r)
    return out

ORIGIN = (0, 0, 0, 0)
for k in (1, 2, 3, 4):
    tmpl = template(k)
    print("\n=== 型紙：行1〜%d（円環 %d 個）===" % (k, len(tmpl)))
    ctr = hits(tmpl, 0, ORIGIN)
    print("検定1 中心の五芒星に m=0 で当たる向き: %s → %s"
          % (ctr, "OK" if ctr else "NG"))
    for m in range(0, 4):
        s = phipow(m)
        good = [(g, hits(tmpl, m, g)) for g in STARS]
        good = [(g, h) for g, h in good if h]
        if not good:
            print("  m=%d (φ^%d=%.4f 倍)  当たる五芒星 0 個" % (m, m, PHI ** m)); continue
        rad = sorted({round(math.hypot(*B.xy(g)), 4) for g, _ in good})
        # 中心方向との差
        dif = Counter()
        for g, h in good:
            gx, gy = B.xy(g)
            if abs(gx) < 1e-9 and abs(gy) < 1e-9: dif["中心そのもの"] += 1; continue
            outward = math.degrees(math.atan2(gy, gx))       # g から見て外向き
            for r in h:
                axis = (90 + 36 * r) % 360                    # 型紙の軸（画面方位）
                d = round((axis - outward) % 360)
                dif[d] += 1
        print("  m=%d (φ^%d=%7.4f 倍)  当たる五芒星 %3d 個  半径 %s"
              % (m, m, PHI ** m, len(good), " ".join("%.3f" % r for r in rad[:6])))
        print("      型紙の軸と外向きの差:", dict(sorted(dif.items(), key=lambda x: str(x[0]))))
