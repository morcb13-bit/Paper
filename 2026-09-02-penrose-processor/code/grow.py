# 担体の拡張（窓から直接生やす）
#
#   expand3 は円環を一つずつ幾何判定で足すので、半径を倍にすると重い。
#   ここでは今ある 1245環の担体から
#     (1) 五角形の向きごとの窓（垂直空間での五角形）
#     (2) 隣り合う五角形の差ベクトルと向きの変化
#   を取り出し、窓の中に入るかどうかだけで新しい五角形を足す。
#
#   検算  出来た担体が、半径110の内側で元の担体と一字一致すること
import json, math, pickle
from collections import Counter, defaultdict
import b13_chain_units as U

TARGET = 260.0          # 目標半径（φ¹⁰ = 122.99 の外側に十分な余裕）

cells = {tuple(int(x) for x in k.split(",")): a
         for k, a in json.load(open("carrier_1245.json"))["cells"].items()}
def conj(v):
    out = (0, 0, 0, 0)
    for i, co in enumerate(v):
        if co: out = U.zadd(out, tuple(co * x for x in U.zt((3 * i) % 10)))
    return out
XY = U.xy

# --- 窓：向きごとに、垂直空間の点の凸包を五角形として取る ---
def hull(pts):
    pts = sorted(set(pts))
    def cross(o, a, b): return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    lo = []
    for p in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0: lo.pop()
        lo.append(p)
    up = []
    for p in reversed(pts):
        while len(up) >= 2 and cross(up[-2], up[-1], p) <= 0: up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]

WIN = {}
for a in range(10):
    pts = [XY(conj(q)) for q, b in cells.items() if b == a]
    h = hull(pts)
    WIN[a] = h
    print("向き %d  窓の頂点 %d 個  半径 %.4f" % (a, len(h), max(math.hypot(*p) for p in h)))

# 窓を法線形式にする（少しだけ外に広げて、有限標本で角が削れる分を戻す）
EPS = 1e-3
LINES = {}
for a, h in WIN.items():
    L = []
    for i in range(len(h)):
        x1, y1 = h[i]; x2, y2 = h[(i + 1) % len(h)]
        nx, ny = y2 - y1, x1 - x2
        n = math.hypot(nx, ny); nx, ny = nx / n, ny / n
        c = nx * x1 + ny * y1
        L.append((nx, ny, c + EPS))
    LINES[a] = L
def inwin(a, p):
    return all(nx * p[0] + ny * p[1] <= c for nx, ny, c in LINES[a])

# --- 隣接の型：差ベクトルと向きの変化 ---
XYc = {q: XY(q) for q in cells}
byloc = defaultdict(list)
for q in cells:
    x, y = XYc[q]; byloc[(round(x / 2), round(y / 2))].append(q)
OFF = set()
for q in cells:
    x, y = XYc[q]
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for r in byloc.get((round(x / 2) + dx, round(y / 2) + dy), []):
                if r == q: continue
                d = math.hypot(XYc[r][0] - x, XYc[r][1] - y)
                if d < 2.0: OFF.add((U.zsub(r, q), cells[q], cells[r]))
print("隣接の型 %d 種" % len(OFF))
OFFBY = defaultdict(list)
for d, a, b in OFF: OFFBY[a].append((d, b))

# --- 広げる ---
new = dict(cells)
frontier = [q for q in cells if math.hypot(*XYc[q]) > 100.0]
it = 0
while frontier:
    it += 1
    nxt = []
    for q in frontier:
        for d, b in OFFBY[new[q]]:
            c = U.zadd(q, d)
            if c in new: continue
            x, y = XY(c)
            if math.hypot(x, y) > TARGET: continue
            if not inwin(b, XY(conj(c))): continue
            new[c] = b; nxt.append(c)
    frontier = nxt
    if it % 5 == 0 or not nxt:
        print("  %2d 巡目  五角形 %d 個" % (it, len(new)))
print("拡張後 五角形 %d 個（元 %d）" % (len(new), len(cells)))

# --- 検算：半径110の内側で元と一致するか ---
old_in = {q for q in cells if math.hypot(*XY(q)) < 110.0}
new_in = {q for q in new if math.hypot(*XY(q)) < 110.0}
print("半径110の内側  元 %d / 新 %d  一致 %s  向きも一致 %s"
      % (len(old_in), len(new_in), old_in == new_in,
         all(new[q] == cells[q] for q in old_in & new_in)))
json.dump({"cells": {",".join(map(str, k)): v for k, v in new.items()}},
          open("carrier_big.json", "w"))
print("carrier_big.json を書き出し")
