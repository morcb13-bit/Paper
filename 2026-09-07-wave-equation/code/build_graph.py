"""expanded.json から担体グラフ（頂点＝五角形の角、辺＝五角形の辺）を組む。
   v245 §2-7 の素性を再現できるかを最初に確かめる。"""
import json, cmath, math
from b13_two_tilings import Zeta, zadd, zt, to_xy, Figure
from b13_layers import LayeredFigure

D = json.load(open('expanded.json'))
R = [tuple(int(x) for x in c) for c in D['rings']]
Fg = LayeredFigure()
for c in R:
    Fg.add_ring(c)
print(f"円環{len(Fg.rings)} 五角形{len(Fg.cells)}")

# 角の取り方を決める：向き o の五角形の角は zt(2j+o)
def corners(q, o):
    return [zadd(q, zt(2 * j + o)) for j in range(5)]

for flip in (0, 1):
    pent = {}
    for q, a in Fg.cells.items():
        pent[q] = corners(q, (a + flip) % 2)
    # 隣接（中心間距離 φ）の五角形が角を2つ共有するか
    cs = list(Fg.cells.items())
    xy = {q: to_xy(q) for q, _ in cs}
    share2 = share_other = 0
    keys = list(xy)
    grid = {}
    for q in keys:
        p = xy[q]
        grid.setdefault((int(p.real // 2), int(p.imag // 2)), []).append(q)
    for q in keys:
        p = xy[q]
        gx, gy = int(p.real // 2), int(p.imag // 2)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for r in grid.get((gx + dx, gy + dy), []):
                    if r == q:
                        continue
                    if abs(xy[r] - p) < 1.7:      # φ=1.618 の隣接だけ
                        n = len(set(pent[q]) & set(pent[r]))
                        if n == 2:
                            share2 += 1
                        else:
                            share_other += 1
    print(f"flip={flip}: 角2つ共有 {share2//2} 組 / それ以外 {share_other//2} 組")
    if share_other == 0 and share2 > 0:
        good = flip
        break

pent = {q: corners(q, (a + good) % 2) for q, a in Fg.cells.items()}

# 頂点と辺
V = {}
adj = {}
def vid(v):
    if v not in V:
        V[v] = len(V)
        adj[V[v]] = set()
    return V[v]

for q, cor in pent.items():
    ids = [vid(c) for c in cor]
    for i in range(5):
        a, b = ids[i], ids[(i + 1) % 5]
        adj[a].add(b)
        adj[b].add(a)

deg = {}
for i in adj:
    deg[len(adj[i])] = deg.get(len(adj[i]), 0) + 1
directed = sum(len(adj[i]) for i in adj)
xyv = [to_xy(v) for v in V]
rad = max(abs(p) for p in xyv)
# 辺の長さ
a0 = next(iter(adj))
b0 = next(iter(adj[a0]))
elen = abs(xyv[a0] - xyv[b0])

print(f"頂点{len(V)} 有向辺{directed} 次数分布{dict(sorted(deg.items()))} "
      f"辺長{elen:.4f} 最大半径{rad:.1f}")

json.dump({'xy': [[p.real, p.imag] for p in xyv],
           'adj': [sorted(adj[i]) for i in range(len(V))]},
          open('carrier_1245_graph.json', 'w'))
