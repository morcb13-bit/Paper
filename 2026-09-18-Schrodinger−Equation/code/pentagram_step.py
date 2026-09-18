"""
pentagram_step.py --- 五芒星グラフの上の一歩を整数だけで組む

事前登録（判別法6）

  検定P1 正十二面体20頂点を Z[phi] の整数対で作り、二乗距離を厳密に分類する
      OK なら：稜が 30本（d^2 = 8 - 4phi）、面の対角線が 60本（d^2 = 4）で、
               その二種しか出ない。長さの比が phi であることは二乗距離から出る
      NG なら：座標か分類が誤り

  検定P2 次数と一筆書き
      OK なら：五芒星グラフは 20頂点・60辺・次数6一様・連結・奇数次数0（一筆書き可）。
               稜のグラフは次数3・奇数次数20（不可）。v253 の記述が再現する
      NG なら：v253 の記述かこの構成のどちらかが誤り

  検定P3 一歩の作用素（有向辺120本）
      psi'(u->v) = (2L/d)*sum_{w~u} psi(w->u) - L*psi(v->u)、d = 6、L = 3
      OK なら：2L/d = 1 で、係数が 1 と -3 の二つだけ。割り算が式から消える。
               さらに T^T T = 9*I（非対角がすべて 0、対角がすべて 9）
      NG なら：整数で書けないか、直交でない

  検定P4 二乗和の倍率（平衡13進の桁列のまま）
      OK なら：一歩ごとに二乗和がちょうど 9 倍。12歩まで桁列の等式として成立
      NG なら：倍率が一定でないか、9 でない

  検定P5 負の対照
      (a) 係数 -3 を -2 に変えると倍率が一定でなくなること
      (b) 担体側の次数集合（v253 §5 の R^2=16、7種）では L がどれだけ大きくなるか
      OK なら：(a) が崩れ、(b) が桁数として大きい。単純さの代金の所在が言える（判別法1）
      NG なら：(a) が崩れないなら検定P4 は NG を返せない検査
"""

from itertools import combinations
from b13num import B13, Zphi

results = []


def report(name, ok, note=""):
    results.append((name, ok, note))
    print(f"{'OK ' if ok else 'NG '} {name}  {note}")


# ---------------- 検定P1 ----------------
# 座標はすべて Z[phi]。(a, b) は a + b*phi
Z = lambda a, b: (a, b)
I1 = Z(1, 0)
IM = Z(-1, 0)
PH = Z(0, 1)
PM = Z(0, -1)
IN = Z(-1, 1)      # 1/phi = phi - 1
INM = Z(1, -1)
O = Z(0, 0)


def zadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def zsub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def zmul(x, y):
    a, b = x
    c, d = y
    return (a * c + b * d, a * d + b * c + b * d)


verts = []
for sx in (1, -1):
    for sy in (1, -1):
        for sz in (1, -1):
            verts.append((Z(sx, 0), Z(sy, 0), Z(sz, 0)))
for s1 in (1, -1):
    for s2 in (1, -1):
        inv = Z(-s1, s1)        # s1 * (phi - 1)
        ph = Z(0, s2)           # s2 * phi
        verts.append((O, inv, ph))
        verts.append((inv, ph, O))
        verts.append((ph, O, inv))

assert len(verts) == 20


def d2(p, q):
    s = O
    for i in range(3):
        t = zsub(p[i], q[i])
        s = zadd(s, zmul(t, t))
    return s


dist_count = {}
for i, j in combinations(range(20), 2):
    dist_count[d2(verts[i], verts[j])] = dist_count.get(d2(verts[i], verts[j]), 0) + 1

EDGE = (8, -4)      # 8 - 4 phi
DIAG = (4, 0)       # 4
p1_ok = (dist_count.get(EDGE) == 30 and dist_count.get(DIAG) == 60)
report("検定P1 二乗距離の分類", p1_ok,
       f"稜 {dist_count.get(EDGE)}本(8-4phi) / 対角線 {dist_count.get(DIAG)}本(4) / 全 {len(dist_count)}種")

# 比が phi であること：DIAG = phi^2 * EDGE を整数対のまま確認
phi2 = (1, 1)       # phi^2 = 1 + phi
ratio_ok = zmul(EDGE, phi2) == DIAG
report("検定P1b 長さの比", ratio_ok, "対角線^2 = phi^2 * 稜^2（整数対のまま一致）")


# ---------------- 検定P2 ----------------
def build(target):
    adj = {i: [] for i in range(20)}
    for i, j in combinations(range(20), 2):
        if d2(verts[i], verts[j]) == target:
            adj[i].append(j)
            adj[j].append(i)
    return adj


star = build(DIAG)
edge_graph = build(EDGE)


def connected(adj):
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == len(adj)


star_deg = sorted({len(v) for v in star.values()})
edge_deg = sorted({len(v) for v in edge_graph.values()})
star_odd = sum(1 for v in star.values() if len(v) % 2)
edge_odd = sum(1 for v in edge_graph.values() if len(v) % 2)

p2_ok = (star_deg == [6] and edge_deg == [3]
         and star_odd == 0 and edge_odd == 20
         and connected(star) and connected(edge_graph))
report("検定P2 次数と一筆書き", p2_ok,
       f"五芒星 次数{star_deg} 奇数{star_odd}個 連結{connected(star)} / "
       f"稜 次数{edge_deg} 奇数{edge_odd}個")


# ---------------- 検定P3 ----------------
D = 6
L = 3
COEF = 2 * L // D          # = 1
assert 2 * L == COEF * D   # 割り切れる（割り算が式から消える）

darcs = []
index = {}
for u in range(20):
    for v in star[u]:
        index[(u, v)] = len(darcs)
        darcs.append((u, v))
N = len(darcs)


def step(psi):
    """psi'(u->v) = COEF * sum_{w~u} psi(w->u) - L * psi(v->u)"""
    out = [0] * N
    insum = [0] * 20
    for u in range(20):
        s = 0
        for w in star[u]:
            s += psi[index[(w, u)]]
        insum[u] = s
    for k, (u, v) in enumerate(darcs):
        out[k] = COEF * insum[u] - L * psi[index[(v, u)]]
    return out


gram_diag = set()
gram_off = set()
images = []
for k in range(N):
    e = [0] * N
    e[k] = 1
    images.append(step(e))
for a in range(N):
    for b in range(a, N):
        s = sum(x * y for x, y in zip(images[a], images[b]))
        if a == b:
            gram_diag.add(s)
        else:
            gram_off.add(s)

p3_ok = (N == 120 and COEF == 1 and gram_diag == {9} and gram_off == {0})
report("検定P3 T^T T = 9I", p3_ok,
       f"有向辺 {N}本 / 係数 (+{COEF}, -{L}) / 対角 {sorted(gram_diag)} / 非対角 {sorted(gram_off)}")


# ---------------- 検定P4 ----------------
B_ZERO = B13.from_int(0)


def bstep(psi):
    insum = [B_ZERO] * 20
    for u in range(20):
        s = B_ZERO
        for w in star[u]:
            s = s.add(psi[index[(w, u)]])
        insum[u] = s
    out = []
    for (u, v) in darcs:
        t = psi[index[(v, u)]]
        out.append(insum[u].sub(t.times_small(L)))
    return out


def bsumsq(psi):
    s = B_ZERO
    for x in psi:
        s = s.add(x.mul(x))
    return s


import random
random.seed(13)

cases = {
    "単一の有向辺に1": [1 if k == 0 else 0 for k in range(N)],
    "三値を配る": [random.choice((-1, 0, 1)) for _ in range(N)],
}

p4_ok = True
p4_notes = []
for name, init in cases.items():
    psi = [B13.from_int(x) for x in init]
    prev = bsumsq(psi)
    ratios = []
    for n in range(12):
        psi = bstep(psi)
        cur = bsumsq(psi)
        # cur == 9 * prev を桁列のまま照合（9倍は繰り返し加算で作る）
        nine = prev.times_small(6).add(prev.times_small(3))
        ratios.append(cur.eq(nine))
        prev = cur
    ok = all(ratios)
    p4_ok = p4_ok and ok
    p4_notes.append(f"{name}: 12歩中 {sum(ratios)}歩で 9倍ちょうど / 最終二乗和 {prev.to_int()}")
report("検定P4 二乗和の倍率（平衡13進）", p4_ok, " | ".join(p4_notes))


# ---------------- 検定P5 ----------------
def step_broken(psi, bad_L):
    insum = [0] * 20
    for u in range(20):
        insum[u] = sum(psi[index[(w, u)]] for w in star[u])
    return [COEF * insum[u] - bad_L * psi[index[(v, u)]] for (u, v) in darcs]


psi = [1 if k == 0 else 0 for k in range(N)]
bad_ratios = []
prev = sum(x * x for x in psi)
for _ in range(8):
    psi = step_broken(psi, 2)
    cur = sum(x * x for x in psi)
    bad_ratios.append(cur / prev if prev else None)
    prev = cur
p5a_ok = len(set(round(r, 6) for r in bad_ratios)) > 1
report("検定P5a 負の対照（係数を -2 に）", p5a_ok,
       "倍率 " + ", ".join(f"{r:.3f}" for r in bad_ratios))


def lcm(a, b):
    x, y = a, b
    while y:
        x, y = y, x % y
    return a * b // x


carrier_degs = [24, 26, 33, 34, 38, 45, 120]     # v253 §5 の R^2=16、7種
m = 1
for d in carrier_degs:
    m = lcm(m, d)
L_carrier = m if m % 2 else m // 2
L_star = L
p5b_ok = L_carrier > L_star
report("検定P5b 次数が不揃いなとき要る L", p5b_ok,
       f"担体 L = {L_carrier}（{len(str(L_carrier))}桁、一歩あたり L^2 = {L_carrier**2}） / "
       f"五芒星 L = {L_star}")

print()
n_ok = sum(1 for _, ok, _ in results if ok)
print(f"{n_ok}/{len(results)} OK")
