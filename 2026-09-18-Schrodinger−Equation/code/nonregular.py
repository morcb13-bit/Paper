"""
nonregular.py --- 次数が不揃いなグラフでも (S, O) が閉じるか

対象を一行で書く（判別法9）：
  いま見ているのは「担体側の一歩」。次数が番地ごとに違う場合に、
  crossing4.py で五芒星グラフ（次数6一様）について出した閉じ方が、そのまま立つか。

事前登録（判別法6）

  検定H0  A' = L^2 A
      次数が不揃いでも c_u = k_u(d_u k_u - 2L) = 0 なので二乗和はちょうど L^2 倍のはず
      OK なら：v253 §6-2 の「T^T T = L^2 I は次数が不揃いでも成り立つ」が再現する
      NG なら：どちらかが誤り

  検定H1  (S, O) の一歩
      K = diag(k_u)、k_u = 2L/d_u として
      OK なら：S' = A K S - L O、O' = L S が厳密成立。O' から次数が消える
      NG なら：非正則では閉じない。五芒星グラフの結果は一様な場合に限る

  検定H2  二階漸化式
      OK なら：S_{n+2} = A K S_{n+1} - L^2 S_n が厳密成立
      NG なら：二階にならない

  検定H3  境界に乗る成分
      OK なら：正則グラフでは全番地一様なベクトルが A K の固有ベクトル。
               非正則では固有ベクトルにならない（境界に乗る成分が別物になる）
      NG なら：両者が同じ成分

  検定H4  負の対照
      OK なら：係数を 1 ずらすと H1 の等式が崩れる
      NG なら：検査になっていない
"""

from itertools import combinations
import random

random.seed(13)
results = []


def report(name, ok, note=""):
    results.append((name, ok, note))
    print(f"{'OK ' if ok else 'NG '} {name}  {note}")


def zmul(x, y):
    a, b = x
    c, d = y
    return (a * c + b * d, a * d + b * c + b * d)


def zsub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def zadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


O0 = (0, 0)
verts = []
for sx in (1, -1):
    for sy in (1, -1):
        for sz in (1, -1):
            verts.append(((sx, 0), (sy, 0), (sz, 0)))
for s1 in (1, -1):
    for s2 in (1, -1):
        inv, ph = (-s1, s1), (0, s2)
        verts += [(O0, inv, ph), (inv, ph, O0), (ph, O0, inv)]


def d2(p, q):
    s = O0
    for i in range(3):
        t = zsub(p[i], q[i])
        s = zadd(s, zmul(t, t))
    return s


def build(target):
    adj = {i: [] for i in range(20)}
    for i, j in combinations(range(20), 2):
        if d2(verts[i], verts[j]) == target:
            adj[i].append(j)
            adj[j].append(i)
    return adj


star = build((4, 0))


def drop_matching(adj, size):
    """辺を size 本、頂点を共有しないように落として次数を不揃いにする"""
    g = {u: list(v) for u, v in adj.items()}
    used, dropped = set(), []
    for u in sorted(g):
        if len(dropped) >= size:
            break
        if u in used:
            continue
        for v in g[u]:
            if v not in used and v > u:
                g[u].remove(v)
                g[v].remove(u)
                used.add(u)
                used.add(v)
                dropped.append((u, v))
                break
    return g


def random_graph(n, seed):
    rnd = random.Random(seed)
    g = {i: [] for i in range(n)}
    for i in range(1, n):                      # まず木で連結にする
        j = rnd.randrange(i)
        g[i].append(j)
        g[j].append(i)
    for _ in range(n):
        a, b = rnd.randrange(n), rnd.randrange(n)
        if a != b and b not in g[a]:
            g[a].append(b)
            g[b].append(a)
    return g


def lcm(a, b):
    x, y = a, b
    while y:
        x, y = y, x % y
    return a * b // x


class Walk:
    def __init__(self, adj):
        self.adj = adj
        self.deg = {u: len(adj[u]) for u in adj}
        m = 1
        for d in set(self.deg.values()):
            m = lcm(m, d)
        self.L = m if m % 2 else m // 2
        assert all(2 * self.L % d == 0 for d in self.deg.values())
        self.k = {u: 2 * self.L // self.deg[u] for u in adj}
        self.darcs, self.index = [], {}
        for u in sorted(adj):
            for v in adj[u]:
                self.index[(u, v)] = len(self.darcs)
                self.darcs.append((u, v))
        self.N = len(self.darcs)

    def S(self, psi):
        return {u: sum(psi[self.index[(w, u)]] for w in self.adj[u]) for u in self.adj}

    def Oo(self, psi):
        return {u: sum(psi[self.index[(u, w)]] for w in self.adj[u]) for u in self.adj}

    def step(self, psi):
        S = self.S(psi)
        return [self.k[u] * S[u] - self.L * psi[self.index[(v, u)]] for (u, v) in self.darcs]

    def AK(self, vec):
        """(A K vec)_u = sum_{w~u} k_w * vec_w"""
        return {u: sum(self.k[w] * vec[w] for w in self.adj[u]) for u in self.adj}

    def A_(self, vec):
        return {u: sum(vec[w] for w in self.adj[u]) for u in self.adj}


tests = {
    "五芒星（6一様）": Walk(star),
    "五芒星から5本落とす（5と6）": Walk(drop_matching(star, 5)),
    "乱数グラフ12頂点": Walk(random_graph(12, 7)),
    "乱数グラフ16頂点": Walk(random_graph(16, 11)),
}

print("  グラフ                        次数            L     k_u")
for name, W in tests.items():
    degs = sorted(set(W.deg.values()))
    ks = sorted(set(W.k.values()))
    print(f"  {name:<28} {str(degs):<14} {W.L:<5} {ks}")
print()


# ---------------- 検定H0 ----------------
h0_ok = True
for name, W in tests.items():
    for _ in range(5):
        psi = [random.randint(-4, 4) for _ in range(W.N)]
        if sum(x * x for x in W.step(psi)) != W.L ** 2 * sum(x * x for x in psi):
            h0_ok = False
report("検定H0 A' = L^2 A", h0_ok, "次数が不揃いでも二乗和はちょうど L^2 倍")


# ---------------- 検定H1 ----------------
h1_ok = True
cnt = 0
for name, W in tests.items():
    for _ in range(5):
        psi = [random.randint(-4, 4) for _ in range(W.N)]
        p2 = W.step(psi)
        S, Ov = W.S(psi), W.Oo(psi)
        AKS = W.AK(S)
        S2_pred = {u: AKS[u] - W.L * Ov[u] for u in W.adj}
        O2_pred = {u: W.L * S[u] for u in W.adj}
        cnt += 1
        if W.S(p2) != S2_pred or W.Oo(p2) != O2_pred:
            h1_ok = False
report("検定H1 (S, O) の一歩", h1_ok,
       f"S' = A K S - L O / O' = L S を {cnt} 件で厳密確認。O' から次数が消える")


# ---------------- 検定H2 ----------------
h2_ok = True
for name, W in tests.items():
    psi = [random.randint(-4, 4) for _ in range(W.N)]
    seq = [W.S(psi)]
    for _ in range(8):
        psi = W.step(psi)
        seq.append(W.S(psi))
    for n in range(len(seq) - 2):
        AKS = W.AK(seq[n + 1])
        pred = {u: AKS[u] - W.L ** 2 * seq[n][u] for u in W.adj}
        if pred != seq[n + 2]:
            h2_ok = False
report("検定H2 二階漸化式", h2_ok, "S_{n+2} = A K S_{n+1} - L^2 S_n（8歩、4グラフ）")


# ---------------- 検定H3 ----------------
print()
print("  一様ベクトル (すべて1) に A K を当てた結果")
h3_rows = []
for name, W in tests.items():
    ones = {u: 1 for u in W.adj}
    img = W.AK(ones)
    vals = sorted(set(img.values()))
    regular = len(set(W.deg.values())) == 1
    is_eigen = len(vals) == 1
    h3_rows.append((name, regular, is_eigen, vals))
    print(f"  {name:<28} 正則={regular}  値の種類={len(vals)}  {vals[:5]}")
h3_ok = all(reg == eig for _, reg, eig, _ in h3_rows)
report("検定H3 境界に乗る成分", h3_ok,
       "正則なら一様ベクトルは A K の固有ベクトル、非正則ならそうでない")


# ---------------- 検定H4 ----------------
W = tests["五芒星から5本落とす（5と6）"]
broke = 0
for _ in range(20):
    psi = [random.randint(-4, 4) for _ in range(W.N)]
    p2 = W.step(psi)
    S, Ov = W.S(psi), W.Oo(psi)
    AKS = W.AK(S)
    if W.S(p2) != {u: AKS[u] - (W.L + 1) * Ov[u] for u in W.adj}:
        broke += 1
report("検定H4 負の対照", broke == 20, f"L を 1 ずらすと 20件中 {broke} 件で崩れる")

print()
n_ok = sum(1 for _, ok, _ in results if ok)
print(f"{n_ok}/{len(results)} OK")
