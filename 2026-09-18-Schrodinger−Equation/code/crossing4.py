"""
crossing4.py --- Q が閉じない件の決着

検定F3 が NG を返した。Q' は A と Q の一次結合にならない。
足りないものを式から探す。

  S_u = sum_{w~u} psi(w->u)      番地 u に入る辺の和
  O_u = sum_{w~u} psi(u->w)      番地 u から出る辺の和

一歩 psi'(u->v) = k S_u - L psi(v->u) を代入すると

  S'_u = k * (A S)_u - L * O_u        A はグラフの隣接
  O'_u = (k d - L) * S_u

事前登録（判別法6）

  検定G1  (S, O) が閉じるか
      OK なら：S' = k A S - L O、O' = (kd - L) S が乱数の psi と複数の (d,k,L) で厳密成立。
               120本の辺の上の一歩が、番地あたり2個・計40個の整数で閉じる
      NG なら：式が違う。さらに量が要る

  検定G2  S だけの二階漸化式
      OK なら：S_{n+2} = k A S_{n+1} - L(kd - L) S_n が厳密成立
      NG なら：二階にならない

  検定G3  隣接の固有値と判別式
      D(lambda) = lambda^2 - 4 L (kd - L)
      OK なら：五芒星グラフの固有値 6, 2, 1, -2, -3 に対し、
               k=1, L=3, d=6 では lambda=6 で D=0（境界）、他の四つで D<0（回る側）。
               根の積が L(kd-L) なので、回る側の半径は sqrt(L(kd-L))
      NG なら：分類がそうならない

  検定G4  保存と半径の一致
      OK なら：kd = 2L のとき、かつそのときに限り L(kd-L) = L^2、
               すなわち回る側の半径が L に一致する
      NG なら：二つの条件が別物

  検定G5  負の対照
      OK なら：係数を 1 ずらすと検定G1・G2 の等式が崩れる
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


graphs = {
    "d=2 円": {i: [(i - 1) % 20, (i + 1) % 20] for i in range(20)},
    "d=3 稜": build((8, -4)),
    "d=6 五芒星": build((4, 0)),
}


class Walk:
    def __init__(self, adj, k, L):
        self.adj, self.k, self.L = adj, k, L
        self.d = len(adj[next(iter(adj))])
        self.c = self.d * k * k - 2 * L * k
        self.darcs, self.index = [], {}
        for u in sorted(adj):
            for v in adj[u]:
                self.index[(u, v)] = len(self.darcs)
                self.darcs.append((u, v))
        self.N = len(self.darcs)

    def S(self, psi):
        return [sum(psi[self.index[(w, u)]] for w in self.adj[u]) for u in sorted(self.adj)]

    def Oo(self, psi):
        return [sum(psi[self.index[(u, w)]] for w in self.adj[u]) for u in sorted(self.adj)]

    def step(self, psi):
        S = self.S(psi)
        return [self.k * S[u] - self.L * psi[self.index[(v, u)]] for (u, v) in self.darcs]

    def adjmul(self, vec):
        return [sum(vec[w] for w in self.adj[u]) for u in sorted(self.adj)]


# ---------------- 検定G1 ----------------
g1_ok = True
count = 0
for name, adj in graphs.items():
    for k in range(-2, 5):
        for L in (1, 2, 3, 6):
            W = Walk(adj, k, L)
            for _ in range(3):
                psi = [random.randint(-4, 4) for _ in range(W.N)]
                p2 = W.step(psi)
                S, Ov = W.S(psi), W.Oo(psi)
                S2_pred = [k * x - L * o for x, o in zip(W.adjmul(S), Ov)]
                O2_pred = [(k * W.d - L) * x for x in S]
                count += 1
                if W.S(p2) != S2_pred or W.Oo(p2) != O2_pred:
                    g1_ok = False
report("検定G1 (S, O) が閉じる", g1_ok,
       f"{count} 件で厳密。120本の辺 → 番地あたり2個・計40個の整数")


# ---------------- 検定G2 ----------------
g2_ok = True
for name, adj in graphs.items():
    for k, L in ((1, 3), (2, 3), (1, 6), (3, 3)):
        W = Walk(adj, k, L)
        psi = [random.randint(-4, 4) for _ in range(W.N)]
        seq = [W.S(psi)]
        for _ in range(10):
            psi = W.step(psi)
            seq.append(W.S(psi))
        coef = L * (k * W.d - L)
        for n in range(len(seq) - 2):
            pred = [k * x - coef * y for x, y in zip(W.adjmul(seq[n + 1]), seq[n])]
            if pred != seq[n + 2]:
                g2_ok = False
report("検定G2 二階漸化式", g2_ok,
       "S_{n+2} = k A S_{n+1} - L(kd-L) S_n（10歩、3グラフ × 4組）")


# ---------------- 検定G3 ----------------
star = graphs["d=6 五芒星"]
n = 20
Amat = [[1 if j in star[i] else 0 for j in range(n)] for i in range(n)]


def matmul(X, Y):
    return [[sum(X[i][t] * Y[t][j] for t in range(n)) for j in range(n)] for i in range(n)]


def matsub_scalar(X, s):
    return [[X[i][j] - (s if i == j else 0) for j in range(n)] for i in range(n)]


eigs = [6, 2, 1, -2, -3]
P = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
for lam in eigs:
    P = matmul(P, matsub_scalar(Amat, lam))
zero_ok = all(all(x == 0 for x in row) for row in P)

k, L, d = 1, 3, 6
disc = {lam: lam * lam - 4 * L * (k * d - L) for lam in eigs}
g3_ok = (zero_ok and disc[6] == 0 and all(disc[l] < 0 for l in eigs if l != 6))
print()
print("  五芒星グラフ  k=1, L=3, d=6 のとき 4L(kd-L) = 36")
for lam in eigs:
    kind = "境界（重根）" if disc[lam] == 0 else ("回る" if disc[lam] < 0 else "伸びる")
    print(f"    lambda = {lam:>2}   D = {disc[lam]:>4}   {kind}")
report("検定G3 判別式による分類", g3_ok,
       f"Π(A - λI) = 0 は {zero_ok} / λ=6 のみ D=0、他は D<0")


# ---------------- 検定G4 ----------------
g4_rows = []
g4_ok = True
for d in (2, 3, 6, 12):
    for L in (1, 2, 3, 6):
        for k in range(0, 5):
            prod = L * (k * d - L)
            cond_cons = (k * d == 2 * L)
            cond_rad = (prod == L * L)
            if cond_cons != cond_rad:
                g4_ok = False
            if cond_cons:
                g4_rows.append((d, k, L))
report("検定G4 保存と半径", g4_ok,
       f"kd = 2L ⟺ L(kd-L) = L^2 を全組で確認。該当 {len(g4_rows)} 組")


# ---------------- 検定G5 ----------------
W = Walk(star, 1, 3)
broke1 = broke2 = 0
for _ in range(20):
    psi = [random.randint(-4, 4) for _ in range(W.N)]
    p2 = W.step(psi)
    S, Ov = W.S(psi), W.Oo(psi)
    if W.S(p2) != [W.k * x - (W.L + 1) * o for x, o in zip(W.adjmul(S), Ov)]:
        broke1 += 1
    if W.Oo(p2) != [(W.k * W.d - W.L + 1) * x for x in S]:
        broke2 += 1
report("検定G5 負の対照", broke1 == 20 and broke2 == 20,
       f"係数を 1 ずらすと 20/20 件で崩れる（S: {broke1}, O: {broke2}）")

print()
n_ok = sum(1 for _, ok, _ in results if ok)
print(f"{n_ok}/{len(results)} OK")
