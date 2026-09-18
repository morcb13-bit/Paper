"""
crossing3.py --- 根の外での二乗和の行方（整数だけで書けるか）

検定E4 で出た形
    T^T T = L^2 I + c(k)*G,   c(k) = d k^2 - 2 L k,  G は同じ頭を持つ有向辺の対
から、二乗和の一歩が読める。G の二次形式は

    psi^T G psi = sum_u ( sum_{w~u} psi(w->u) )^2

で、番地ごとの「入る辺の和」の二乗和。これを Q と書く。Q は整数。

事前登録（判別法6）

  検定F1  二乗和の一歩
      A = sum psi^2、Q = sum_u S_u^2（S_u は番地 u に入る辺の和）として
      OK なら：A' = L^2 * A + c(k) * Q が、乱数の psi と複数の (d, k, L) で厳密に成立
      NG なら：式が違う

  検定F2  固有空間の分解
      OK なら：T^T T の固有値は L^2（重複 有向辺数 - 番地数）と L^2 + c*d（重複 番地数）の二つだけ。
               星の和が 0 のベクトルは k に依らず厳密に L^2 倍、
               星の上で一様なベクトルは L^2 + c*d 倍
      NG なら：固有値が二つでない

  検定F3  組 (A, Q) で閉じるか
      OK なら：Q' が A と Q の一次結合で書ける（係数は d, k, L だけで決まる）
      NG なら：閉じない。もう一つ量が要る

  検定F4  負の対照
      OK なら：c(k) を正しい値から 1 ずらすと検定F1 の等式が崩れる
      NG なら：検定F1 は NG を返せない検査

  検定F5  偶数次数でも L を上げれば内側が開くか
      d=6 は L=3 のとき根が 0 と 1 で内側に整数が無い。L=6 なら根は 0 と 2。
      OK なら：k=1, L=6 で c < 0 となり、A の倍率が L^2 = 36 を下回る
      NG なら：L を上げても内側は開かない
"""

from itertools import combinations
from fractions import Fraction
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


O = (0, 0)
verts = []
for sx in (1, -1):
    for sy in (1, -1):
        for sz in (1, -1):
            verts.append(((sx, 0), (sy, 0), (sz, 0)))
for s1 in (1, -1):
    for s2 in (1, -1):
        inv, ph = (-s1, s1), (0, s2)
        verts += [(O, inv, ph), (inv, ph, O), (ph, O, inv)]


def d2(p, q):
    s = O
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
        self.adj = adj
        self.k = k
        self.L = L
        self.d = len(adj[next(iter(adj))])
        self.c = self.d * k * k - 2 * L * k
        self.darcs = []
        self.index = {}
        for u in sorted(adj):
            for v in adj[u]:
                self.index[(u, v)] = len(self.darcs)
                self.darcs.append((u, v))
        self.N = len(self.darcs)

    def starsums(self, psi):
        return {u: sum(psi[self.index[(w, u)]] for w in self.adj[u]) for u in self.adj}

    def step(self, psi):
        S = self.starsums(psi)
        return [self.k * S[u] - self.L * psi[self.index[(v, u)]] for (u, v) in self.darcs]

    def A(self, psi):
        return sum(x * x for x in psi)

    def Q(self, psi):
        return sum(s * s for s in self.starsums(psi).values())


# ---------------- 検定F1 ----------------
cases = []
for name, adj in graphs.items():
    for k in range(-2, 5):
        for L in (1, 2, 3, 6):
            cases.append((name, adj, k, L))

f1_ok = True
checked = 0
for name, adj, k, L in cases:
    W = Walk(adj, k, L)
    for _ in range(4):
        psi = [random.randint(-4, 4) for _ in range(W.N)]
        lhs = W.A(W.step(psi))
        rhs = L * L * W.A(psi) + W.c * W.Q(psi)
        checked += 1
        if lhs != rhs:
            f1_ok = False
report("検定F1 二乗和の一歩", f1_ok,
       f"A' = L^2 A + c(k) Q を {checked} 件で厳密確認（3グラフ × k=-2..4 × L=1,2,3,6）")


# ---------------- 検定F2 ----------------
def eig_counts(W):
    """T^T T = L^2 I + c G。G は番地ごとの全1ブロック（大きさ d）"""
    # 星の和が 0 のベクトル（各番地で差を取る）
    zero_ok = True
    for _ in range(20):
        psi = [0] * W.N
        for u in W.adj:
            arcs = [W.index[(w, u)] for w in W.adj[u]]
            coef = [random.randint(-3, 3) for _ in arcs]
            coef[-1] = -sum(coef[:-1])
            for a, cc in zip(arcs, coef):
                psi[a] += cc
        if W.A(W.step(psi)) != W.L * W.L * W.A(psi):
            zero_ok = False
    # 星の上で一様なベクトル
    uni_ok = True
    for _ in range(20):
        psi = [0] * W.N
        for u in W.adj:
            cc = random.randint(-3, 3)
            for w in W.adj[u]:
                psi[W.index[(w, u)]] += cc
        if W.A(W.step(psi)) != (W.L * W.L + W.c * W.d) * W.A(psi):
            uni_ok = False
    return zero_ok, uni_ok


f2_ok = True
f2_note = []
for name, adj in graphs.items():
    for k, L in ((1, 3), (2, 3), (3, 3), (1, 6)):
        W = Walk(adj, k, L)
        z, u = eig_counts(W)
        f2_ok = f2_ok and z and u
    W = Walk(adj, 1, 3)
    f2_note.append(f"{name}: L^2 の重複 {W.N - len(adj)} / L^2+cd の重複 {len(adj)}")
report("検定F2 固有空間の分解", f2_ok, " | ".join(f2_note))


# ---------------- 検定F3 ----------------
f3_rows = []
f3_ok = True
for name, adj in graphs.items():
    for k, L in ((1, 3), (2, 3), (1, 6)):
        W = Walk(adj, k, L)
        samples = []
        for _ in range(6):
            psi = [random.randint(-4, 4) for _ in range(W.N)]
            p2 = W.step(psi)
            samples.append((W.A(psi), W.Q(psi), W.Q(p2)))
        # Q' = alpha*A + beta*Q を 2件で解き、残り4件で照合
        (a1, q1, r1), (a2, q2, r2) = samples[0], samples[1]
        det = a1 * q2 - a2 * q1
        if det == 0:
            continue
        alpha = Fraction(r1 * q2 - r2 * q1, det)
        beta = Fraction(a1 * r2 - a2 * r1, det)
        good = all(alpha * a + beta * q == r for a, q, r in samples[2:])
        f3_rows.append((name, k, L, alpha, beta, good))
        if not good:
            f3_ok = False

print()
print("  グラフ        k  L   Q' = alpha*A + beta*Q       閉じるか")
for name, k, L, al, be, good in f3_rows:
    print(f"  {name:<12} {k}  {L}   alpha={str(al):<8} beta={str(be):<10} {'はい' if good else 'いいえ'}")
report("検定F3 組 (A, Q) で閉じるか", f3_ok,
       "閉じない場合はもう一つ量が要る" if not f3_ok else "一次結合で閉じる")


# ---------------- 検定F4 ----------------
W = Walk(graphs["d=6 五芒星"], 3, 3)
broken = 0
for _ in range(20):
    psi = [random.randint(-4, 4) for _ in range(W.N)]
    lhs = W.A(W.step(psi))
    if lhs != W.L * W.L * W.A(psi) + (W.c + 1) * W.Q(psi):
        broken += 1
report("検定F4 負の対照", broken == 20,
       f"c を 1 ずらすと 20件中 {broken} 件で等式が崩れる")


# ---------------- 検定F5 ----------------
print()
print("  d=6 五芒星。L を上げると根が動く")
f5_rows = []
for L in (3, 6, 9, 12):
    roots = (0, Fraction(2 * L, 6))
    ks = [k for k in range(0, 6) if 0 < k < roots[1]]
    f5_rows.append((L, roots[1], ks, [6 * k * k - 2 * L * k for k in ks]))
    print(f"  L={L:<3} 根 0 と {roots[1]}   内側の整数 k={ks if ks else 'なし'}   c(k)={[6*k*k-2*L*k for k in ks]}")

W = Walk(graphs["d=6 五芒星"], 1, 6)
psi = [1 if t == 0 else 0 for t in range(W.N)]
prev = W.A(psi)
ratios = []
for _ in range(8):
    psi = W.step(psi)
    cur = W.A(psi)
    ratios.append(cur / prev)
    prev = cur
f5_ok = (W.c < 0 and all(r < 36 for r in ratios[:4]))
print("  k=1, L=6 （内側、c = " + str(W.c) + "）の倍率 : " + ", ".join(f"{r:.3f}" for r in ratios))
report("検定F5 偶数次数の内側", f5_ok, f"L=6 で根は 0 と 2、k=1 が内側に入る。L^2=36 を下回る")

print()
n_ok = sum(1 for _, ok, _ in results if ok)
print(f"{n_ok}/{len(results)} OK")
