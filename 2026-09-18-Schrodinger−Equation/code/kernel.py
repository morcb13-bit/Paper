"""
kernel.py --- (S, O) に落ちる前後で何が捨てられているか

対象を一行で書く（判別法9）：
  いま見ているのは「写像 Phi: psi -> (S, O) の階数と核」。
  120本の辺の上の psi と、番地あたり2個・計40個の (S, O) の関係。

事前登録（判別法6）

  検定K1  階数
      sum_u S_u = sum_u O_u = sum psi なので階数は 39 以下のはず
      OK なら：階数がちょうど 39、核の次元が 81
      NG なら：別の値。その場合 (S, O) の独立な個数を言い直す

  検定K2  核が一歩で保たれるか
      核 = S も O も全番地で 0 の psi
      OK なら：核の psi は一歩のあとも核に留まる（部分空間が保たれる）
      NG なら：核は一歩で混ざる

  検定K3  核の上での一歩の形
      S = 0 なら psi'(u->v) = -L psi(v->u) で、向きを裏返して L 倍するだけのはず
      OK なら：核の上で T^2 = L^2 I（二歩で元に戻る。倍率のほかに何も起きない）
      NG なら：核の上でも別の動きがある

  検定K4  (S, O) が状態か射影か
      OK なら：核が 0 でない（＝(S, O) は psi を一意に決めない。射影である）。
               かつ像は sum S = sum O を満たす対で尽きる
      NG なら：全単射。その場合 (S, O) は状態と呼べる

  検定K5  負の対照
      OK なら：核でない psi を入れると K3 の等式が崩れる
      NG なら：検査になっていない
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


adj = {i: [] for i in range(20)}
for i, j in combinations(range(20), 2):
    if d2(verts[i], verts[j]) == (4, 0):
        adj[i].append(j)
        adj[j].append(i)

k, L, d = 1, 3, 6
darcs, index = [], {}
for u in sorted(adj):
    for v in adj[u]:
        index[(u, v)] = len(darcs)
        darcs.append((u, v))
N = len(darcs)


def Svec(psi):
    return [sum(psi[index[(w, u)]] for w in adj[u]) for u in sorted(adj)]


def Ovec(psi):
    return [sum(psi[index[(u, w)]] for w in adj[u]) for u in sorted(adj)]


def step(psi):
    S = Svec(psi)
    return [k * S[u] - L * psi[index[(v, u)]] for (u, v) in darcs]


# Phi の行列（40 x 120）
Phi = []
for u in sorted(adj):
    row = [0] * N
    for w in adj[u]:
        row[index[(w, u)]] = 1
    Phi.append(row)
for u in sorted(adj):
    row = [0] * N
    for w in adj[u]:
        row[index[(u, w)]] = 1
    Phi.append(row)


def rank_and_kernel(M, ncols):
    """有理数で掃き出し、階数と核の基底を返す"""
    A = [[Fraction(x) for x in row] for row in M]
    nrows = len(A)
    pivots = []
    r = 0
    for c in range(ncols):
        piv = None
        for i in range(r, nrows):
            if A[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(nrows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        pivots.append(c)
        r += 1
        if r == nrows:
            break
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free:
        v = [Fraction(0)] * ncols
        v[fc] = Fraction(1)
        for i, pc in enumerate(pivots):
            v[pc] = -A[i][fc]
        basis.append(v)
    return r, basis


rank, kbasis = rank_and_kernel(Phi, N)
report("検定K1 階数", rank == 39 and len(kbasis) == 81,
       f"階数 {rank} / 核の次元 {len(kbasis)}（有向辺 {N}本、番地 {len(adj)}個）")


# ---------------- 検定K2・K3 ----------------
def to_int_vec(v):
    den = 1
    for x in v:
        den = den * x.denominator // __import__("math").gcd(den, x.denominator)
    return [int(x * den) for x in v]


kernel_ints = [to_int_vec(v) for v in kbasis]
k2_ok = True
k3_ok = True
for _ in range(30):
    psi = [0] * N
    for _ in range(4):
        b = random.choice(kernel_ints)
        cc = random.randint(-3, 3)
        psi = [a + cc * x for a, x in zip(psi, b)]
    if any(x != 0 for x in Svec(psi)) or any(x != 0 for x in Ovec(psi)):
        k2_ok = False
    p1 = step(psi)
    if any(x != 0 for x in Svec(p1)) or any(x != 0 for x in Ovec(p1)):
        k2_ok = False
    # 一歩は「裏返して -L 倍」か
    rev = [-L * psi[index[(v, u)]] for (u, v) in darcs]
    if p1 != rev:
        k3_ok = False
    # 二歩で L^2 倍の元の姿か
    p2 = step(p1)
    if p2 != [L * L * x for x in psi]:
        k3_ok = False
report("検定K2 核が保たれる", k2_ok, "核の psi は一歩のあとも核に留まる")
report("検定K3 核の上での一歩", k3_ok, "psi' = -L * (向きを裏返したもの) / T^2 = L^2 I")


# ---------------- 検定K4 ----------------
img_ok = True
for _ in range(20):
    psi = [random.randint(-5, 5) for _ in range(N)]
    if sum(Svec(psi)) != sum(Ovec(psi)) or sum(Svec(psi)) != sum(psi):
        img_ok = False
k4_ok = (len(kbasis) > 0 and img_ok)
report("検定K4 状態か射影か", k4_ok,
       f"核が {len(kbasis)} 次元あるので (S, O) は psi を一意に決めない ── 射影。"
       f"像は sum S = sum O = sum psi を満たす対")


# ---------------- 検定K5 ----------------
broke = 0
for _ in range(20):
    psi = [random.randint(-4, 4) for _ in range(N)]
    if any(x != 0 for x in Svec(psi)):
        rev = [-L * psi[index[(v, u)]] for (u, v) in darcs]
        if step(psi) != rev:
            broke += 1
report("検定K5 負の対照", broke == 20,
       f"核でない psi では 20件中 {broke} 件で K3 の形が崩れる")


# ---------------- 参考 ----------------
print()
print("  次元の内訳")
print(f"    有向辺           {N}")
print(f"    核（S=0 かつ O=0） {len(kbasis)}   一歩は裏返しと -L 倍のみ")
print(f"    (S, O) が担う     {rank}   二階漸化式が走る側")
print(f"    合計             {len(kbasis) + rank}")
print()
print("  F2 の 100次元（S=0 のみ）との関係: 120 - 20 - 20 + 1 = 81")

n_ok = sum(1 for _, ok, _ in results if ok)
print()
print(f"{n_ok}/{len(results)} OK")
