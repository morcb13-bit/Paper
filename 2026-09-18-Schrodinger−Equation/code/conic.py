"""
conic.py --- 対 (S_n, S_{n+1}) に乗る二次形式

対象を一行で書く（判別法9）：
  二階漸化式 s_{n+2} = lambda s_{n+1} - L^2 s_n の、隣り合う二項の対が
  平面のどこに乗るか。

事前登録（判別法6）

  検定J1  二次形式の倍率
      q(a, b) = b^2 - lambda a b + L^2 a^2 として
      OK なら：q(s_{n+1}, s_{n+2}) = L^2 * q(s_n, s_{n+1}) が整数のまま厳密成立
      NG なら：この形は乗らない。別の形を探す

  検定J2  D の符号と曲線の種類
      q の判別式は lambda^2 - 4 L^2 で、二階の D と同じ
      OK なら：D<0 で楕円、D=0 で平行二直線（退化）、D>0 で双曲線に分かれる。
               五芒星グラフでは lambda=6 が D=0、他の四つが D<0
      NG なら：分類がそうならない

  検定J3  グラフの側から出た数列で確かめる
      OK なら：実際に psi を走らせて S を取り、固有ベクトルへ射影した成分が
               J1 の等式を満たす
      NG なら：固有値ごとの分解が効いていない

  検定J4  負の対照
      OK なら：lambda を 1 ずらすと J1 が崩れる
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


L = 3
EIGS = [6, 2, 1, -2, -3]


def q(a, b, lam):
    return b * b - lam * a * b + L * L * a * a


# ---------------- 検定J1 ----------------
j1_ok = True
for lam in EIGS:
    for _ in range(30):
        s0, s1 = random.randint(-9, 9), random.randint(-9, 9)
        a, b = s0, s1
        for _ in range(10):
            a2, b2 = b, lam * b - L * L * a
            if q(a2, b2, lam) != L * L * q(a, b, lam):
                j1_ok = False
            a, b = a2, b2
report("検定J1 二次形式の倍率", j1_ok,
       "q(s_{n+1}, s_{n+2}) = L^2 q(s_n, s_{n+1}) を 5固有値 × 30初期値 × 10歩で厳密確認")


# ---------------- 検定J2 ----------------
print()
print("  lambda   D = lambda^2 - 4L^2   曲線")
kinds = {}
for lam in EIGS:
    D = lam * lam - 4 * L * L
    kind = "楕円" if D < 0 else ("平行二直線（退化）" if D == 0 else "双曲線")
    kinds[lam] = kind
    print(f"    {lam:>3}         {D:>4}           {kind}")
j2_ok = (kinds[6] == "平行二直線（退化）"
         and all(kinds[l] == "楕円" for l in EIGS if l != 6))
report("検定J2 D の符号と曲線", j2_ok, "五芒星グラフの5固有値の分類")

# 退化の形を明示：lambda=6 では q = (b - 3a)^2
deg_ok = all(q(a, b, 6) == (b - 3 * a) ** 2
             for a in range(-5, 6) for b in range(-5, 6))
report("検定J2b 退化の形", deg_ok, "lambda=6 では q(a,b) = (b - 3a)^2")


# ---------------- 検定J3 ----------------
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


def d2(p, qq):
    s = O0
    for i in range(3):
        t = zsub(p[i], qq[i])
        s = zadd(s, zmul(t, t))
    return s


adj = {i: [] for i in range(20)}
for i, j in combinations(range(20), 2):
    if d2(verts[i], verts[j]) == (4, 0):
        adj[i].append(j)
        adj[j].append(i)

darcs, index = [], {}
for u in sorted(adj):
    for v in adj[u]:
        index[(u, v)] = len(darcs)
        darcs.append((u, v))
N = len(darcs)


def Svec(psi):
    return [sum(psi[index[(w, u)]] for w in adj[u]) for u in sorted(adj)]


def step(psi):
    S = Svec(psi)
    return [S[u] - L * psi[index[(v, u)]] for (u, v) in darcs]


def adjmul(vec):
    return [sum(vec[w] for w in adj[u]) for u in sorted(adj)]


# 固有ベクトルへの射影を使わず、(A - mu I) を他の固有値ぶん掛けて成分を取り出す
def project(vec, lam):
    out = list(vec)
    for mu in EIGS:
        if mu == lam:
            continue
        out = [x - mu * y for x, y in zip(adjmul(out), out)]
    return out


psi = [random.randint(-3, 3) for _ in range(N)]
seq = [Svec(psi)]
for _ in range(8):
    psi = step(psi)
    seq.append(Svec(psi))

j3_ok = True
j3_note = []
for lam in EIGS:
    comps = [project(s, lam) for s in seq]
    nz = next((i for i in range(20) if any(c[i] != 0 for c in comps)), None)
    if nz is None:
        j3_note.append(f"λ={lam}:成分0")
        continue
    ser = [c[nz] for c in comps]
    ok = all(q(ser[n + 1], ser[n + 2], lam) == L * L * q(ser[n], ser[n + 1], lam)
             for n in range(len(ser) - 2))
    j3_ok = j3_ok and ok
    j3_note.append(f"λ={lam}:{'OK' if ok else 'NG'}")
report("検定J3 走らせた数列で確認", j3_ok, " ".join(j3_note))


# ---------------- 検定J4 ----------------
broke = 0
total = 0
skipped = []
for lam in EIGS:
    for _ in range(10):
        a, b = random.randint(-9, 9), random.randint(-9, 9)
        if (a, b) == (0, 0) or (lam == 6 and b == 3 * a):
            skipped.append((lam, a, b))   # lambda=6 の不変直線 b=3a は倍率 3 で回らない
            continue
        total += 1
        aa, bb = a, b
        fails = False
        for _ in range(5):                # 一歩では弱い。5歩のどこかで崩れるかを見る
            a2, b2 = bb, lam * bb - L * L * aa
            if q(a2, b2, lam + 1) != L * L * q(aa, bb, lam + 1):
                fails = True
            aa, bb = a2, b2
        if fails:
            broke += 1
report("検定J4 負の対照", broke == total,
       f"lambda を 1 ずらすと {total}件中 {broke} 件で崩れる（5歩。除外した退化点 {len(skipped)}件）")


# ---------------- 図のためのデータ ----------------
print()
print("  図2 用の軌道（L^n で割る前の整数。s0=1, s1=0 から8歩）")
for lam in EIGS:
    a, b = 1, 0
    ser = [a, b]
    for _ in range(8):
        a, b = b, lam * b - L * L * a
        ser.append(b)
    print(f"    λ={lam:>3}  q={q(1, 0, lam):>3}  {ser}")

print()
n_ok = sum(1 for _, ok, _ in results if ok)
print(f"{n_ok}/{len(results)} OK")
