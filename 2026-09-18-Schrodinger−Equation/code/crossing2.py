"""
crossing2.py --- k を動かしたとき T^T T がどうなるか

一歩  psi'(u->v) = k * sum_{w~u} psi(w->u) - L * psi(v->u)
の k を、2L/d に固定せず独立に動かす。

事前登録（測る前に書く）

  検定E4  T^T T の成分を k の関数として読む
      OK なら：T^T T = L^2 I + c(k) * G （G は同じ星に属する有向辺の対で 1、他は 0）
               で、c(k) = d*k^2 - 2*L*k。根は k = 0 と k = 2L/d のちょうど二つ
      NG なら：二次でない、あるいは根が二つでない

  検定E5  c(k) の符号
      OK なら：二根の内側で c < 0（二乗和が縮む）、外側で c > 0（伸びる）、
               根の上でちょうど L^2 倍。楕円／放物線／双曲線と同じ三領域になる
      NG なら：符号が三領域に分かれない

  検定E6  内側が整数で届く d があるか
      d = 6, L = 3 では根が 0 と 1 で隣り合い、内側に整数が無い。
      d = 3, L = 3 では根が 0 と 2 で、k = 1 が内側に入る。
      OK なら：d = 3 の k = 1 で二乗和が毎歩縮み、倍率が L^2 未満になる
      NG なら：内側の整数点が存在しないか、倍率が一定でない

  検定E7（負の対照）
      c(k) = 0 を満たさない k で T^T T が非スカラーであること（非対角に 0 以外が出る）
      OK なら：検定E4 は NG を返せる検査
      NG なら：どの k でもスカラーになり、検査になっていない
"""

from itertools import combinations

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


def arcs_of(adj):
    darcs, index = [], {}
    for u in sorted(adj):
        for v in adj[u]:
            index[(u, v)] = len(darcs)
            darcs.append((u, v))
    return darcs, index


def gram_full(adj, k, L):
    darcs, index = arcs_of(adj)
    N = len(darcs)

    def step(psi):
        insum = {u: sum(psi[index[(w, u)]] for w in adj[u]) for u in adj}
        return [k * insum[u] - L * psi[index[(v, u)]] for (u, v) in darcs]

    images = []
    for t in range(N):
        e = [0] * N
        e[t] = 1
        images.append(step(e))
    M = [[sum(x * y for x, y in zip(images[a], images[b])) for b in range(N)]
         for a in range(N)]
    return darcs, index, M


def Lmin(d):
    return d // 2 if d % 2 == 0 else d


# ---------------- 検定E4 ----------------
adj = graphs["d=6 五芒星"]
d, L = 6, 3
darcs, index, _ = gram_full(adj, 1, L)
N = len(darcs)

# G: 同じ頭 u を共有する有向辺の対（w->u と v->u が同じ星）
def same_star(a, b):
    return darcs[a][1] == darcs[b][1]


rows = []
form_ok = True
for k in range(-2, 5):
    _, _, M = gram_full(adj, k, L)
    c_pred = d * k * k - 2 * L * k
    ok = True
    for a in range(N):
        for b in range(N):
            want = (L * L if a == b else 0) + (c_pred if same_star(a, b) else 0)
            if M[a][b] != want:
                ok = False
    form_ok = form_ok and ok
    offvals = sorted({M[a][b] for a in range(N) for b in range(N) if a != b})
    rows.append((k, M[0][0], c_pred, offvals, ok))

print("   k   対角    c(k)=6k^2-6k   非対角の値      式と一致")
for k, diag, c, offvals, ok in rows:
    print(f"  {k:>2}   {diag:>4}    {c:>6}        {offvals}      {'はい' if ok else 'いいえ'}")

roots = [k for k, _, c, _, _ in rows if c == 0]
e4_ok = form_ok and roots == [0, 1]
report("検定E4 二次の形", e4_ok,
       f"T^T T = L^2 I + c(k) G / c(k) = d k^2 - 2 L k / 根 k = {roots}")


# ---------------- 検定E5 ----------------
def c_of(d, L, k):
    return d * k * k - 2 * L * k


sign_rows = []
for name, g in graphs.items():
    dd = len(g[next(iter(g))])
    LL = Lmin(dd)
    r2 = 2 * LL / dd
    inside = [k for k in range(-3, 6) if 0 < k < r2]
    outside = [k for k in range(-3, 6) if k < 0 or k > r2]
    sign_rows.append((name, dd, LL, r2,
                      [c_of(dd, LL, k) for k in inside],
                      [c_of(dd, LL, k) for k in outside]))

e5_ok = all(all(c < 0 for c in ins) and all(c > 0 for c in outs)
            for _, _, _, _, ins, outs in sign_rows)
print()
print("  グラフ        d  L  根 0 と      内側の c(k)    外側の c(k)")
for name, dd, LL, r2, ins, outs in sign_rows:
    print(f"  {name:<12} {dd}  {LL}  {r2:<10}  {ins if ins else '整数なし'}   {outs[:4]}…")
report("検定E5 符号の三領域", e5_ok, "内側 c<0（縮む） / 根の上 c=0（L^2 倍ちょうど） / 外側 c>0（伸びる）")


# ---------------- 検定E6 ----------------
def run_sumsq(adj, k, L, steps, init=None):
    darcs, index = arcs_of(adj)
    N = len(darcs)
    psi = init[:] if init else [1 if t == 0 else 0 for t in range(N)]
    out = []
    for _ in range(steps):
        insum = {u: sum(psi[index[(w, u)]] for w in adj[u]) for u in adj}
        psi = [k * insum[u] - L * psi[index[(v, u)]] for (u, v) in darcs]
        out.append(sum(x * x for x in psi))
    return out


g3 = graphs["d=3 稜"]
seq_in = run_sumsq(g3, 1, 3, 8)              # 内側（k=1、根は 0 と 2）
seq_root = run_sumsq(g3, 2, 3, 8)            # 根の上（k=2）
seq_out = run_sumsq(g3, 3, 3, 8)             # 外側（k=3）

prev = 1
ratios_in = []
for s in seq_in:
    ratios_in.append(s / prev)
    prev = s
prev = 1
ratios_root = []
for s in seq_root:
    ratios_root.append(s / prev)
    prev = s
prev = 1
ratios_out = []
for s in seq_out:
    ratios_out.append(s / prev)
    prev = s

e6_ok = (all(abs(r - 9) < 1e-9 for r in ratios_root)
         and any(r < 9 for r in ratios_in)
         and all(r > 9 for r in ratios_out[2:]))
print()
print("  d=3, L=3。根は k=0 と k=2")
print("  内側 k=1 の倍率 : " + ", ".join(f"{r:.3f}" for r in ratios_in))
print("  根上 k=2 の倍率 : " + ", ".join(f"{r:.3f}" for r in ratios_root))
print("  外側 k=3 の倍率 : " + ", ".join(f"{r:.3f}" for r in ratios_out))
report("検定E6 内側の整数点", e6_ok, "k=1 は L^2=9 を下回り、k=3 は上回る")


# ---------------- 検定E7 ----------------
nonscalar = []
for k in range(-2, 5):
    if c_of(6, 3, k) == 0:
        continue
    _, _, M = gram_full(graphs["d=6 五芒星"], k, 3)
    off = {M[a][b] for a in range(len(M)) for b in range(len(M)) if a != b}
    nonscalar.append(off != {0})
e7_ok = all(nonscalar)
report("検定E7 負の対照", e7_ok,
       f"根でない k ({len(nonscalar)}件) はすべて非スカラー ── 検定E4 は NG を返せる")

print()
n_ok = sum(1 for _, ok, _ in results if ok)
print(f"{n_ok}/{len(results)} OK")
