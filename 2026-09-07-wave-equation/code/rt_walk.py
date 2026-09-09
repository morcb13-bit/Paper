#  rt_walk.py  ── 黄金のひし形30面体（RT）を担体にしたオートマトン
#
#  すべて Z[phi] の厳密整数。浮動小数は最後の表示に一度も使わない。
#
#  検定1 RT の組み立て
#      OK なら：頂点32（次数5が12・次数3が20）・辺60・面30、面は全部ひし形で対角比 phi
#      NG なら：担体が RT になっていない。以降の数値は全部無効
#
#  検定2 面の隣接（30面を担体にしたときの網）
#      OK なら：次数4一様・辺60。一歩の係数は c*d = 2t で c=1,t=2 が最小
#      NG なら：面側を担体にできない。頂点側（32・次数3と5）へ切り替える
#
#  検定3 一歩が振幅を保つか
#      OK なら：T^T T = t^2 I。二乗和は毎歩ちょうど t^2 倍（1環の144倍にあたる量）
#      NG なら：一歩の書き方が間違っている
#
#  検定4 帯（zone）
#      OK なら：帯は6本、各10面、帯の中だけで見ると次数2一様の輪
#      NG なら：Kabai の帯の読みが面側の網に写っていない
#
#  検定5 辺の側 → 面の側の還元（1環の §1-1 にあたる恒等式）
#      OK なら：det(xI-T) = (x^2-t^2)^(E-V) * prod(x^2 - nu x + t^2)、nu は面の網の隣接固有値
#      NG なら：120次元を直接扱うしかない
#      検査は落ちうる：べき和（トレース）を k=1..20 で突き合わせる。1つでも違えば NG
#
#  検定6 檻（T^k = I の最小 k）
#      OK なら：剰余 m ごとに k が出る
#      負の対照：隣接を1本だけ付け替えた網（次数が4でなくなる）では c*d=2t が整数で書けず落ちる
#
#  検定7 判別式に 5 が出るか（引継書 v247 §4-3 の当の問い）
#      5 が出る なら：phi が時間の側に乗る。1環は 17 と 2 で出なかった
#      出ない なら：RT でも phi は出ない

from fractions import Fraction
from itertools import permutations

# ---------------------------------------------------------------- Z[phi]
# x = (a, b)  は a + b*phi、phi^2 = phi + 1

def padd(x, y): return (x[0] + y[0], x[1] + y[1])
def psub(x, y): return (x[0] - y[0], x[1] - y[1])
def pmul(x, y):
    a, b = x; c, d = y
    return (a * c + b * d, a * d + b * c + b * d)
def pneg(x): return (-x[0], -x[1])

PZ = (0, 0); PONE = (1, 0); PHI = (0, 1)

def psign(x):
    # 2*(a+b*phi) = (2a+b) + b*sqrt5
    a, b = x
    u = 2 * a + b
    if b == 0: return (u > 0) - (u < 0)
    if u == 0: return (b > 0) - (b < 0)
    if u > 0 and b > 0: return 1
    if u < 0 and b < 0: return -1
    # 符号が違う: u^2 vs 5b^2
    lhs = u * u; rhs = 5 * b * b
    if lhs == rhs: return 0
    big_u = lhs > rhs
    return 1 if (u > 0) == big_u else -1

def pcmp(x, y): return psign(psub(x, y))

# 3次元ベクトル
def vsub(u, v): return tuple(psub(a, b) for a, b in zip(u, v))
def vadd(u, v): return tuple(padd(a, b) for a, b in zip(u, v))
def vdot(u, v):
    s = PZ
    for a, b in zip(u, v): s = padd(s, pmul(a, b))
    return s
def vnorm2(u): return vdot(u, u)

# ---------------------------------------------------------------- 担体 RT
def cyc(t):
    a, b, c = t
    return [(a, b, c), (c, a, b), (b, c, a)]

def build_rt():
    ico = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            base = (PZ, (0, s1), (s2, 0))          # (0, ±phi, ±1)
            ico += cyc(base)
    ico = sorted(set(ico))
    dod = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            for s3 in (1, -1):
                dod.append(((s1, 0), (s2, 0), (s3, 0)))   # (±1,±1,±1)
    PHIM1 = (-1, 1)                                        # phi-1 = 1/phi
    for s1 in (1, -1):
        for s2 in (1, -1):
            base = (PZ, (-s1, s1), (0, s2))                # (0, ±1/phi, ±phi)
            dod += cyc(base)
    dod = sorted(set(dod))
    return ico, dod

def rt_graph():
    ico, dod = build_rt()
    verts = ico + dod
    n_ico = len(ico)
    # ico-dod の最小距離
    best = None
    for i, u in enumerate(ico):
        for j, v in enumerate(dod):
            d = vnorm2(vsub(u, v))
            if best is None or pcmp(d, best) < 0: best = d
    edges = set()
    for i, u in enumerate(ico):
        for j, v in enumerate(dod):
            if vnorm2(vsub(u, v)) == best:
                edges.add((i, n_ico + j))
    return verts, n_ico, sorted(edges), best

# ---------------------------------------------------------------- 面
def rt_faces(verts, n_ico, edges):
    adj = {i: set() for i in range(len(verts))}
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    faces = []
    # 面 = ico 2個 + dod 2個 の4輪
    for i in range(n_ico):
        for k in range(i + 1, n_ico):
            common = sorted(adj[i] & adj[k])
            if len(common) == 2:
                faces.append((i, common[0], k, common[1]))
    return faces, adj

# ---------------------------------------------------------------- 一歩
def arc_operator(nbr, t):
    """有向辺の上の一歩。 psi'(u->v) = c * sum_{w~u} psi(w->u) - t * psi(v->u)
       c * deg(u) = 2t"""
    arcs = []
    for u in nbr:
        for v in nbr[u]:
            arcs.append((u, v))
    idx = {a: i for i, a in enumerate(arcs)}
    n = len(arcs)
    T = [[0] * n for _ in range(n)]
    for (u, v) in arcs:
        d = len(nbr[u])
        assert (2 * t) % d == 0, f"c が整数にならない: d={d}, t={t}"
        c = (2 * t) // d
        r = idx[(u, v)]
        for w in nbr[u]:
            T[r][idx[(w, u)]] += c
        T[r][idx[(v, u)]] -= t
    return T, arcs, idx

def mat_mul(A, B, mod=None):
    n = len(A); p = len(B[0]); m = len(B)
    C = [[0] * p for _ in range(n)]
    for i in range(n):
        Ai = A[i]; Ci = C[i]
        for k in range(m):
            a = Ai[k]
            if a:
                Bk = B[k]
                for j in range(p): Ci[j] += a * Bk[j]
        if mod:
            for j in range(p): Ci[j] %= mod
    return C

def mat_trace(A): return sum(A[i][i] for i in range(len(A)))

# ---------------------------------------------------------------- 多項式
def poly_mul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i + j] += x * y
    return r

def poly_trim(a):
    while len(a) > 1 and a[-1] == 0: a.pop()
    return a

def charpoly(A):
    """Leverrier-Faddeev。整数行列 → 係数リスト（低次から）"""
    n = len(A)
    I = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    M = [row[:] for row in I]
    cs = [Fraction(1)]
    Mk = [row[:] for row in I]
    coeffs = [Fraction(1)]
    Mcur = None
    c_list = [Fraction(1)]
    Mprev = I
    for k in range(1, n + 1):
        AM = mat_mul(A, Mprev)
        ck = Fraction(-mat_trace(AM), k)
        Mprev = [[AM[i][j] + (ck if i == j else 0) for j in range(n)] for i in range(n)]
        c_list.append(ck)
    # char = x^n + c1 x^{n-1} + ... + cn
    out = [int(c) for c in reversed(c_list)]
    return out   # 低次から

def poly_divmod(a, b):
    a = a[:]; db = len(b) - 1
    q = [0] * max(1, len(a) - db)
    while len(a) - 1 >= db and any(a):
        d = len(a) - 1 - db
        f = Fraction(a[-1], b[-1])
        if f.denominator != 1: return None, None
        f = int(f)
        q[d] = f
        for i, bi in enumerate(b): a[i + d] -= f * bi
        poly_trim(a)
        if len(a) == 1 and a[0] == 0: break
    return q, a

def factor_small(p):
    """根の絶対値が小さいモニック整数多項式を 1次・2次に分解する"""
    facs = []
    cur = p[:]
    for r in range(-8, 9):
        while len(cur) > 1:
            q, rem = poly_divmod(cur, [-r, 1])
            if q is not None and len(rem) == 1 and rem[0] == 0:
                facs.append([-r, 1]); cur = q
            else:
                break
    for b in range(-16, 17):
        for c in range(-40, 41):
            while len(cur) > 2:
                q, rem = poly_divmod(cur, [c, b, 1])
                if q is not None and len(rem) == 1 and rem[0] == 0:
                    facs.append([c, b, 1]); cur = q
                else:
                    break
    return facs, cur

def pstr(p):
    terms = []
    for i in range(len(p) - 1, -1, -1):
        c = p[i]
        if c == 0: continue
        if i == 0: terms.append(f"{c:+d}")
        elif i == 1: terms.append(f"{c:+d}x" if abs(c) != 1 else ("+x" if c > 0 else "-x"))
        else: terms.append((f"{c:+d}" if abs(c) != 1 else ("+" if c > 0 else "-")) + f"x^{i}")
    return "".join(terms).lstrip("+")

# ---------------------------------------------------------------- 剰余体での位数
def poly_mulmod(a, b, f, m):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i + j] = (r[i + j] + x * y) % m
    # f でわる
    df = len(f) - 1
    inv = pow(f[-1], -1, m)
    for i in range(len(r) - 1, df - 1, -1):
        if r[i]:
            fac = r[i] * inv % m
            for j, fj in enumerate(f): r[i - df + j] = (r[i - df + j] - fac * fj) % m
    r = r[:df]
    while len(r) < df: r.append(0)
    return r

def poly_powmod(a, e, f, m):
    r = [1] + [0] * (len(f) - 2)
    while e:
        if e & 1: r = poly_mulmod(r, a, f, m)
        a = poly_mulmod(a, a, f, m); e >>= 1
    return r

def factorize(n):
    f = {}; d = 2
    while d * d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def is_irred_modm(f, m):
    d = len(f) - 1
    if d <= 1: return True
    x = [0, 1] + [0] * (d - 2)
    # x^(m^i) - x の gcd を素朴に：1次・2次の因数の有無を総当り
    for r in range(m):
        v = 0
        for c in reversed(f): v = (v * r + c) % m
        if v == 0: return False
    if d <= 3: return True
    for b in range(m):
        for c in range(m):
            g = [c, b, 1]
            q, rem = poly_divmod_mod(f, g, m)
            if all(t == 0 for t in rem): return False
    return True

def poly_divmod_mod(a, b, m):
    a = [x % m for x in a]; db = len(b) - 1
    inv = pow(b[-1], -1, m)
    q = [0] * max(1, len(a) - db)
    for i in range(len(a) - 1, db - 1, -1):
        if a[i]:
            f = a[i] * inv % m
            q[i - db] = f
            for j, bj in enumerate(b): a[i - db + j] = (a[i - db + j] - f * bj) % m
    return q, a[:db]

def factor_modm(f, m):
    """小さい次数の多項式を F_m 上で完全に因数分解（総当り）"""
    f = [x % m for x in f]
    out = []
    stack = [f]
    while stack:
        g = stack.pop()
        if len(g) - 1 == 0: continue
        if len(g) - 1 == 1: out.append(g); continue
        split = False
        for r in range(m):
            v = 0
            for c in reversed(g): v = (v * r + c) % m
            if v == 0:
                q, _ = poly_divmod_mod(g, [(-r) % m, 1], m)
                stack.append([(-r) % m, 1]); stack.append(q); split = True; break
        if split: continue
        if len(g) - 1 == 2 or len(g) - 1 == 3:
            out.append(g); continue
        for b in range(m):
            for c in range(m):
                cand = [c, b, 1]
                q, rem = poly_divmod_mod(g, cand, m)
                if all(t == 0 for t in rem):
                    stack.append(cand); stack.append(q); split = True; break
            if split: break
        if not split: out.append(g)
    return out

def order_of_x(f, m):
    d = len(f) - 1
    N = m ** d - 1
    if N == 0: return 1
    k = N
    for p, e in factorize(N).items():
        for _ in range(e):
            if k % p: break
            t = k // p
            if poly_powmod([0, 1] + [0] * (d - 2) if d >= 2 else [0], t, f, m) == ([1] + [0] * (d - 1)):
                k = t
            else: break
    return k

def cage(charp, m):
    """T^k = I の最小 k（mod m）。charp は T の特性多項式（低次から）"""
    facs = factor_modm(charp, m)
    seen = {}
    for g in facs:
        key = tuple(g)
        seen[key] = seen.get(key, 0) + 1
    k = 1
    for g, mult in seen.items():
        g = list(g)
        if g == [0, 1]:   # x | charp → 可逆でない
            return None
        o = order_of_x(g, m)
        k = k * o // __import__("math").gcd(k, o)
        if mult > 1:
            e = 1
            while m ** e < mult: e += 1
            k = k * (m ** e) // __import__("math").gcd(k, m ** e)
    return k

# ================================================================ 実行
verts, n_ico, edges, elen2 = rt_graph()
faces, vadj = rt_faces(verts, n_ico, edges)

deg = {i: len(vadj[i]) for i in range(len(verts))}
d5 = sum(1 for i in deg if deg[i] == 5); d3 = sum(1 for i in deg if deg[i] == 3)

# ひし形か・対角比
ok_rh = True; ratios = set()
for (a, b, c, d) in faces:
    P = [verts[a], verts[b], verts[c], verts[d]]
    sides = [vnorm2(vsub(P[(i + 1) % 4], P[i])) for i in range(4)]
    if len(set(sides)) != 1: ok_rh = False
    D1 = vnorm2(vsub(P[2], P[0])); D2 = vnorm2(vsub(P[3], P[1]))
    ratios.add((D1, D2))

print("検定1 RT の組み立て")
print(f"  頂点 {len(verts)}（次数5が{d5}・次数3が{d3}） 辺 {len(edges)} 面 {len(faces)}")
print(f"  全部ひし形か: {ok_rh}   辺の二乗長 {elen2}")
d1s = set(r[0] for r in ratios); d2s = set(r[1] for r in ratios)
print(f"  対角の二乗長 {sorted(d1s)} / {sorted(d2s)}")
print(f"  → 検定1 {'OK' if (len(verts)==32 and len(edges)==60 and len(faces)==30 and ok_rh and d5==12 and d3==20) else 'NG'}")

# ---- 面の隣接
fidx = {frozenset(f): i for i, f in enumerate(faces)}
fedges = {}
for i, f in enumerate(faces):
    for k in range(4):
        e = frozenset((f[k], f[(k + 1) % 4]))
        fedges.setdefault(e, []).append(i)
nbr = {i: [] for i in range(len(faces))}
n_fe = 0
for e, fs in fedges.items():
    if len(fs) == 2:
        nbr[fs[0]].append(fs[1]); nbr[fs[1]].append(fs[0]); n_fe += 1
fdeg = sorted(set(len(nbr[i]) for i in nbr))
print("\n検定2 面の隣接（担体：30面）")
print(f"  頂点 {len(faces)} 辺 {n_fe} 次数 {fdeg}")
print(f"  → 検定2 {'OK' if fdeg==[4] and n_fe==60 else 'NG'}")

# ---- 帯
def direction(u, v):
    d = vsub(verts[u], verts[v])
    # 符号の正規化
    for comp in d:
        s = psign(comp)
        if s: 
            if s < 0: d = tuple(pneg(c) for c in d)
            break
    return d
dirs = {}
for a, b in edges:
    dirs.setdefault(direction(a, b), []).append(frozenset((a, b)))
zones = []
for dvec, es in dirs.items():
    fs = set()
    for e in es: fs.update(fedges[e])
    zones.append(sorted(fs))
print("\n検定4 帯（zone）")
print(f"  辺の方向 {len(dirs)} 種、各 {sorted(set(len(v) for v in dirs.values()))} 本")
print(f"  帯 {len(zones)} 本、各 {sorted(set(len(z) for z in zones))} 面")
z0 = zones[0]
zsub = {f: [g for g in nbr[f] if g in set(z0)] for f in z0}
zdeg = sorted(set(len(zsub[f]) for f in zsub))
print(f"  帯の中だけで見た次数 {zdeg}")
print(f"  → 検定4 {'OK' if len(zones)==6 and all(len(z)==10 for z in zones) and zdeg==[2] else 'NG'}")

# ---- 一歩（RT 30面、t=2, c=1）
t = 2
T, arcs, aidx = arc_operator(nbr, t)
n = len(arcs)
TT = mat_mul([[T[j][i] for j in range(n)] for i in range(n)], T)
ok3 = all(TT[i][j] == (t * t if i == j else 0) for i in range(n) for j in range(n))
print(f"\n検定3 一歩が振幅を保つか（担体30面・t={t}・c={2*t//4}）")
print(f"  有向辺 {n}   T^T T = {t*t} I : {ok3}")
print(f"  → 検定3 {'OK' if ok3 else 'NG'}")

# ---- 面の側の隣接固有値
A = [[0] * len(faces) for _ in faces]
for i in nbr:
    for j in nbr[i]: A[i][j] = 1
cpA = charpoly(A)
facsA, restA = factor_small(cpA)
from collections import Counter
cntA = Counter(tuple(f) for f in facsA)
print("\n面の網（30頂点・次数4）の特性多項式")
for f, k in sorted(cntA.items(), key=lambda kv: (len(kv[0]), kv[0])):
    print(f"  {pstr(list(f)):<24} 重複 {k}")
print(f"  残り: {pstr(restA) if len(restA)>1 else '（無し）'}")

# ---- 恒等式で辺の側へ
def to_arc_factor(g, t):
    """x^d * g((x^2+t^2)/x)"""
    d = len(g) - 1
    out = [0] * (2 * d + 1)
    for j, c in enumerate(g):
        if not c: continue
        term = [1]
        for _ in range(j): term = poly_mul(term, [t * t, 0, 1])
        for _ in range(d - j): term = poly_mul(term, [0, 1])
        for i, v in enumerate(term): out[i] += c * v
    return poly_trim(out)

pred = [1]
for f, k in cntA.items():
    af = to_arc_factor(list(f), t)
    for _ in range(k): pred = poly_mul(pred, af)
extra = poly_mul([1], [1])
base = [-t * t, 0, 1]
for _ in range(n_fe - len(faces)): pred = poly_mul(pred, base)
pred = poly_trim(pred)

# べき和で突き合わせ
def power_sums_from_poly(p, K):
    # ニュートンの公式（モニック、低次から）
    nn = len(p) - 1
    a = [p[nn - i] for i in range(nn + 1)]   # a[0]=1, a[i]= coeff of x^{n-i}
    s = []
    for k in range(1, K + 1):
        v = -k * a[k] if k <= nn else 0
        for i in range(1, min(k - 1, nn) + 1): v -= a[i] * s[k - i - 1]
        s.append(v)
    return s
K = 20
ps_pred = power_sums_from_poly(pred, K)
Tk = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
ps_real = []
for k in range(K):
    Tk = mat_mul(Tk, T)
    ps_real.append(mat_trace(Tk))
ok5 = (len(pred) - 1 == n) and ps_pred == ps_real
print(f"\n検定5 辺の側 → 面の側の還元")
print(f"  予測次数 {len(pred)-1} / 実次数 {n}")
print(f"  べき和 k=1..{K} 一致: {ps_pred == ps_real}")
print(f"  → 検定5 {'OK' if ok5 else 'NG'}")

# ---- 檻
print("\n検定6 檻（T^k = I の最小 k）  担体：RT 30面")
rows = []
for m in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
    k = cage(pred, m)
    rows.append((m, k))
    print(f"  m={m:<3} k={k}")

# ---- 帯だけの担体（10面の輪）
z = zones[0]
zs = set(z)
nbrz = {f: [g for g in nbr[f] if g in zs] for f in z}
tz = 1
Tz, arcz, _ = arc_operator(nbrz, tz)
nz = len(arcz)
Az = [[0] * len(z) for _ in z]
pos = {f: i for i, f in enumerate(z)}
for f in nbrz:
    for g in nbrz[f]: Az[pos[f]][pos[g]] = 1
cpz = charpoly(Az)
predz = [1]
fz, restz = factor_small(cpz)
cz = Counter(tuple(x) for x in fz)
for f, k in cz.items():
    af = to_arc_factor(list(f), tz)
    for _ in range(k): predz = poly_mul(predz, af)
for _ in range(10 - 10): pass
predz = poly_trim(predz)
print(f"\n帯だけを担体にした場合（10面の輪・次数2一様・t={tz}, c={tz}）")
print(f"  有向辺 {nz}")
Tk = [[1 if i == j else 0 for j in range(nz)] for i in range(nz)]
k = 0
while True:
    Tk = mat_mul(Tk, Tz); k += 1
    if all(Tk[i][j] == (1 if i == j else 0) for i in range(nz) for j in range(nz)):
        break
    if k > 200: k = None; break
print(f"  整数のまま T^k = I になる最小 k : {k}")
perm = all(sum(1 for x in row if x != 0) == 1 for row in Tz) and all(abs(x) in (0, 1) for row in Tz for x in row)
print(f"  一歩が置換（各行に非零が1つ・値±1）か : {perm}")

# ---- 負の対照
print("\n負の対照（検査が落ちうるか）")
bad = {i: list(v) for i, v in nbr.items()}
bad[0] = bad[0][:3]
bad[bad[0][0]] = [x for x in bad[bad[0][0]]]
try:
    arc_operator(bad, 2)
    print("  次数3の点を混ぜた網で t=2 → 通ってしまった（検査が甘い）")
except AssertionError as e:
    print(f"  次数3の点を混ぜた網で t=2 → 落ちる: {e}")
