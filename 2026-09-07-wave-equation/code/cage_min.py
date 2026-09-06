"""1環の一歩 T を組み、剰余 m での最小多項式（クリロフ）と最小周期を出す。整数のみ。"""
import sys
from ring1 import build_ring1

def carrier():
    V, E = build_ring1()
    idx = {v: i for i, v in enumerate(V)}
    nb = {i: [] for i in range(len(V))}
    for a, b in E:
        nb[idx[a]].append(idx[b]); nb[idx[b]].append(idx[a])
    darts = []                       # (u,v) = u→v
    for v in nb:
        for u in nb[v]:
            darts.append((u, v))
    di = {d: i for i, d in enumerate(darts)}
    return nb, darts, di

NB, DARTS, DI = carrier()
E = len(DARTS)
C = {2: 12, 3: 8, 4: 6}

def step(psi, m):
    S = [0]*len(NB)
    for v in NB:
        S[v] = sum(psi[DI[(u, v)]] for u in NB[v]) % m
    out = [0]*E
    for i, (u, v) in enumerate(DARTS):        # ψ(u→v) → ψ'(v→u)
        out[DI[(v, u)]] = (C[len(NB[v])]*S[v] - 12*psi[i]) % m
    return out

def minpoly(psi0, m):
    """クリロフ列の最初の一次従属から最小多項式 g（係数リスト、低次から、モニック）"""
    rows, piv, combos = [], [], []
    psi = psi0[:]
    k = 0
    while True:
        v = psi[:]; comb = [0]*(k+1); comb[k] = 1
        for r, p, cb in zip(rows, piv, combos):
            f = v[p] % m
            if f:
                v = [(a - f*b) % m for a, b in zip(v, r)]
                comb = [(a - f*b) % m for a, b in zip(comb, cb + [0]*(k+1-len(cb)))]
        nz = next((j for j, a in enumerate(v) if a), None)
        if nz is None:
            g = comb
            inv = pow(g[-1], -1, m)
            return [a*inv % m for a in g]
        inv = pow(v[nz], -1, m)
        rows.append([a*inv % m for a in v]); piv.append(nz)
        combos.append([a*inv % m for a in comb])
        psi = step(psi, m); k += 1

def polymulmod(a, b, g, m):
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i+j] = (r[i+j] + x*y) % m
    d = len(g)-1
    for i in range(len(r)-1, d-1, -1):
        c = r[i]
        if c:
            r[i] = 0
            for j in range(d):
                r[i-d+j] = (r[i-d+j] - c*g[j]) % m
    return (r[:d] + [0]*d)[:d]

def polypow(a, e, g, m):
    r = [1]+[0]*(len(g)-2)
    while e:
        if e & 1: r = polymulmod(r, a, g, m)
        a = polymulmod(a, a, g, m); e >>= 1
    return r

def order_of_x(g, m):
    """x^N = 1 となる最小 N。N | m^d - 1 の倍数×m冪 の形を素因数から降ろす"""
    import sympy
    d = len(g)-1
    x = ([-g[0] % m] if d == 1 else [0, 1]+[0]*(d-2))
    # 候補上限：|(F_m[x]/g)^*| を割る N を探す。まず N0 = 素因数分解可能な倍数を作る
    # 単純に：m^d - 1 の倍数 × m^ceil(log_m d) を上限にする
    cap = (m**d - 1) * m**(max(1, d)).bit_length()
    N = (m**d - 1) * m**( (d-1).bit_length() )
    assert polypow(x, N, g, m) == ([1]+[0]*(d-1)), "x^N != 1"
    for p, e in sympy.factorint(N).items():
        for _ in range(e):
            if N % p: break
            if polypow(x, N//p, g, m) == ([1]+[0]*(d-1)): N //= p
            else: break
    return N

if __name__ == "__main__":
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    psi0 = [0]*E; psi0[0] = 1
    g = minpoly(psi0, m)
    print(f"m={m}  有向辺={E}  クリロフ次元 deg g = {len(g)-1}")
    print(f"最小周期 = {order_of_x(g, m)}")
