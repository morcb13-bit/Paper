# b13_iconet.py — ペンローズの床の上で正二十面体の向きを繋ぎ、経路の和で推論し、加算だけで学ぶ
# 判定はすべて Z[φ] の整数で行う。浮動小数は担体（床）を作るときの線の交点探しにだけ使い、
# 番地は整数の 5 つ組で持つ。
import itertools, math, random, sys
sys.setrecursionlimit(10000)
from b13_icolearn import (add, neg, sub, mul, sign, Z0, PHI, ICO, VERTS, qmul2, rot, four, VIDX,
                          TABLES, dot)

# ---------------- 床：ペンタグリッドから菱形の床を作る ----------------
GAMMA = [0.1377, 0.2913, -0.3511, 0.0621, -0.1400]   # 和が 0
GAMMA[4] = -sum(GAMMA[:4])
ZETA = [(math.cos(2*math.pi*j/5), math.sin(2*math.pi*j/5)) for j in range(5)]

def build(N):
    edges = set(); verts = set()
    for j, l in itertools.combinations(range(5), 2):
        a, b = ZETA[j]; c, d = ZETA[l]; det = a*d - b*c
        for kj in range(-N, N+1):
            for kl in range(-N, N+1):
                rj, rl = kj - GAMMA[j], kl - GAMMA[l]
                x = (rj*d - b*rl)/det; y = (a*rl - rj*c)/det
                K = [math.ceil(x*ZETA[i][0] + y*ZETA[i][1] + GAMMA[i]) for i in range(5)]
                cs = []
                for e1, e2 in ((0, 0), (1, 0), (1, 1), (0, 1)):
                    k = list(K); k[j] = kj + e1; k[l] = kl + e2; cs.append(tuple(k))
                for i in range(4):
                    u, v = cs[i], cs[(i+1) % 4]
                    edges.add((min(u, v), max(u, v))); verts.add(u)
    return verts, edges

def xpos(k):   # 2x を Z[φ] で：cos72 = (φ−1)/2, cos144 = −φ/2
    return add(add((2*k[0], 0), mul((k[1]+k[4], 0), (-1, 1))), mul((-(k[2]+k[3]), 0), PHI))
def fpos(k):
    return (sum(k[i]*ZETA[i][0] for i in range(5)), sum(k[i]*ZETA[i][1] for i in range(5)))

# ---------------- 回転を 3×3 の表に（2 倍で持つ） ----------------
E3 = [tuple([Z0] + [((1 if i == j else 0), 0) for j in range(3)]) for i in range(3)]
def mat2(q):   # 2R（成分は Z[φ]）
    cols = []
    for e in E3:
        r = rot(q, e)                       # 4R e
        cols.append([(c[0]//2, c[1]//2) for c in r[1:]])
        assert all(c[0] % 2 == 0 and c[1] % 2 == 0 for c in r[1:])
    return [[cols[c][r] for c in range(3)] for r in range(3)]
MAT = {q: mat2(q) for q in ICO}
def apply(M, v):
    return tuple(add(add(mul(M[r][0], v[0]), mul(M[r][1], v[1])), mul(M[r][2], v[2])) for r in range(3))
def vadd(u, v): return tuple(add(a, b) for a, b in zip(u, v))
def vscale(s, v): return tuple(mul((s, 0), a) for a in v)
def vdot(u, v):
    r = Z0
    for a, b in zip(u, v): r = add(r, mul(a, b))
    return r
V3 = [v[1:] for v in VERTS]                 # 12 頂点（3 成分）
ZERO3 = (Z0, Z0, Z0)
IDQ = ((2, 0), Z0, Z0, Z0)

def readout(y):  # 着地：内積が最大の頂点（同点は番号の小さい方）
    best, bv = 0, vdot(y, V3[0])
    for i in range(1, 12):
        d = vdot(y, V3[i])
        if sign(sub(d, bv)) > 0: best, bv = i, d
    return best

# ---------------- 網：入力の列（層 0）から出力の一点へ、段を一つずつ進む ----------------
def make_net(N=5, R=3.6, m=6):
    verts, edges = build(N)
    verts = [v for v in verts if math.hypot(*fpos(v)) < R]
    vs = set(verts)
    adj = {v: [] for v in verts}
    for u, v in edges:
        if u in vs and v in vs: adj[u].append(v); adj[v].append(u)
    # 出力：x がもっとも大きい 1 個（並べるのは Z[φ] の比較）
    import functools
    cmp = lambda a, b: sign(sub(xpos(a), xpos(b)))
    order = sorted(verts, key=functools.cmp_to_key(cmp))
    out = order[-1]
    # 出力からの段数（BFS）。入力は段数 L の番地から m 個（x の小さい順）
    dout = {out: 0}; frontier = [out]
    while frontier:
        nf = []
        for u in frontier:
            for w in adj[u]:
                if w not in dout: dout[w] = dout[u] + 1; nf.append(w)
        frontier = nf
    L = max(d for v, d in dout.items() if sum(1 for w in dout if dout[w] == d) >= m)
    inputs = [v for v in order if dout.get(v) == L][:m]
    dist = {v: L - d for v, d in dout.items() if d <= L}
    # 段 ℓ→ℓ+1 の辺だけを残し、入力から届く番地だけ残す
    fwd = {v: [w for w in adj[v] if w in dist and dist[w] == dist[v] + 1] for v in dist}
    reach = set(inputs)
    for l in range(0, L):
        for v in [v for v in dist if dist[v] == l and v in reach]:
            reach.update(fwd[v])
    nodes = sorted(reach, key=lambda v: (dist[v], v))
    fwd = {v: [w for w in fwd[v] if w in reach] for v in nodes}
    return dict(nodes=nodes, fwd=fwd, dist=dist, inputs=inputs, out=out, L=L)

def count_paths(net):
    c = {net['out']: 1}
    for v in reversed(net['nodes']):
        if v != net['out']: c[v] = sum(c[w] for w in net['fwd'][v])
    return sum(c[v] for v in net['inputs'])

# 入力 i は、固定の頂点 IN_V[i] に符号 s_i（−1, 0, +1）を掛けて入る
IN_V = [0, 1, 2, 3, 4, 5]   # 正反対でない 6 頂点を下で選び直す

def basis(net, orient):
    """各入力から出力までの経路の和（線形）。入力 i ごとに出力のベクトルを返す"""
    res = []
    for i, src in enumerate(net['inputs']):
        val = {src: V3[IN_V[i]]}
        for v in net['nodes']:
            if v not in val: continue
            y = apply(MAT[orient[v]], val[v])          # その番地の向きで回す
            if v == net['out']: res.append(y); break
            for w in net['fwd'][v]:
                val[w] = vadd(val.get(w, ZERO3), y)    # 隣へ足す
        else:
            res.append(ZERO3)
    return res

def outputs(B, X):
    ys = []
    for s in X:
        y = ZERO3
        for si, b in zip(s, B):
            if si: y = vadd(y, b if si > 0 else vscale(-1, b))
        ys.append(y)
    return ys

def evaluate(B, X, T):
    ys = outputs(B, X)
    c = 0; m = Z0
    for y, t in zip(ys, T):
        if readout(y) == t: c += 1
        m = add(m, vdot(y, V3[t]))
    return c, m

def better(a, b):  # (当たり数, 内積の和) の辞書式比較
    if a[0] != b[0]: return a[0] > b[0]
    return sign(sub(a[1], b[1])) > 0

def learn_net(net, X, T, steps, passes=30):
    orient = {v: IDQ for v in net['nodes']}
    cur = evaluate(basis(net, orient), X, T); n = 0
    for _ in range(passes):
        moved = False
        for v in net['nodes']:
            best, bs = None, cur
            for g in steps:
                o2 = dict(orient); o2[v] = qmul2(g, orient[v])
                s2 = evaluate(basis(net, o2), X, T)
                if better(s2, bs): best, bs = o2, s2
            if best is not None: orient, cur = best, bs; moved = True; n += 1
        if not moved: break
    return orient, n

def learn_single(X, T, steps, passes=30):
    """単体：全入力を一つの五芒星に集めて一度だけ回す"""
    q = IDQ
    def ev(q):
        M = MAT[q]
        return evaluate([apply(M, V3[IN_V[i]]) for i in range(len(IN_V))], X, T)
    cur = ev(q)
    for _ in range(passes * 10):
        best, bs = None, cur
        for g in steps:
            q2 = qmul2(g, q); s2 = ev(q2)
            if better(s2, bs): best, bs = q2, s2
        if best is None: break
        q, cur = best, bs
    return q

STEP72 = [q for q in ICO if q[0] == PHI]
STEPALL = [q for q in ICO if q[0] != (2, 0) and q[0] != (-2, 0)]

if __name__ == "__main__":
    # 正反対を含まない 6 頂点
    chosen = []
    for i in range(12):
        if all(VERTS[i][1:] != tuple(neg(c) for c in VERTS[j][1:]) for j in chosen): chosen.append(i)
    IN_V[:] = chosen[:6]
    net = make_net()
    m = len(net['inputs'])
    print(f"担体：番地 {len(net['nodes'])} 個、段数 {net['L']}、入力 {m} 個、"
          f"入力から出力への経路 {count_paths(net)} 本")
    ALLX = [s for s in itertools.product((-1, 0, 1), repeat=m) if any(s)]
    rng = random.Random(13)
    NTR = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    TRIALS = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    STEPS = STEPALL if (len(sys.argv) > 3 and sys.argv[3] == 'all') else STEP72
    agg = dict(before=0, net=0, single=0, maj=0, tot=0, ctrl=0, ctrl_tot=0, ctrl_maj=0, netall=0)
    for t in range(TRIALS):
        teacher = {v: rng.choice(ICO) for v in net['nodes']}
        TB = basis(net, teacher)
        labels = [readout(y) for y in outputs(TB, ALLX)]
        idx = list(range(len(ALLX))); rng.shuffle(idx)
        tr, te = idx[:NTR], idx[NTR:]
        Xtr = [ALLX[i] for i in tr]; Ttr = [labels[i] for i in tr]
        Xte = [ALLX[i] for i in te]; Tte = [labels[i] for i in te]
        maj = max(Ttr, key=Ttr.count); majhit = sum(1 for y in Tte if y == maj)
        before = evaluate(basis(net, {v: IDQ for v in net['nodes']}), Xte, Tte)[0]
        o, n = learn_net(net, Xtr, Ttr, STEPS)
        hit = evaluate(basis(net, o), Xte, Tte)[0]
        q1 = learn_single(Xtr, Ttr, STEPS)
        s1 = evaluate([apply(MAT[q1], V3[IN_V[i]]) for i in range(m)], Xte, Tte)[0]
        # 負の対照：札を入れ替える
        Tsh = list(Ttr); rng.shuffle(Tsh)
        Tte_sh = [rng.choice(Ttr) for _ in Tte]   # 試験の札も同じ分布から無関係に
        oc, _ = learn_net(net, Xtr, Tsh, STEPS)
        ch = evaluate(basis(net, oc), Xte, Tte_sh)[0]
        cmaj = sum(1 for y in Tte_sh if y == max(Tsh, key=Tsh.count))
        agg['before'] += before; agg['net'] += hit; agg['single'] += s1; agg['maj'] += majhit
        agg['tot'] += len(Tte); agg['ctrl'] += ch; agg['ctrl_maj'] += cmaj
        print(f"[{t}] 札の種類 {len(set(labels))}  試験 {len(Tte)}：学習前 {before}  網 {hit}（{n} 歩）  "
              f"単体 {s1}  多数決 {majhit}  ｜ 対照 {ch}（多数決 {cmaj}）", flush=True)
    T_ = agg['tot']
    print(f"合計（試験 {T_}）：学習前 {agg['before']}  網 {agg['net']}  単体 {agg['single']}  "
          f"多数決 {agg['maj']}  ｜ 対照 {agg['ctrl']}（多数決 {agg['ctrl_maj']}、12 分の 1 は {T_//12}）")
