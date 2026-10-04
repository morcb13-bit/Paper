# b13_fannet.py — 10 枚の扇に着地する番地を、ペンローズの床の繋がりで網にして学ばせる
#
# 番地の値：指し ζ^k（k = 0..9、36° 刻みの扇）か、黙っている（−1）。
# 一歩：届いた指しを数え（10 方向の個数）、和がいちばん近い扇に着地し、自分の向き r を足して隣へ渡す。
#       いちばん近い扇が二つ以上並んだら（打ち消し合い・ちょうど境目）黙る。
# 判定は Z[φ] の整数だけ：2cos(36°·m) = (2,0) (0,1) (−1,1) (1,−1) (0,−1) (−2,0) …
# 学習：外れが減るときだけ、番地の向きを扇一枚（±1）回す。
#
# 網は三通り＋対照：
#   A 一方向   入力から出力へ段を一つずつ進む辺だけ
#   B 床の隣   床の辺を両向きに（戻る・横に渡る辺を含む）
#   C 隣＋φ²   B に、φ² 離れた番地どうしの繋がりを足す
#   L 線形     A と同じ辺で、途中では着地せず個数のまま回して渡し、出力でだけ着地する
import itertools, random, sys
import numpy as np
from b13_iconet import build, xpos, fpos
from b13_icolearn import sign as zsign, sub as zsub

CA = np.array([2, 0, -1, 1, 0, -2, 0, 1, -1, 0], dtype=np.int64)   # 2cos(36m) = CA + CB φ
CB = np.array([0, 1, 1, -1, -1, 0, -1, -1, 1, 1], dtype=np.int64)

def vsign(a, b):
    """a + bφ の符号（配列）。整数だけ"""
    n = a*a + a*b - b*b
    s = np.where((a >= 0) & (b >= 0), np.sign(a + b), 0)
    s = np.where((a <= 0) & (b <= 0), -np.sign(-a - b), s)
    s = np.where((a > 0) & (b < 0), np.sign(n), s)
    s = np.where((a < 0) & (b > 0), -np.sign(n), s)
    return s

def land(C):
    """C: (..., 10) の個数 → 着地した扇（並んだら −1）"""
    A = np.stack([(C * np.roll(CA, k)).sum(-1) for k in range(10)], -1)
    B = np.stack([(C * np.roll(CB, k)).sum(-1) for k in range(10)], -1)
    best = np.zeros(C.shape[:-1], dtype=np.int64)
    ba, bb = A[..., 0].copy(), B[..., 0].copy(); ties = np.zeros_like(best)
    for k in range(1, 10):
        s = vsign(A[..., k] - ba, B[..., k] - bb)
        up = s > 0; eq = s == 0
        best = np.where(up, k, best); ties = np.where(up, 0, np.where(eq, ties + 1, ties))
        ba = np.where(up, A[..., k], ba); bb = np.where(up, B[..., k], bb)
    return np.where(ties > 0, -1, best)

# ---------------- 床と網 ----------------
def zeta_pos(k):   # 番地（5 整数）→ Z[ζ5] の 4 整数（1, ζ, ζ², ζ³）
    return tuple(k[i] - k[4] for i in range(4))
def zpow(j):
    j %= 5
    return (-1, -1, -1, -1) if j == 4 else tuple(1 if i == j else 0 for i in range(4))
def vadd(u, v): return tuple(a + b for a, b in zip(u, v))
def vneg(u): return tuple(-a for a in u)
PHI2_DIRS = set()
for j in range(5):
    phz = vneg(vadd(zpow(j+2), zpow(j+3)))       # φζ^j = −(ζ^{j+2} + ζ^{j+3})
    d = vadd(phz, zpow(j))                       # φ²ζ^j = φζ^j + ζ^j
    PHI2_DIRS.add(d); PHI2_DIRS.add(vneg(d))

def make_floor(R=3.6, m=6):
    import math, functools
    verts, edges = build(5)
    verts = [v for v in verts if math.hypot(*fpos(v)) < R]
    vs = set(verts)
    E = [(u, v) for u, v in edges if u in vs and v in vs]
    adj = {v: set() for v in verts}
    for u, v in E: adj[u].add(v); adj[v].add(u)
    # 床の繋がりが一つにまとまった部分だけ使う
    order = sorted(verts, key=functools.cmp_to_key(lambda a, b: zsign(zsub(xpos(a), xpos(b)))))
    out = order[-1]
    seen = {out}; fr = [out]
    while fr:
        nf = []
        for u in fr:
            for w in adj[u]:
                if w not in seen: seen.add(w); nf.append(w)
        fr = nf
    verts = [v for v in order if v in seen]
    dout = {out: 0}; fr = [out]
    while fr:
        nf = []
        for u in fr:
            for w in adj[u]:
                if w not in dout: dout[w] = dout[u] + 1; nf.append(w)
        fr = nf
    L = max(d for d in dout.values() if sum(1 for w in dout if dout[w] == d) >= m)
    inputs = [v for v in verts if dout[v] == L][:m]
    idx = {v: i for i, v in enumerate(verts)}
    pos = {v: zeta_pos(v) for v in verts}
    byp = {pos[v]: v for v in verts}
    phi2 = set()
    for v in verts:
        for d in PHI2_DIRS:
            w = byp.get(vadd(pos[v], d))
            if w is not None: phi2.add((min(idx[v], idx[w]), max(idx[v], idx[w])))
    nb = set((min(idx[u], idx[v]), max(idx[u], idx[v])) for u, v in E)
    # 一方向：入力からの段数 ℓ → ℓ+1 で、出力へ届くものだけ
    dist = {v: L - dout[v] for v in verts if dout[v] <= L}
    fwd = [(idx[u], idx[w]) for u in dist for w in adj[u] if w in dist and dist[w] == dist[u] + 1]
    return dict(N=len(verts), inputs=[idx[v] for v in inputs], out=idx[out], L=L,
                nb=sorted(nb), phi2=sorted(phi2), fwd=fwd)

def directed(pairs): return [(a, b) for a, b in pairs] + [(b, a) for a, b in pairs]

# ---------------- 走らせる ----------------
def run(net, edges, r, X, T, linear=False):
    """X: (E, m) の −1/0/+1。出力の番地の扇（−1 は黙る）を返す"""
    N = net['N']; Ex = X.shape[0]
    src = np.array([a for a, b in edges], dtype=np.int64)
    dst = np.array([b for a, b in edges], dtype=np.int64)
    inp = np.array(net['inputs'])
    # 入力の指し：+1 → 扇 0、−1 → 扇 5、0 → 黙る
    xin = np.where(X > 0, 0, np.where(X < 0, 5, -1))
    def onehot_in():
        H = np.zeros((Ex, len(inp), 10), dtype=np.int64)
        for i in range(len(inp)):
            on = xin[:, i] >= 0
            H[on, i, xin[on, i]] = 1
        return H
    if linear:   # 途中では着地しない：個数のまま、送り手の向きで回して渡す
        S = np.zeros((Ex, N, 10), dtype=np.int64); S[:, inp, :] = onehot_in()
        for t in range(T):
            C = np.zeros((Ex, N, 10), dtype=np.int64)
            for a, b in edges:
                C[:, b, :] += np.roll(S[:, a, :], int(r[a]), axis=-1)
            C[:, inp, :] = onehot_in()
            S = C
        return land(S[:, net['out'], :])
    P = -np.ones((Ex, N), dtype=np.int64)          # 各番地が送り出す指し（向き込み）
    P[:, inp] = np.where(xin >= 0, (xin + r[inp][None, :]) % 10, -1)   # 入力は外から留める
    for t in range(T):
        C = np.zeros((Ex, N, 10), dtype=np.int64)
        on = P[:, src] >= 0
        e_idx, k_idx = np.nonzero(on)
        np.add.at(C, (e_idx, dst[k_idx], P[e_idx, src[k_idx]]), 1)
        K = land(C)
        P = np.where(K >= 0, (K + r[None, :]) % 10, -1)
        P[:, inp] = np.where(xin >= 0, (xin + r[inp][None, :]) % 10, -1)
    return K[:, net['out']]

POS = np.array([1, 1, 1, 0, 0, 0, 0, 0, 1, 1])     # 扇 0,1,2,8,9 は「はい」（ζ^0 との内積が正）
READ = {"mode": "符号"}     # 符号：扇の向きで読む／発火：黙ったか着地したかで読む
STEPS = {"d": (1, -1)}
def score(out, y):
    if READ["mode"] == "発火":
        pred = (out >= 0).astype(np.int64)
        hit = int((pred == y).sum())
        # 内積の代わりに、当たり数だけで比べる（余白は 0）
        return hit, (0, 0)
    pred = np.where(out >= 0, POS[np.maximum(out, 0)], 0)
    hit = int((pred == y).sum())
    sgn = np.where(y == 1, 1, -1)
    a = np.where(out >= 0, CA[np.maximum(out, 0)], 0) * sgn
    b = np.where(out >= 0, CB[np.maximum(out, 0)], 0) * sgn
    return hit, (int(a.sum()), int(b.sum()))

def better(s, t):
    if s[0] != t[0]: return s[0] > t[0]
    return zsign(zsub(s[1], t[1])) > 0

def learn(net, edges, X, y, T, rng, linear=False, passes=12):
    N = net['N']
    r = np.array([rng.randrange(10) for _ in range(N)], dtype=np.int64)
    cur = score(run(net, edges, r, X, T, linear), y)
    for _ in range(passes):
        moved = False
        order = list(range(N)); rng.shuffle(order)
        for v in order:
            for d in STEPS['d']:
                r2 = r.copy(); r2[v] = (r2[v] + d) % 10
                s2 = score(run(net, edges, r2, X, T, linear), y)
                if better(s2, cur): r, cur, moved = r2, s2, True; break
        if not moved: break
    return r

TASKS = {
    "同符号": lambda X: (X[:, 0] * X[:, 1] == 1).astype(np.int64),     # 線形では解けない
    "一つ目が正": lambda X: (X[:, 0] == 1).astype(np.int64),           # 線形で解ける（確認用）
}

if __name__ == "__main__":
    net = make_floor()
    m = len(net['inputs'])
    print(f"床：番地 {net['N']} 個、入力 {m} 個、入力から出力まで {net['L']} 段、"
          f"隣の辺 {len(net['nb'])}、φ² の繋がり {len(net['phi2'])}、一方向の辺 {len(net['fwd'])}")
    ALL = np.array([s for s in itertools.product((-1, 0, 1), repeat=m) if s[0] != 0 and s[1] != 0])
    NETS = {
        "A 一方向": (net['fwd'], net['L'], False),
        "B 床の隣": (directed(net['nb']), net['L'] + 4, False),
        "C 隣＋φ²": (directed(net['nb'] + net['phi2']), net['L'] + 4, False),
        "L 線形": (net['fwd'], net['L'], True),
    }
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    for tname, f in TASKS.items():
        y_all = f(ALL)
        print(f"\n題「{tname}」 入力 {len(ALL)} 通り（はい {int(y_all.sum())}）、学習用 {NTR}、試験 {len(ALL) - NTR}、{seeds} 回の合計")
        for nname, (edges, T, lin) in NETS.items():
            tr_hit = te_hit = ctl = tot = 0
            for s in range(seeds):
                rng = random.Random(1000 + s)
                idx = list(range(len(ALL))); rng.shuffle(idx)
                tr, te = np.array(idx[:NTR]), np.array(idx[NTR:])
                r = learn(net, edges, ALL[tr], y_all[tr], T, rng, lin)
                tr_hit += score(run(net, edges, r, ALL[tr], T, lin), y_all[tr])[0]
                te_hit += score(run(net, edges, r, ALL[te], T, lin), y_all[te])[0]
                ysh = y_all.copy(); perm = list(range(len(ALL))); rng.shuffle(perm); ysh = ysh[perm]
                rc = learn(net, edges, ALL[tr], ysh[tr], T, rng, lin)
                ctl += score(run(net, edges, rc, ALL[te], T, lin), ysh[te])[0]
                tot += len(te)
            print(f"  {nname}：学習用 {tr_hit}/{NTR*seeds}  試験 {te_hit}/{tot}  ｜ 対照（札を混ぜる）試験 {ctl}/{tot}", flush=True)
