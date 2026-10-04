# b13_icolearn.py — 五芒星に貼った正二十面体（イコシアン）の向きを、例から加算だけで学ぶ
# 数は Z[φ] の (a, b) = a + bφ、φ² = φ + 1。浮動小数・割り算は使わない。
import itertools, random

# ---- Z[φ] ----
def add(x, y): return (x[0]+y[0], x[1]+y[1])
def neg(x): return (-x[0], -x[1])
def sub(x, y): return add(x, neg(y))
def mul(x, y):
    a, b = x; c, d = y
    # (a+bφ)(c+dφ) = ac + (ad+bc)φ + bd(φ+1)
    return (a*c + b*d, a*d + b*c + b*d)
def sign(x):
    a, b = x
    if a >= 0 and b >= 0: return 0 if (a, b) == (0, 0) else 1
    if a <= 0 and b <= 0: return -1
    if a > 0:            # a + bφ, b<0 : a > cφ ⇔ a² − ac − c² > 0
        c = -b; t = a*a - a*c - c*c
        return 1 if t > 0 else -1
    d = -a               # b>0, a<0 : bφ > d ⇔ d² − db − b² < 0
    t = d*d - d*b - b*b
    return 1 if t < 0 else -1
Z0, Z1, PHI = (0, 0), (1, 0), (0, 1)
PHI_1 = (-1, 1)          # φ − 1 = 1/φ

# ---- 四元数（2 倍で持つ） ----
def qmul(p, q):
    a1, b1, c1, d1 = p; a2, b2, c2, d2 = q
    m = mul
    return (sub(sub(sub(m(a1,a2), m(b1,b2)), m(c1,c2)), m(d1,d2)),
            add(sub(add(m(a1,b2), m(b1,a2)), m(d1,c2)), m(c1,d2)),
            add(add(sub(m(a1,c2), m(b1,d2)), m(c1,a2)), m(d1,b2)),
            add(sub(add(m(a1,d2), m(b1,c2)), m(c1,b2)), m(d1,a2)))
def conj(q): return (q[0], neg(q[1]), neg(q[2]), neg(q[3]))
def half(q):             # 2 倍で持つ四元数どうしの積は 4 倍になる → 2 で割って 2 倍に戻す
    out = []
    for a, b in q:
        if a % 2 or b % 2: return None   # 割り切れなければ外に出た
        out.append((a//2, b//2))
    return tuple(out)
def qmul2(p, q): return half(qmul(p, q))

def icosians():
    out = set()
    two = (2, 0)
    for i in range(4):
        for s in (1, -1):
            v = [Z0]*4; v[i] = (2*s, 0); out.add(tuple(v))
    for ss in itertools.product((1, -1), repeat=4):
        out.add(tuple((s, 0) for s in ss))
    base = [Z0, Z1, PHI, PHI_1]
    even = [p for p in itertools.permutations(range(4))
            if sum(1 for i in range(4) for j in range(i+1, 4) if p[i] > p[j]) % 2 == 0]
    for p in even:
        for ss in itertools.product((1, -1), repeat=3):
            vals = [Z0, Z1, PHI, PHI_1]
            vals = [vals[0]] + [vals[k] if s > 0 else neg(vals[k]) for k, s in zip((1, 2, 3), ss)]
            v = [None]*4
            for pos, src in enumerate(p): v[pos] = vals[src]
            out.add(tuple(v))
    return sorted(out)

ICO = icosians()
assert len(ICO) == 120
S = set(ICO)
assert all(qmul2(p, q) in S for p in ICO for q in ICO)   # 掛け算で閉じる

# ---- 正二十面体の 12 頂点（(0, φ, 1) の向き、純四元数） ----
VERTS = []
for s1 in (1, -1):
    for s2 in (1, -1):
        a = PHI if s1 > 0 else neg(PHI); b = (s2, 0)
        for t in [(Z0, a, b), (b, Z0, a), (a, b, Z0)]:
            VERTS.append((Z0,) + t)

def rot(q, v):           # (2q) v (2q̄) = 4 · 回したもの
    return qmul(qmul(q, v), conj(q))
def four(v): return tuple(mul((4, 0), c) for c in v)
VIDX = {four(v): i for i, v in enumerate(VERTS)}
assert all(rot(q, v) in VIDX for q in ICO for v in VERTS)  # 頂点は頂点に着地

def table(q): return tuple(VIDX[rot(q, v)] for v in VERTS)
TABLES = {}
for q in ICO: TABLES.setdefault(table(q), q)
assert len(TABLES) == 60                                  # 回転は 60 通り

def dot(x, y):           # 純四元数の内積（Z[φ]）
    r = Z0
    for i in (1, 2, 3): r = add(r, mul(x[i], y[i]))
    return r

# ---- 学習の一歩：72° の回転（実部 = φ）だけ ----
STEPS = [q for q in ICO if q[0] == PHI]

def score(q, examples):  # Σ 内積（着地 · 正解）、大きいほどよい
    s = Z0
    for vi, ti in examples:
        s = add(s, dot(rot(q, VERTS[vi]), VERTS[ti]))
    return s

def learn(examples, limit=100):
    q = ((2, 0), Z0, Z0, Z0)      # 恒等から始める
    cur = score(q, examples); n = 0
    while n < limit:
        best, bs = None, cur
        for g in STEPS:
            q2 = qmul2(g, q); s2 = score(q2, examples)
            if sign(sub(s2, bs)) > 0: best, bs = q2, s2
        if best is None: break
        q, cur = best, bs; n += 1
    return q, n

def hits(q, pairs):
    t = table(q)
    return sum(1 for vi, ti in pairs if t[vi] == ti)

if __name__ == "__main__":
    print("イコシアン", len(ICO), "個／回転", len(TABLES), "通り／一歩の候補", len(STEPS), "個")
    rng = random.Random(13)
    targets = list(TABLES.keys())
    for k in (1, 2, 3):
        subsets = list(itertools.combinations(range(12), k))
        tot = hit = full = before = runs = 0; maxstep = 0
        for tt in targets:
            for sub_ in subsets:
                train = [(i, tt[i]) for i in sub_]
                test = [(i, tt[i]) for i in range(12) if i not in sub_]
                q0 = ((2, 0), Z0, Z0, Z0)
                q, n = learn(train)
                h = hits(q, test); h0 = hits(q0, test)
                tot += len(test); hit += h; before += h0; runs += 1
                full += (h == len(test)); maxstep = max(maxstep, n)
        print(f"k={k}: 試行 {runs}  試験の当たり 学習前 {before}/{tot}  学習後 {hit}/{tot}  "
              f"全問正解の試行 {full}/{runs}  最大の歩数 {maxstep}")
    # 負の対照：正解の札を回転と無関係にでたらめに付ける
    for k in (2, 3):
        tot = hit = runs = 0
        for _ in range(2000):
            lab = [rng.randrange(12) for _ in range(12)]
            sub_ = rng.sample(range(12), k)
            train = [(i, lab[i]) for i in sub_]
            test = [(i, lab[i]) for i in range(12) if i not in sub_]
            q, n = learn(train)
            tot += len(test); hit += hits(q, test); runs += 1
        print(f"負の対照 k={k}: 試行 {runs}  試験の当たり {hit}/{tot}（偶然の水準は 1/12 ≒ {tot//12}）")
