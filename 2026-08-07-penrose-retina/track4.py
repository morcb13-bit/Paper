#  検定FL8  動く荷だけが井戸を作る（膜に緩和の時定数τを持たせる）
#  検定FL9  必ず落ちる設定
#
#  膜の休み形 S は荷 L へ 1/τ の速さで追いつく。井戸を作るのは追いつけていない分 W = L - S。
#  動かない背景は τ 枚のうちに均されて消える。動く荷だけが凹みとして残り、
#  通り過ぎた後ろには盛り上がりが残る（前へ押す）。
#
#      OK なら：対象を置かない画面で視線が動かない（FL9 が落ちる）。
#               かつ追跡が FL6（膜をならさない場合）以上を保つ
#      NG なら：FL9 が落ちない、または追跡が悪化する
#
#  対照（判別法2）：τ を大きくすると膜はならされず、FL6 の結果に戻ること。
#
#  膜   V(x) = - Σ_{励起セル c} exp(-|x-c|^2 / 2σ^2)      局所の核。比較も候補も無い
#  視線 x'' = -grad V - γ x'                              落ちるだけ
#
#      OK なら：減衰を小さくすると一定速度の対象に併走でき（流し撮り）、
#               大きくすると止まった対象の底に落ち着く（凝視）。
#               検定CS（見失うまで 9〜17枚・ずれ 5.2〜9.0）を超える
#      NG なら：どちらの減衰でも見失う。井戸では追えない
#
#  検定FL7：対象を置かず背景だけにすると、視線は行き先を持たないこと。
#           ここが落ちなければ検査が反証不能である。
#
#  担体は走らせ役が作った切断投影の点集合であって、監督の五角形の担体ではない。
#  最近接の間隔を 1.618 に揃えてあるので一歩の尺度は揃うが、個数は直接は比べられない。

import math, random

TAU = 2 * math.pi

# ---- 担体（ペンタグリッドの双対でひし形の頂点を出す） -------------------
def penrose_cells(K=7, R=22.0, seed=3):
    rng = random.Random(seed)
    e = [(math.cos(TAU * j / 5), math.sin(TAU * j / 5)) for j in range(5)]
    g = [rng.uniform(0.1, 0.9) for _ in range(4)]
    g.append(-sum(g))                                   # Σγ = 0
    V = set()
    for j in range(5):
        for l in range(j + 1, 5):
            det = e[j][0] * e[l][1] - e[j][1] * e[l][0]
            for kj in range(-K, K + 1):
                for kl in range(-K, K + 1):
                    a, b = kj - g[j], kl - g[l]
                    px = (a * e[l][1] - b * e[j][1]) / det
                    py = (b * e[j][0] - a * e[l][0]) / det
                    if px * px + py * py > R * R:
                        continue
                    idx = []
                    for p in range(5):
                        if p == j:   idx.append(kj)
                        elif p == l: idx.append(kl)
                        else:        idx.append(math.ceil(px * e[p][0] + py * e[p][1] + g[p]))
                    vx = sum(idx[p] * e[p][0] for p in range(5))
                    vy = sum(idx[p] * e[p][1] for p in range(5))
                    V.add((round(vx, 6), round(vy, 6)))
    V = [v for v in V if v[0] ** 2 + v[1] ** 2 <= R * R]
    # 最近接の間隔を 1.618 に揃える
    s = 1.6180339887
    return [(s * x, s * y) for x, y in V]

CELLS = penrose_cells()
NC = len(CELLS)

def nearest(x, y):
    bi, bd = 0, 1e18
    for i, c in enumerate(CELLS):
        d = (c[0] - x) ** 2 + (c[1] - y) ** 2
        if d < bd:
            bi, bd = i, d
    return bi

# ---- 場面 ----------------------------------------------------------------
def triangle24(size=7.0):
    pts = []
    V = [(0.0, size), (-size * 0.866, -size * 0.5), (size * 0.866, -size * 0.5)]
    for i in range(3):
        a, b = V[i], V[(i + 1) % 3]
        for t in range(8):
            f = t / 8.0
            pts.append((a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1])))
    return pts

SHAPE = triangle24()

def scene(frames=60, step=0.45, turn=None, bg=100, flick=8, seed=11, with_obj=True):
    rng = random.Random(seed)
    base = rng.sample(range(NC), bg)
    th = rng.uniform(0, TAU)
    pos = [-6.0, -4.0]
    out = []
    for t in range(frames):
        if turn and t and t % turn == 0:
            th += rng.uniform(-1.2, 1.2)
        exc = set(base)
        for _ in range(flick):
            exc.add(rng.randrange(NC))
        if with_obj:
            for p in SHAPE:
                exc.add(nearest(pos[0] + p[0], pos[1] + p[1]))
        out.append(([CELLS[i] for i in sorted(exc)], tuple(pos)))
        pos[0] += step * math.cos(th)
        pos[1] += step * math.sin(th)
        if abs(pos[0]) > 12 or abs(pos[1]) > 12:            # 画角に留める
            th += math.pi * 0.5
    return out

# ---- 膜と視線 ------------------------------------------------------------
SIG = 2.6
CUT = 3.0 * SIG

def grad(exc, x, y):
    gx = gy = 0.0
    n = 0
    inv = 1.0 / (SIG * SIG)
    for c, wgt in exc:
        dx, dy = x - c[0], y - c[1]
        if abs(dx) > CUT or abs(dy) > CUT:
            continue
        r2 = dx * dx + dy * dy
        if r2 > CUT * CUT:
            continue
        n += 1
        w = wgt * math.exp(-0.5 * r2 * inv) * inv
        gx += w * dx
        gy += w * dy
    return gx, gy, n                                       # V = -Σexp → -gradV = -(w dx)

def track(frames, gamma, sub=24, k=6.0, lost=6.0, tau=6.0):
    x, y = frames[0][1]
    vx = vy = 0.0
    dt = 1.0 / sub
    dev, touched, lost_at = [], 0, None
    x0, y0, path, px, py = x, y, 0.0, x, y
    S = {}
    for t, (exc, true) in enumerate(frames):
        cur = set(exc)
        W = []
        for c in cur:
            w = 1.0 - S.get(c, 0.0)
            if abs(w) > 0.02:
                W.append((c, w))
        for c, s0 in list(S.items()):
            if c not in cur and abs(s0) > 0.02:
                W.append((c, -s0))
        for _ in range(sub):
            gx, gy, n = grad(W, x, y)
            touched += n
            ax = -k * gx - gamma * vx
            ay = -k * gy - gamma * vy
            vx += dt * ax; vy += dt * ay
            x += dt * vx;  y += dt * vy
        for c in cur:
            S[c] = S.get(c, 0.0) + (1.0 - S.get(c, 0.0)) / tau
        for c in list(S.keys()):
            if c not in cur:
                S[c] -= S[c] / tau
                if S[c] < 0.02:
                    del S[c]
        path += math.hypot(x - px, y - py); px, py = x, y
        d = math.hypot(x - true[0], y - true[1])
        dev.append(d)
        if lost_at is None and d > lost:
            lost_at = t
    return ((lost_at if lost_at is not None else len(frames)), sum(dev) / len(dev),
            touched / len(frames) / sub, math.hypot(x - x0, y - y0), path)

# ---- 走らせる ------------------------------------------------------------
print("担体 %d セル（最近接 1.618）／背景100・ちらつき8／60枚・一歩0.45" % NC)
print()
print("検定FL8  膜の緩和 τ（γ=1.5 と γ=12.0）")
print("%-14s %11s %9s %9s %9s %11s" %
      ("τ / 減衰γ", "向き不変", "8枚ごと", "4枚ごと", "2枚ごと", "1枚ごと"))
for tau in (2.0, 6.0, 20.0, 1e6):
    for gamma in (1.5, 12.0):
        row = []
        for turn in (None, 8, 4, 2, 1):
            L, D = [], []
            for sd in (11, 23, 37):
                fr = scene(turn=turn, seed=sd)
                l, d, tt, _, _ = track(fr, gamma, tau=tau)
                L.append(l); D.append(d)
            row.append("%.0f枚/%.1f" % (sum(L)/3, sum(D)/3))
        lab = ("ならさず" if tau > 1e5 else "%.0f" % tau) + " / γ=%.1f" % gamma
        print("%-14s %11s %9s %9s %9s %11s" % (lab, *row))

print()
print("検定FL9  視線自身がどれだけ動いたか（始点からの距離 / 通った道のり）")
print("%-12s %22s %22s" % ("τ / 減衰γ", "対象あり（追う）", "対象なし（背景だけ）"))
for tau in (2.0, 6.0, 20.0, 1e6):
    for gamma in (1.5, 12.0):
        a = track(scene(seed=11), gamma, tau=tau)
        b = track(scene(seed=11, with_obj=False), gamma, tau=tau)
        lab = ("ならさず" if tau > 1e5 else "%.0f" % tau) + " / γ=%.1f" % gamma
        print("%-12s %10.1f / %-10.1f %10.1f / %-10.1f" % (lab, a[3], a[4], b[3], b[4]))
print("  対象は60枚で約20動く。対象なしで同じだけ動くなら、落ちる設定が落ちていない")

print()
print("一枚あたりに触ったセル数： %.1f 個" % track(scene(seed=11), 1.5, tau=6.0)[2])
