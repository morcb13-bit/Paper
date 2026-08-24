#  検定FL10  流し撮りで切り出す
#      膜が追いつけていない分 W = L - S の正の側を、動いているものとして取る。
#      OK なら：適合率が高い（背景を拾わない）。再現率は動きの速さで決まる
#      NG なら：適合率が背景の割合そのもの＝切り出せていない
#      比較対象（既測定・検定M1〜M3）：同じ一歩で続くセルを集める形は
#      2枚 0.44 / 3枚 0.78 / 4枚 0.92 / 5枚 1.00（適合率）、再現率 0.89→0.47
#
#  検定FL11  止まったら凝視
#      対象が途中で止まる。動く荷が消えるまでの枚数と、そのとき視線が対象の上にいるか。
#      OK なら：止まってから数枚で動く荷が消え、視線は対象の上に残る
#      NG なら：消えない（止まったことが出ない）、または視線が対象から外れる
#      必ず落ちる設定：止まらない対象では動く荷が消えないこと

import math, random
_t4 = open("track4.py").read().split("# ---- 走らせる")[0]
exec(_t4.split("import math, random")[1])

TAU2 = 2 * math.pi

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

def scene2(frames=60, step=0.45, bg=100, flick=8, seed=11, stop_at=None):
    rng = random.Random(seed)
    base = rng.sample(range(NC), bg)
    th = rng.uniform(0, TAU2)
    pos = [-6.0, -4.0]
    out = []
    for t in range(frames):
        exc = set(base)
        for _ in range(flick):
            exc.add(rng.randrange(NC))
        obj = set(nearest(pos[0] + p[0], pos[1] + p[1]) for p in SHAPE)
        exc |= obj
        out.append(([CELLS[i] for i in sorted(exc)],
                    tuple(pos),
                    set(CELLS[i] for i in obj)))
        s = 0.0 if (stop_at is not None and t >= stop_at) else step
        pos[0] += s * math.cos(th)
        pos[1] += s * math.sin(th)
        if abs(pos[0]) > 12 or abs(pos[1]) > 12:
            th += math.pi * 0.5
    return out

def relax(frames, tau):
    """毎枚の W = L - S を返す（膜が追いつけていない分）"""
    S, out = {}, []
    for exc, true, obj in frames:
        cur = set(exc)
        W = {}
        for c in cur:
            W[c] = 1.0 - S.get(c, 0.0)
        for c, s0 in S.items():
            if c not in cur:
                W[c] = -s0
        out.append((W, true, obj, cur))
        for c in cur:
            S[c] = S.get(c, 0.0) + (1.0 - S.get(c, 0.0)) / tau
        for c in list(S.keys()):
            if c not in cur:
                S[c] -= S[c] / tau
                if S[c] < 0.02:
                    del S[c]
    return out

print("検定FL10  動いているものの切り出し（W > 0.5 の側）")
print("%-8s %-6s %10s %10s %10s" % ("一歩", "τ", "適合率", "再現率", "取った枚数"))
for step in (0.45, 0.9, 1.8):
    for tau in (2.0, 6.0, 20.0):
        P, R, K = [], [], []
        for sd in (11, 23, 37):
            for W, true, obj, cur in relax(scene2(step=step, seed=sd), tau):
                sel = set(c for c, w in W.items() if w > 0.5)
                if not sel:
                    continue
                hit = len(sel & obj)
                P.append(hit / len(sel)); R.append(hit / len(obj)); K.append(len(sel))
        print("%-8.2f %-6.0f %10.2f %10.2f %10.1f" %
              (step, tau, sum(P)/len(P), sum(R)/len(R), sum(K)/len(K)))
print("  対照：画面の励起に占める対象の割合 ≒ %.2f（切り出さずに取ったときの適合率）"
      % (24 / 132.0))

print()
print("検定FL11  止まったら凝視（30枚目で停止・τ=2・γ=1.5＝電車のつまみ）")

def run_stop(seed, stop_at, tau=2.0, gamma=1.5, k=6.0, sub=24, rad=9.0):
    fr = scene2(frames=60, seed=seed, stop_at=stop_at)
    rl = relax(fr, tau)
    x, y = fr[0][1]
    vx = vy = 0.0
    dt = 1.0 / sub
    sig, dev, cover = [], [], []
    for (W, true, obj, cur) in rl:
        wl = [(c, w) for c, w in W.items() if abs(w) > 0.02]
        for _ in range(sub):
            gx, gy, n = grad(wl, x, y)
            vx += dt * (-k * gx - gamma * vx)
            vy += dt * (-k * gy - gamma * vy)
            x += dt * vx; y += dt * vy
        s = sum(abs(w) for c, w in wl if math.hypot(c[0]-x, c[1]-y) < rad)
        sig.append(s)
        dev.append(math.hypot(x - true[0], y - true[1]))
        near = [c for c in cur if math.hypot(c[0]-x, c[1]-y) < rad]
        cover.append((len(set(near) & obj) / len(obj), len(set(near) & obj) / max(len(near), 1)))
    return sig, dev, cover

for label, stop_at in (("30枚目で止まる", 30), ("止まらない（対照）", None)):
    S3 = [run_stop(sd, stop_at) for sd in (11, 23, 37)]
    mov = [sum(s[0][20:30]) / 10 for s in S3]
    for d in (2, 4, 6, 10):
        aft = [sum(s[0][30 + d - 2:30 + d]) / 2 for s in S3]
        print("  %-16s 停止%2d枚後の動く荷 %.2f（動いている間 %.2f）"
              % (label if d == 2 else "", d, sum(aft)/3, sum(mov)/3))
    d30 = [sum(s[1][25:30]) / 5 for s in S3]
    print("  %-16s 停止直前（25〜30枚）の視線のずれ %.1f" % ("", sum(d30)/3))
    dv = [sum(s[1][40:60]) / 20 for s in S3]
    cv = [sum(c[0] for c in s[2][40:60]) / 20 for s in S3]
    cp = [sum(c[1] for c in s[2][40:60]) / 20 for s in S3]
    print("  %-16s 40〜60枚の視線のずれ %.1f／凝視の範囲に対象が %.2f 入る（その範囲の純度 %.2f）"
          % ("", sum(dv)/3, sum(cv)/3, sum(cp)/3))
    print()
