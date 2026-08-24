#  落ちる模型 ── 沈みの浅い方へ落ちる。候補も表も持たない。
#
#  沈み h(u) = max_i (p_i - c)·u        c は点の重心、u は下向き
#  落ちる  : いま接している点 p* の接線成分だけ u を傾ける
#  止まる  : h の局所最小。そこが番地
#
#  検定FL0 落ちるか
#      OK なら：どの向きから始めても h が単調に減り、有限歩で止まる。
#               止まった先を揺らすと必ず h が上がる（局所最小である）
#      NG なら：止まらない／揺らすと下がる。落ちる模型として成立しない
#
#  検定FL1 止まり先はいくつあるか
#      OK なら：物体ごとに止まり先の個数と深さが決まる
#      NG なら：どの物体でも同じ。番地にならない
#
#  検定FL2 向きが揃うか（引き込み域）
#      OK なら：ばらまいた初期向きが少数の番地に落ちる。120通りを読む必要が消える
#      NG なら：初期向きの数だけ止まり先がある。落としても揃わない
#
#  検定FL3 必ず落ちる設定
#      丸い立体では、止まり先が全部同じ深さで、引き込み域が等分になること。
#      ここが落ちなければ検査が反証不能である。

import math, random

PHI = (1 + 5 ** 0.5) / 2

def nrm(p):
    r = math.sqrt(sum(x * x for x in p))
    return tuple(x / r for x in p)

def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

def ico():
    v = set()
    for s1 in (1, -1):
        for s2 in (1, -1):
            v.add(nrm((0.0, s1 * 1.0, s2 * PHI)))
            v.add(nrm((s1 * 1.0, s2 * PHI, 0.0)))
            v.add(nrm((s2 * PHI, 0.0, s1 * 1.0)))
    return sorted(v)

def sphere42():
    v = ico()
    pts = list(v)
    for i in range(len(v)):
        for j in range(i + 1, len(v)):
            if dot(v[i], v[j]) > 0.4:                      # 最近接のみ（1/sqrt5=0.447）
                pts.append(nrm(tuple(v[i][k] + v[j][k] for k in range(3))))
    return pts

D0 = nrm((0.31, 0.77, 0.56))                               # 削る向き（対称軸に乗らない）

def drop_top(pts, k):
    order = sorted(range(len(pts)), key=lambda i: -dot(pts[i], D0))
    cut = set(order[:k])
    return [p for i, p in enumerate(pts) if i not in cut]

def centered(pts):
    n = len(pts)
    c = [sum(p[k] for p in pts) / n for k in range(3)]
    return [tuple(p[k] - c[k] for k in range(3)) for p in pts]

# ---- 物体 ----------------------------------------------------------------
S = sphere42()
OBJ = {
    "丸（対照）":   S,
    "丸-2点":      drop_top(S, 2),
    "丸-5点":      drop_top(S, 5),
    "丸-9点":      drop_top(S, 9),
    "丸+外2点":    S + [tuple(1.30 * x for x in D0),
                        tuple(1.30 * x for x in nrm((-0.62, 0.14, 0.77)))],
    "ばらの12点":  ico(),
    "平たい丸":    [(p[0], p[1], 0.40 * p[2]) for p in S],
    "細長い丸":    [(p[0], p[1], 2.50 * p[2]) for p in S],
}
OBJ = {k: centered(v) for k, v in OBJ.items()}

# ---- 沈みと落下 ----------------------------------------------------------
def sink(pts, u):
    return max(dot(p, u) for p in pts)

def contact(pts, u, tol=1e-6):
    h = sink(pts, u)
    return frozenset(i for i, p in enumerate(pts) if dot(p, u) > h - tol)

def minnorm2(Q):
    """2次元の点集合の凸包のうち原点に最も近い点。原点が内側なら (0,0)。"""
    k = len(Q)
    if k == 1:
        return Q[0]
    # 原点が凸包の内側か（2次元カラテオドリ：三角形・線分・点で尽きる）
    for i in range(k):
        if Q[i][0] ** 2 + Q[i][1] ** 2 < 1e-24:
            return (0.0, 0.0)
    for i in range(k):
        for j in range(i + 1, k):
            a, b = Q[i], Q[j]
            cr = a[0] * b[1] - a[1] * b[0]
            if abs(cr) < 1e-18 and a[0] * b[0] + a[1] * b[1] < 0:
                return (0.0, 0.0)
            for l in range(j + 1, k):
                c = Q[l]
                d1 = a[0] * b[1] - a[1] * b[0]
                d2 = b[0] * c[1] - b[1] * c[0]
                d3 = c[0] * a[1] - c[1] * a[0]
                if (d1 >= -1e-18 and d2 >= -1e-18 and d3 >= -1e-18) or \
                   (d1 <= 1e-18 and d2 <= 1e-18 and d3 <= 1e-18):
                    return (0.0, 0.0)
    best, bn = None, float("inf")
    for i in range(k):
        n = Q[i][0] ** 2 + Q[i][1] ** 2
        if n < bn:
            best, bn = Q[i], n
    for i in range(k):
        for j in range(i + 1, k):
            a, b = Q[i], Q[j]
            dx, dy = b[0] - a[0], b[1] - a[1]
            dd = dx * dx + dy * dy
            if dd < 1e-24:
                continue
            t = -(a[0] * dx + a[1] * dy) / dd
            if 0.0 < t < 1.0:
                p = (a[0] + t * dx, a[1] + t * dy)
                n = p[0] ** 2 + p[1] ** 2
                if n < bn:
                    best, bn = p, n
    return best

def basis(u):
    a = (1.0, 0.0, 0.0) if abs(u[0]) < 0.9 else (0.0, 1.0, 0.0)
    e1 = nrm(tuple(a[k] - dot(a, u) * u[k] for k in range(3)))
    e2 = (u[1] * e1[2] - u[2] * e1[1],
          u[2] * e1[0] - u[0] * e1[2],
          u[0] * e1[1] - u[1] * e1[0])
    return e1, e2

def fall_core(pts, u, rounds=60):
    """沈みの浅い方へ落ちる。接している点が複数あるときは、その組が許す
       いちばん急な向きへ傾ける。沈みが増える一歩は踏まない。"""
    h = sink(pts, u)
    eps, delta = 0.3, 1e-2
    steps, mono = 0, True
    for _ in range(rounds):
        for _ in range(200):
            A = [p for p in pts if dot(p, u) > h - delta]
            e1, e2 = basis(u)
            Q = [(dot(p, e1), dot(p, e2)) for p in A]
            m = minnorm2(Q)
            if m[0] ** 2 + m[1] ** 2 < 1e-20:
                break                                   # この粗さでは止まっている
            d = tuple(m[0] * e1[k] + m[1] * e2[k] for k in range(3))
            moved = False
            for _ in range(40):
                v = nrm(tuple(u[k] - eps * d[k] for k in range(3)))
                hv = sink(pts, v)
                if hv < h - 1e-16:
                    u, h, moved, steps = v, hv, True, steps + 1
                    break
                eps *= 0.5
            if not moved:
                break
        eps, delta = max(eps, 1e-3) , delta * 0.25
        if delta < 1e-13:
            break
    return u, steps, mono

def fall(pts, u):
    u, steps, mono = fall_core(pts, u)
    for _ in range(6):                                    # 稜の上（接点2）は釣り合いだが安定でない
        if len(contact(pts, u, 1e-9)) >= 3:
            break
        e1, e2 = basis(u)
        a = 1e-7
        best, bh = u, sink(pts, u)
        for th in (0, 1, 2, 3):
            c, sn = math.cos(th * math.pi / 2), math.sin(th * math.pi / 2)
            v = nrm(tuple(u[k] + a * (c * e1[k] + sn * e2[k]) for k in range(3)))
            v, st, _ = fall_core(pts, v)
            if sink(pts, v) < bh - 1e-15:
                best, bh = v, sink(pts, v)
        if best is u:
            break
        u = best
    return u, steps, mono

def rand_dir(rng):
    while True:
        p = (rng.gauss(0, 1), rng.gauss(0, 1), rng.gauss(0, 1))
        if sum(x * x for x in p) > 1e-6:
            return nrm(p)

def is_local_min(pts, u, rng, trials=200, ang=2e-3):
    h = sink(pts, u)
    for _ in range(trials):
        d = rand_dir(rng)
        dt = tuple(d[k] - dot(d, u) * u[k] for k in range(3))
        n = math.sqrt(sum(x * x for x in dt))
        if n < 1e-12:
            continue
        v = nrm(tuple(u[k] + ang * dt[k] / n for k in range(3)))
        if sink(pts, v) < h - 1e-12:
            return False
    return True

# ---- 走らせる ------------------------------------------------------------
N = 400
rng = random.Random(13)

print("検定FL0/FL1/FL2  初期向き %d 本" % N)
print()
print("%-11s %4s %6s %6s %8s %8s %10s %8s %6s" %
      ("物体", "点数", "止まり先", "深さの種類", "最深", "最浅",
       "最深へ落ちる", "一番深い一つ", "平均歩数"))

rows = {}
diag = {}
for name, pts in OBJ.items():
    rng2 = random.Random(13)
    rest = {}
    steps_all, mono_all, lmin_all = [], True, True
    for _ in range(N):
        u0 = rand_dir(rng2)
        u, st, mono = fall(pts, u0)
        steps_all.append(st)
        mono_all &= mono
        key = contact(pts, u, 1e-5)
        rest.setdefault(key, [0, sink(pts, u), u])
        rest[key][0] += 1
    bad, sizes = 0, {}
    for key, (cnt, h, u) in rest.items():
        sizes[len(key)] = sizes.get(len(key), 0) + 1
        if not is_local_min(pts, u, rng):
            lmin_all = False
            bad += 1
    diag[name] = (bad, sizes)
    depths = sorted(v[1] for v in rest.values())
    kinds = len(set(round(d, 6) for d in depths))
    dmin = depths[0]
    deep = sum(v[0] for v in rest.values() if v[1] < dmin + 1e-6) / N
    top1 = max(v[0] for v in rest.values()) / N
    rows[name] = (len(rest), kinds, dmin, depths[-1], deep, mono_all, lmin_all, top1)
    print("%-11s %4d %6d %8d %10.4f %8.4f %10.1f%% %9.1f%% %6.0f" %
          (name, len(pts), len(rest), kinds, dmin, depths[-1],
           100 * deep, 100 * top1, sum(steps_all) / N))

print()
print("単調に減ったか / 止まり先は局所最小か")
for name, r in rows.items():
    print("  %-12s 単調 %s   局所最小 %s   非最小 %d 個   接点数 %s" %
          (name, "OK" if r[5] else "NG", "OK" if r[6] else "NG",
           diag[name][0], sorted(diag[name][1].items())))

# ---- 検定FL3 必ず落ちる設定 ---------------------------------------------
print()
print("検定FL3  丸い立体は揃わないこと")
n_rest, kinds, dmin, dmax, deep, _, _, top1 = rows["丸（対照）"]
print("  止まり先 %d 個 / 深さの種類 %d / 深さの幅 %.5f / 一番深い一つ %.1f%%（等分なら %.1f%%）"
      % (n_rest, kinds, dmax - dmin, 100 * top1, 100.0 / n_rest))
print("  判定：", "落ちている（揃わない）" if kinds <= 3 and top1 < 3.0 / n_rest
      else "落ちていない ── 検査を疑うこと")
