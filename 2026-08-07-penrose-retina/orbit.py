#  検定FL5  初速を与えた軌道は閉じるか（抵抗なし）
#
#  沈み V(x) = h(u(x)) - h(u0)   u0 は止まり先、x は u0 のまわりの接平面
#  運動  x'' = -grad V           抵抗を入れない
#  測る量：近点から次の近点までの回り角（頂角）。
#
#      OK（＝軌道が閉じる・楕円）なら：近点から近点まで 180 度ちょうど（沈みが二次。
#                                     井戸の底を中心とする楕円）または 360 度ちょうど
#                                     （沈みが逆二乗。焦点を中心とする楕円）で、
#                                     初速にも物体にも依らない
#      NG なら：頂角がそれ以外の値になり、初速や物体で動く。軌道は閉じない
#
#  対照（判別法2）：二次の井戸と逆二乗の井戸を同じ積分器で通し、
#                   90 度・180 度が出ることを先に確かめる。出なければ積分器の欠陥。

import math, random
exec(open("fall.py").read().split("# ---- 走らせる")[0])

def tangent(u0):
    return basis(u0)

def make_sink_V(pts, u0):
    e1, e2 = tangent(u0)
    h0 = sink(pts, u0)
    def V(x):
        w = tuple(u0[k] + x[0] * e1[k] + x[1] * e2[k] for k in range(3))
        n = math.sqrt(dot(w, w))
        return max(dot(p, w) for p in pts) / n - h0
    def G(x):
        w = tuple(u0[k] + x[0] * e1[k] + x[1] * e2[k] for k in range(3))
        n = math.sqrt(dot(w, w))
        ps = max(pts, key=lambda p: dot(p, w))
        pw = dot(ps, w)
        return (dot(ps, e1) / n - pw * x[0] / n ** 3,
                dot(ps, e2) / n - pw * x[1] / n ** 3)
    return V, G

def V_log(x):
    return math.log(math.hypot(*x))          # ゴム膜：中央の錘による凹み
def G_log(x):
    r2 = x[0] ** 2 + x[1] ** 2
    return (x[0] / r2, x[1] / r2)

def V_quad(x):  return 0.5 * (x[0] ** 2 + x[1] ** 2)
def G_quad(x):  return (x[0], x[1])
def V_kep(x):
    r = math.hypot(*x); return -1.0 / r
def G_kep(x):
    r = math.hypot(*x); return (x[0] / r ** 3, x[1] / r ** 3)

def apsidal(V, G, x0, v0, dt=2e-5, nmax=40_000_000, want=6):
    """近点から近点までの回り角を返す"""
    x = list(x0); v = list(v0)
    a = [-g for g in G(x)]
    E0 = 0.5 * (v[0] ** 2 + v[1] ** 2) + V(x)
    rprev, rprev2 = math.hypot(*x), None
    angs, tot, last = [], 0.0, None
    th_prev = math.atan2(x[1], x[0])
    dE = 0.0
    for it in range(nmax):
        for k in (0, 1):
            v[k] += 0.5 * dt * a[k]
            x[k] += dt * v[k]
        a = [-g for g in G(x)]
        for k in (0, 1):
            v[k] += 0.5 * dt * a[k]
        th = math.atan2(x[1], x[0])
        d = th - th_prev
        if d > math.pi:  d -= 2 * math.pi
        if d < -math.pi: d += 2 * math.pi
        tot += d
        th_prev = th
        r = math.hypot(*x)
        if rprev2 is not None and rprev < rprev2 and rprev < r:      # 近点
            if last is not None:
                angs.append(abs(tot - last))
            last = tot
            if len(angs) >= want:
                break
        rprev2, rprev = rprev, r
        dE = max(dE, abs(0.5 * (v[0] ** 2 + v[1] ** 2) + V(x) - E0))
    return angs, dE, E0

def report(name, V, G, r0, frac):
    g = G((r0, 0.0))
    vc = math.sqrt(r0 * math.hypot(*g))          # 円軌道の速さ（どの力でも v^2 = r|F|）
    v0 = (0.0, frac * vc)
    angs, dE, E0 = apsidal(V, G, (r0, 0.0), v0)
    if not angs:
        print("  %-16s 近点が取れない" % name); return
    d = [math.degrees(a) for a in angs]
    print("  %-16s 頂角 %7.2f 度 (幅 %.3f)   エネルギー保存 %.1e" %
          (name, sum(d) / len(d), max(d) - min(d), dE))

print("対照（積分器の確認：閉じる二つ）")
report("二次の井戸",   V_quad, G_quad, 1.0, 0.6)
report("逆二乗の井戸", V_kep,  G_kep,  1.0, 0.8)
print()
print("ゴム膜（中央の錘が作る凹み）")
for f in (0.9, 0.7, 0.5):
    report("対数の井戸 初速%.1f" % f, V_log, G_log, 1.0, f)

print()
print("沈みの井戸（止まり先のまわり）")
rng = random.Random(13)
for name in ("丸-9点", "平たい丸", "ばらの12点", "丸-2点"):
    pts = OBJ[name]
    u0, _, _ = fall(pts, rand_dir(random.Random(7)))
    V, G = make_sink_V(pts, u0)
    for frac in (0.9, 0.7, 0.5):
        report("%s 初速%.1f" % (name, frac), V, G, 0.02, frac)
