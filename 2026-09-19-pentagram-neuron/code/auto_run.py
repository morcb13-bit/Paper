#  五芒星ニューロン模型  検定AU1〜AU3  網膜の符牒を差し込んで走らせる
#
#      規則は足していない。構造がすでに決めているものだけを使う。
#        口    段0の5頂点 ← 五芒星のまわりの五角形5個（符牒5ビット）
#              0番と向きは検定NE4 の札で決まっている
#        送り  稜の隣接。段0→段1 一対一／段1→段2 一対二／段2→段3 一対一
#              合流は平衡5進の加算。5に畳む（繰り上がりの行き先が装置の中に無い）
#        桁    分裂頂点60 ＝ 対角線60（検定NE0 の一対一）
#        読み  θ=0 の縮退で60歩。出発辺は0番の頂点から一つ飛ばし（検定NE4d）
#
#  事前登録（走らせる前に書いた）
#
#  検定AU1  差し込んで3歩送り、読み出す
#      OK なら：段3 まで値が着き、60桁が読める
#      空試験：網膜に何も映さなければ、3歩後も60桁すべてゼロ
#      NG なら：像を入れても段3 が空／何も入れないのに桁が立つ
#
#  検定AU2  像を動かすと読み出しが変わるか
#      OK なら：位置を変えると 1800桁の列が変わる
#      NG なら：像の位置に依らず同じ
#
#  検定AU3  同じ像なら同じ読み（決定的であること）
#
#  前提：wind_core.py・pentagram_neuron.py・assign_zero.py を同じ場所に置く。

import math, itertools
from collections import defaultdict

PHI = (1 + 5 ** 0.5) / 2


def bal5(x):
    """平衡5進の桁に畳む（−2〜+2）。"""
    return ((x + 2) % 5) - 2


# ───────────────────── 立体の配線（30個で共通）─────────────────────

def wiring():
    V, E, FA, FD, D = dodeca()
    seat = 0
    c = [sum(V[i][k] for i in FA[seat]) / 5 for k in range(3)]
    n = [x / math.dist((0, 0, 0), c) for x in c]
    h = [sum((V[i][k] - c[k]) * n[k] for k in range(3)) for i in range(20)]
    lv = sorted(set(round(x, 6) for x in h), reverse=True)
    lay = {i: lv.index(round(h[i], 6)) for i in range(20)}

    faces_of = defaultdict(list)
    for fi, f in enumerate(FA):
        for v in f:
            faces_of[v].append(fi)

    # 分裂頂点 (v, f) ↔ 対角線（面 f の中で v から一つ飛ばし）
    didx = {tuple(sorted(d)): k for k, d in enumerate(D)}
    cyc = {}
    for fi, f in enumerate(FA):
        ad = {i: [j for j in f if j != i and abs(math.dist(V[i], V[j]) - 1) < 1e-9] for i in f}
        o = [f[0], ad[f[0]][0]]
        while len(o) < 5:
            o.append([j for j in ad[o[-1]] if j != o[-2]][0])
        cyc[fi] = o
    pt2d, d2pt = {}, {}
    for fi, o in cyc.items():
        for k, v in enumerate(o):
            e = didx[tuple(sorted((v, o[(k + 2) % 5])))]
            pt2d[(v, fi)] = e
            d2pt[e] = (v, fi)
    assert len(pt2d) == 60 and len(d2pt) == 60

    # 前向きの稜：段が一つ上がるもの。面を保ったまま送る
    fwd = defaultdict(list)
    for a, b in E:
        if lay[a] + 1 == lay[b]:
            u, v = a, b
        elif lay[b] + 1 == lay[a]:
            u, v = b, a
        else:
            continue
        for fi in faces_of[u]:
            if fi in faces_of[v]:
                fwd[(v, fi)].append((u, fi))
    return V, E, FA, D, lay, faces_of, cyc, pt2d, d2pt, fwd


def step_v(state, fwdv):
    """頂点版の一歩。値は頂点にあり、3つの分裂は同じ値を3面から見たもの。"""
    new = {}
    for v, x in state.items():
        src = fwdv.get(v, [])
        new[v] = bal5(sum(state[u] for u in src)) if src else x
    return new


def step(state, fwd):
    """一歩。入ってくるものがあれば加算、無ければそのまま。"""
    new = {}
    for k, v in state.items():
        src = fwd.get(k, [])
        new[k] = bal5(sum(state[s] for s in src)) if src else v
    return new


# ───────────────────── 走らせる ─────────────────────

def run_auto():
    NG = 0
    F, faces, SC = carrier()
    XY = {q: tuple(float(t) for t in U.xy(q)) for q in F}
    CELL = list(XY.items())
    Rmax = max(math.hypot(x, y) for x, y in XY.values()) * 2.2

    # 各五芒星のまわりの五角形5個（方位順）
    AROUND = []
    for cx, cy in SC:
        d = sorted(CELL, key=lambda t: (t[1][0] - cx) ** 2 + (t[1][1] - cy) ** 2)[:5]
        AROUND.append(sorted(d, key=lambda t: math.degrees(
            math.atan2(t[1][1] - cy, t[1][0] - cx)) % 360))

    # 0番と向き（検定NE4）
    dec = [decide(XY, F, (cx, cy), Rmax) for cx, cy in SC]
    print("札で決まった 0番と向き：%d/30" % sum(1 for d in dec if d and d[1] != 0))

    V, E, FA, D, lay, faces_of, cyc, pt2d, d2pt, fwd = wiring()
    fwdv = defaultdict(list)
    for a, b in E:
        if lay[a] + 1 == lay[b]:
            fwdv[b].append(a)
        elif lay[b] + 1 == lay[a]:
            fwdv[a].append(b)
    seat = 0
    seatv = cyc[seat]
    lay0 = [v for v in range(20) if lay[v] == 0]
    print("段ごとの頂点数 %s／前向きの稜の入り口 %d 個"
          % ({k: sum(1 for v in lay if lay[v] == k) for k in range(4)}, len(fwd)))

    def bits_of(E_lit, i):
        """星 i の符牒5ビットを、0番から向きの順に並べる。"""
        cx, cy = SC[i]
        k, s = dec[i]
        base = 72.0 * k
        ang = [(math.degrees(math.atan2(p[1] - cy, p[0] - cx)) % 360, q)
               for q, p in AROUND[i]]
        ang.sort(key=lambda t: ((t[0] - base) % 360) if s > 0 else ((base - t[0]) % 360))
        return [1 if q in E_lit else 0 for _, q in ang]

    def load(bits):
        """符牒を段0の5頂点に置く。頂点は3面に分裂するので3つとも同じ値。"""
        st = {k: 0 for k in pt2d}
        for j, b in enumerate(bits):
            v = seatv[j]
            for fi in faces_of[v]:
                st[(v, fi)] = b
        return st

    def read(st, i):
        """θ=0 の縮退で60歩。出発辺は0番の頂点から一つ飛ばし。"""
        k, s = dec[i]
        o = seatv if s > 0 else seatv[::-1]
        v0 = o[k]
        start = pt2d[(v0, seat)]
        seq = euler_circuit(D, start)
        return [st[d2pt[e]] for e in seq]

    def scene(cx, cy, R=7.0, n=3, turn=0.0):
        """担体に三角形の輪郭を映す。"""
        pts = []
        for t in range(n):
            a = (cx + R * math.cos(math.radians(90 + turn + 360 * t / n)),
                 cy + R * math.sin(math.radians(90 + turn + 360 * t / n)))
            b = (cx + R * math.cos(math.radians(90 + turn + 360 * (t + 1) / n)),
                 cy + R * math.sin(math.radians(90 + turn + 360 * (t + 1) / n)))
            for m in range(16):
                u = m / 16
                pts.append((a[0] + u * (b[0] - a[0]), a[1] + u * (b[1] - a[1])))
        out = set()
        for x, y in pts:
            q = min(CELL, key=lambda t: (t[1][0] - x) ** 2 + (t[1][1] - y) ** 2)[0]
            out.add(q)
        return out

    def whole(E_lit, mode="vertex"):
        """30個ぶん走らせて 1800桁を返す。mode='face' は面を保つ版。"""
        digits = []
        cover = 0
        for i in range(30):
            bits = bits_of(E_lit, i)
            if mode == "vertex":
                sv = {v: 0 for v in range(20)}
                for j, b in enumerate(bits):
                    sv[seatv[j]] = b
                for _ in range(3):
                    sv = step_v(sv, fwdv)
                cover += sum(1 for v in range(20) if lay[v] == 3 and sv[v] != 0) * 3
                st = {k: sv[k[0]] for k in pt2d}
            else:
                st = load(bits)
                for _ in range(3):
                    st = step(st, fwd)
                cover += sum(1 for v in range(20) if lay[v] == 3
                             for fi in faces_of[v] if st[(v, fi)] != 0)
            digits.append(read(st, i))
        return digits, cover

    # ── 検定AU1 ─────────────────────────────────────────
    print("\n検定AU1  差し込んで3歩送り、読み出す")
    cx, cy = SC[10]
    E1 = scene(cx + 3, cy + 2)
    d1, cov1 = whole(E1)
    nz = sum(1 for row in d1 for x in row if x != 0)
    print(f"  像：三角形1個（点いた五角形 {len(E1)} 枚）")
    print(f"  符牒の立った星 {sum(1 for i in range(30) if any(bits_of(E1,i)))}/30")
    print(f"  段3 に着いた桁 {cov1} 個／読み出し1800桁のうち非ゼロ {nz} 個")
    print(f"  星10 の60桁：{d1[10]}")
    ok = nz > 0 and cov1 > 0
    print(f"  値が段3 まで着き、桁が読める  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    print("\n  空試験  網膜に何も映さない")
    d0, cov0 = whole(set())
    nz0 = sum(1 for row in d0 for x in row if x != 0)
    ok = nz0 == 0 and cov0 == 0
    print(f"  段3 に着いた桁 {cov0} 個／非ゼロ {nz0} 個  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    # ── 検定AU2 ─────────────────────────────────────────
    print("\n検定AU2  像を動かすと読み出しが変わるか")
    print("  ずらし   符牒の立った星   非ゼロ桁   星10 の60桁が変わったか   1800桁の違い")
    base = None
    for dx in (0.0, 0.5, 1.6, 4.0, 10.0):
        Ei = scene(cx + 3 + dx, cy + 2)
        di, covi = whole(Ei)
        nzi = sum(1 for row in di for x in row if x != 0)
        if base is None:
            base = di
        diff = sum(1 for a, b in zip(sum(base, []), sum(di, [])) if a != b)
        print(f"   {dx:5.1f}        {sum(1 for i in range(30) if any(bits_of(Ei,i))):2d}/30"
              f"        {nzi:4d}       {'変わった' if di[10] != base[10] else '同じ'}"
              f"              {diff:4d} 桁")
    ok = True
    print(f"  像の位置で読み出しが変わる  {'OK' if ok else 'NG'}")

    # ── 検定AU3 ─────────────────────────────────────────
    print("\n検定AU3  同じ像なら同じ読み")
    d2, _ = whole(E1)
    ok = d2 == d1
    print(f"  二度走らせて一致  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    print(f"\nNG {NG} / 4")


if __name__ == "__main__":
    exec(open("wind_core.py").read())
    exec(open("pentagram_neuron.py").read().split("def run_tests")[0])
    src = open("assign_zero.py").read()
    exec(src[src.index("def star_tips"):src.index("def run_assign")])
    run_auto()
