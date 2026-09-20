#  検定CG1  六角格子とペンローズ担体を、同じ規則で突き合わせる
#
#      発端（引継書 v255 §7-3 / v256 §8-3）：同じ受容器数・同じ課題・同じ雑音・
#                    同じ検定画像で六角とペンローズを回す。六角側の追跡規則が未定義。
#
#  規則は一つしか書かない。担体だけ差し替える。
#      着地      いちばん近い受容器
#      一歩の候補 中心の受容器から半径6以内にある受容器への変位（＋ゼロ）
#      支持セル  一歩 d で動かした先が次の画面でも点いている受容器
#      適合率    支持セルのうち実際の図形に属する割合
#      再現率    実際の図形のうち支持セルに入った割合
#
#  揃えるもの
#      受容器数 628（ペンローズ担体と同数）／視野はペンローズ担体の凸包
#      六角の間隔は、その視野に628個ちょうど入るように探して決める
#      図形（三角形 R=7.0・四角形 R=6.0 を15°回したもの）・雑音の枚数・乱数の種
#
#  事前登録（走らせる前に書いた）
#
#  検定CG1a  正しい一歩が山の1位に出るか
#      OK なら：両方の担体で1位
#      NG なら：どちらかが雑音に埋もれる
#
#  検定CG1b  枚数を増やしたときの適合率・再現率（測定。両方を並べる）
#
#  検定CG1c  必ず落ちる設定：二つが同じ一歩で動くと、両方の担体で分けられない
#      OK なら：両方の担体で、二つの図形が同じ支持セルに乗る
#      NG なら：分けられてしまう＝この検査は NG を返せない
#
#  自己検査  平面座標で書き直したペンローズ側が、スキルの motion_group.py と
#            同じ数（三角形 支持39個・適合率0.44・再現率0.89）を返すか
#
#  前提：wind_core.py を exec する。

import math, random
from collections import defaultdict

GAP = 1.618034


# ───────────────────── 担体 ─────────────────────

class Retina:
    """受容器の集まり。着地と一歩だけを持つ。"""

    def __init__(self, name, pts):
        self.name = name
        self.P = list(pts)
        self.g = 3.0
        self.B = defaultdict(list)
        for i, (x, y) in enumerate(self.P):
            self.B[(int(x // self.g), int(y // self.g))].append(i)
        cx = sum(x for x, y in self.P) / len(self.P)
        cy = sum(y for x, y in self.P) / len(self.P)
        self.c = min(range(len(self.P)),
                     key=lambda i: (self.P[i][0] - cx) ** 2 + (self.P[i][1] - cy) ** 2)
        ox, oy = self.P[self.c]
        self.STEPS = [(0.0, 0.0)]
        for x, y in self.P:
            d2 = (x - ox) ** 2 + (y - oy) ** 2
            if 0 < d2 <= 36:
                self.STEPS.append((x - ox, y - oy))

    def land(self, x, y):
        gx, gy = int(x // self.g), int(y // self.g)
        best, bd = None, None
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for i in self.B.get((gx + dx, gy + dy), ()):
                    px, py = self.P[i]
                    d = (px - x) ** 2 + (py - y) ** 2
                    if bd is None or d < bd:
                        bd, best = d, i
        return best

    def put(self, shape, c):
        out = set()
        for px, py in shape:
            i = self.land(c[0] + px, c[1] + py)
            if i is not None:
                out.add(i)
        return out

    def noise(self, n, seed):
        return set(random.Random(seed).sample(range(len(self.P)), n))

    def support(self, E1, E2, d):
        out = set()
        for i in E1:
            x, y = self.P[i]
            j = self.land(x + d[0], y + d[1])
            if j is not None and j in E2:
                out.add(i)
        return out

    def step_near(self, L):
        """長さが L にいちばん近い一歩（ゼロを除く）。"""
        return min(self.STEPS[1:], key=lambda v: abs(math.hypot(*v) - L))


def hull(P):
    """凸包（単調鎖）。"""
    Q = sorted(set(P))
    def half(Q):
        h = []
        for p in Q:
            while len(h) >= 2 and ((h[-1][0] - h[-2][0]) * (p[1] - h[-2][1])
                                   - (h[-1][1] - h[-2][1]) * (p[0] - h[-2][0])) <= 0:
                h.pop()
            h.append(p)
        return h
    return half(Q)[:-1] + half(Q[::-1])[:-1]


def inside(H, p):
    for k in range(len(H)):
        a, b = H[k], H[(k + 1) % len(H)]
        if (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]) < -1e-9:
            return False
    return True


def hex_in(H, s, cx, cy):
    """間隔 s の六角格子のうち、凸包 H の中に入る点。"""
    out = []
    R = 120
    for j in range(-R, R + 1):
        y = cy + j * s * math.sqrt(3) / 2
        for i in range(-R, R + 1):
            x = cx + (i + (j % 2) * 0.5) * s
            if inside(H, (x, y)):
                out.append((x, y))
    return out


def poly(R, n, turn=0.0, M=8):
    V = [(R * math.cos(math.radians(90 + turn + 360.0 * k / n)),
          R * math.sin(math.radians(90 + turn + 360.0 * k / n))) for k in range(n)]
    pts = []
    for k in range(n):
        a, b = V[k], V[(k + 1) % n]
        for i in range(M):
            t = i / M
            pts.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
    return pts


TRI_P = poly(7.0, 3)
SQR_P = poly(6.0, 4, turn=15.0)


def pr(got, truth):
    if not got:
        return 0.0, 0.0
    return len(got & truth) / len(got), len(got & truth) / len(truth)


# ───────────────────── 検定 ─────────────────────

def case(R, vA, vB, nz=80, seed=5):
    """二枚の画面。三角形は vA、四角形は vB で動く。"""
    cA = (R.P[R.c][0] - 7, R.P[R.c][1] + 6)
    cB = (R.P[R.c][0] + 8, R.P[R.c][1] - 7)
    A1, B1 = R.put(TRI_P, cA), R.put(SQR_P, cB)
    A2 = R.put(TRI_P, (cA[0] + vA[0], cA[1] + vA[1]))
    B2 = R.put(SQR_P, (cB[0] + vB[0], cB[1] + vB[1]))
    E1 = A1 | B1 | R.noise(nz, seed)
    E2 = A2 | B2 | R.noise(nz, seed + 1)
    sc = sorted(((len(R.support(E1, E2, d)), k) for k, d in enumerate(R.STEPS)),
                reverse=True)
    rank = [k for _, k in sc].index(R.STEPS.index(vA)) + 1
    g = R.support(E1, E2, vA)
    p, r = pr(g, A1)
    gb = R.support(E1, E2, vB)
    pb, rb = pr(gb, B1)
    return dict(top=sc[0][0], med=sorted(s for s, _ in sc)[len(sc) // 2],
                rank=rank, n=len(g), p=p, r=r,
                nb=len(gb), pb=pb, rb=rb, mix=len(g & gb))


def frames_r(R, T, vA, vB, nz=80, seed=5):
    cA = (R.P[R.c][0] - 7, R.P[R.c][1] + 6)
    cB = (R.P[R.c][0] + 8, R.P[R.c][1] - 7)
    Es, truth = [], []
    for t in range(T):
        a = R.put(TRI_P, (cA[0] + t * vA[0], cA[1] + t * vA[1]))
        b = R.put(SQR_P, (cB[0] + t * vB[0], cB[1] + t * vB[1]))
        truth.append(a)
        Es.append(a | b | R.noise(nz, seed + t))
    return Es, truth


def track_r(R, Es, d):
    keep = set(Es[0])
    for t in range(len(Es) - 1):
        nxt = set()
        for i in keep:
            x, y = R.P[i]
            j = R.land(x + (t + 1) * d[0], y + (t + 1) * d[1])
            if j is not None and j in Es[t + 1]:
                nxt.add(i)
        keep = nxt
    return keep


def run_cg1(RET):
    NG = 0
    print("担体の条件")
    print("  名前        受容器数   一歩の候補   中心から半径6以内")
    for R in RET:
        print(f"    {R.name:8s}    {len(R.P):4d}      {len(R.STEPS):3d} 通り      "
              f"{len(R.STEPS)-1} 個")

    # 図形の一歩は、各担体が自分で持っている候補のうち長さが近いものを取る
    steps = {R.name: (R.step_near(3.078), R.step_near(4.980)) for R in RET}
    print("\n  図形の一歩（画面座標の長さ）")
    for R in RET:
        a, b = steps[R.name]
        print(f"    {R.name:8s}  三角形 {math.hypot(*a):.3f}　四角形 {math.hypot(*b):.3f}")

    # ── 検定CG1a ─────────────────────────────────────────
    print("\n検定CG1a  正しい一歩が山の1位に出るか（二枚・雑音80）")
    print("  担体        山の1位  中央値  順位   支持セル  適合率  再現率")
    ok = True
    for R in RET:
        vA, vB = steps[R.name]
        d = case(R, vA, vB)
        print(f"    {R.name:8s}   {d['top']:4d}    {d['med']:4d}   {d['rank']:2d} 位"
              f"    {d['n']:4d}    {d['p']:.2f}    {d['r']:.2f}")
        ok = ok and d['rank'] == 1
    print(f"  両方とも1位  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    # ── 検定CG1b ─────────────────────────────────────────
    print("\n検定CG1b  枚数を増やす（雑音80・三角形の一歩を追う）")
    print("  担体        枚数  支持セル  適合率  再現率   ほかの一歩の最大")
    for R in RET:
        vA, vB = steps[R.name]
        for T in (2, 3, 4, 5):
            Es, truth = frames(R, T, vA, vB)
            g = track(R, Es, vA)
            p, r = pr(g, truth[0])
            other = max(len(track(R, Es, d)) for d in R.STEPS if d != vA)
            print(f"    {R.name:8s}    {T}     {len(g):4d}    {p:.2f}    {r:.2f}"
                  f"       {other:4d} 個")

    # ── 検定CG1c ─────────────────────────────────────────
    print("\n検定CG1c  必ず落ちる設定：二つが同じ一歩で動く")
    print("  担体        三角形の支持  四角形の支持  重なり  分けられたか")
    ok = True
    for R in RET:
        vA, vB = steps[R.name]
        d = case(R, vB, vB)
        sep = d['mix'] == 0
        print(f"    {R.name:8s}      {d['n']:4d}        {d['nb']:4d}     {d['mix']:4d}"
              f"    {'分けられた' if sep else '分けられない'}")
        ok = ok and not sep
    print(f"  両方とも分けられない  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    print(f"\nNG {NG} / 2")


if __name__ == "__main__":
    # ペンローズ側は元の mg_lib.py（4次元の格子ベクトルで一歩を打つ）をそのまま使う。
    # 平面座標で書き直すと一歩が格子並進でなくなり、規則が壊れる（自己検査で確認した）。
    exec(open("mg_lib.py").read())   # TRI/SQR はここで4次元のものに上書きされる

    penXY = [tuple(float(t) for t in U.xy(q)) for q in CELLS]
    H = hull(penXY)
    cx = sum(x for x, y in penXY) / len(penXY)
    cy = sum(y for x, y in penXY) / len(penXY)
    lo, hi = 1.0, 4.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if len(hex_in(H, mid, cx, cy)) > 628:
            lo = mid
        else:
            hi = mid
    best = None
    for sgap in (lo, hi, (lo + hi) / 2):
        pts = hex_in(H, sgap, cx, cy)
        if best is None or abs(len(pts) - 628) < abs(len(best) - 628):
            best = pts
    hexr = Retina("六角", best)

    pvA = min(STEPS[1:], key=lambda v: abs(math.hypot(*U.xy(v)) - 3.1))
    pvB = max(STEPS[1:], key=lambda v: math.hypot(*U.xy(v)) if math.hypot(*U.xy(v)) < 5.2 else 0)
    LA, LB = math.hypot(*U.xy(pvA)), math.hypot(*U.xy(pvB))
    hvA, hvB = hexr.step_near(LA), hexr.step_near(LB)

    print("\n担体の条件")
    print("  名前        受容器数   一歩の候補   三角形の一歩   四角形の一歩")
    print("    ペンローズ      %4d      %3d 通り     %.3f        %.3f"
          % (len(CELLS), len(STEPS), LA, LB))
    print("    六角         %4d      %3d 通り     %.3f        %.3f"
          % (len(hexr.P), len(hexr.STEPS), math.hypot(*hvA), math.hypot(*hvB)))

    print("\n自己検査  元の motion_group.py と同じ数が出るか（三角形 39個・0.44・0.89・1位）")
    got, truth = run_case(pvA, pvB, label="■ ペンローズ（元の規則）")

    print("\n検定CG1a  正しい一歩が山の1位に出るか（二枚・雑音80）")
    dh = case(hexr, hvA, hvB)
    print("  六角　山の1位 %d／中央値 %d／順位 %d 位" % (dh['top'], dh['med'], dh['rank']))
    print("        三角形 支持 %d 個・適合率 %.2f・再現率 %.2f" % (dh['n'], dh['p'], dh['r']))
    print("        四角形 支持 %d 個・適合率 %.2f・再現率 %.2f" % (dh['nb'], dh['pb'], dh['rb']))

    print("\n検定CG1b  枚数を増やす（雑音80・三角形の一歩を追う）")
    print("  担体        枚数  支持セル  適合率  再現率   ほかの一歩の最大")
    for T in (2, 3, 4, 5):
        Es, tr = frames(T, pvA, pvB)
        g = track(Es, pvA)
        p, r = pr(g, tr[0])
        other = max(len(track(Es, d)) for d in STEPS if d != pvA)
        print("    ペンローズ       %d     %4d    %.2f    %.2f       %4d 個"
              % (T, len(g), p, r, other))
    for T in (2, 3, 4, 5):
        Es, tr = frames_r(hexr, T, hvA, hvB)
        g = track_r(hexr, Es, hvA)
        p, r = pr(g, tr[0])
        other = max(len(track_r(hexr, Es, d)) for d in hexr.STEPS if d != hvA)
        print("    六角          %d     %4d    %.2f    %.2f       %4d 個"
              % (T, len(g), p, r, other))

    print("\n検定CG1c  必ず落ちる設定：二つが同じ一歩で動く")
    run_case(pvA, pvA, label="■ ペンローズ")
    dc = case(hexr, hvB, hvB)
    print("  六角　三角形の支持 %d 個／四角形の支持 %d 個／重なり %d 個 → %s"
          % (dc['n'], dc['nb'], dc['mix'], '分けられない' if dc['mix'] else '分けられた'))
