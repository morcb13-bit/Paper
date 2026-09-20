#  検定CG1d  種を振って、正しい一歩が山の1位に出る回数を数える
#
#      検定CG1a は種 5 の一枚の結果でしかない。担体の差か、その種のたまたまかを分ける。
#
#  OK なら：どちらかの担体が一貫して1位を取る（回数の差が種のばらつきを超える）
#  NG なら：どちらも種ごとにばらつき、差が出ない
#
#  前提：wind_core.py・mg_lib.py・cg1_hex.py を同じ場所に置く。

import math

src = open("cg1_hex.py").read()
exec(src[:src.index('if __name__ == "__main__":')])
exec(open("mg_lib.py").read())

TRI_P = poly(7.0, 3)
SQR_P = poly(6.0, 4, turn=15.0)

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
for s in (lo, hi, (lo + hi) / 2):
    pts = hex_in(H, s, cx, cy)
    if best is None or abs(len(pts) - 628) < abs(len(best) - 628):
        best = pts
hexr = Retina("六角", best)

pvA = min(STEPS[1:], key=lambda v: abs(math.hypot(*U.xy(v)) - 3.1))
pvB = max(STEPS[1:], key=lambda v: math.hypot(*U.xy(v)) if math.hypot(*U.xy(v)) < 5.2 else 0)
LA, LB = math.hypot(*U.xy(pvA)), math.hypot(*U.xy(pvB))
hvA, hvB = hexr.step_near(LA), hexr.step_near(LB)


def pen_rank(seed, nz=80):
    cA, cB = spot(-7, 6), spot(8, -7)
    A1, B1 = put(TRI, cA), put(SQR, cB)
    A2, B2 = put(TRI, U.zadd(cA, pvA)), put(SQR, U.zadd(cB, pvB))
    E1 = A1 | B1 | noise(nz, seed)
    E2 = A2 | B2 | noise(nz, seed + 1)
    sc = sorted(((len(support(E1, E2, d)), i) for i, d in enumerate(STEPS)), reverse=True)
    rank = [i for _, i in sc].index(STEPS.index(pvA)) + 1
    g = support(E1, E2, pvA)
    p, r = pr(g, A1)
    return rank, p, r


def hex_rank(seed, nz=80):
    d = case(hexr, hvA, hvB, nz=nz, seed=seed)
    return d['rank'], d['p'], d['r']


print("検定CG1d  種を振る（二枚・雑音80）")
print("  種    ペンローズ 順位/適合率/再現率     六角 順位/適合率/再現率")
pw = hw = 0
for seed in range(5, 25, 2):
    a = pen_rank(seed)
    b = hex_rank(seed)
    pw += a[0] == 1
    hw += b[0] == 1
    print("   %2d      %2d 位  %.2f  %.2f            %2d 位  %.2f  %.2f"
          % (seed, a[0], a[1], a[2], b[0], b[1], b[2]))
print("  1位を取った回数　ペンローズ %d/10　六角 %d/10" % (pw, hw))

print("\n  雑音を変える（種5・三角形の一歩を追う）")
print("  雑音   ペンローズ 順位   六角 順位")
for nz in (0, 40, 80, 150, 300):
    a = pen_rank(5, nz)
    b = hex_rank(5, nz)
    print("   %4d        %2d 位        %2d 位" % (nz, a[0], b[0]))


print("\n検定CG1e  雑音を増やして種を振る（二枚・1位を取った回数）")
print("  雑音    ペンローズ   六角")
for nz in (80, 150, 300, 450):
    pw = sum(1 for s in range(5, 25, 2) if pen_rank(s, nz)[0] == 1)
    hw = sum(1 for s in range(5, 25, 2) if hex_rank(s, nz)[0] == 1)
    print("   %4d      %2d/10     %2d/10" % (nz, pw, hw))

print("\n検定CG1f  枚数を増やしたとき、正しい一歩がほかを抜くか（雑音80・種を振る）")
print("  担体      枚数   抜いた回数  支持セルの中央値  ほかの最大の中央値")
for name in ("ペンローズ", "六角"):
    for T in (3, 5):
        win, mine, oth = 0, [], []
        for seed in range(5, 25, 2):
            if name == "ペンローズ":
                Es, tr = frames(T, pvA, pvB, seed=seed)
                g = len(track(Es, pvA))
                o = max(len(track(Es, d)) for d in STEPS if d != pvA)
            else:
                Es, tr = frames_r(hexr, T, hvA, hvB, seed=seed)
                g = len(track_r(hexr, Es, hvA))
                o = max(len(track_r(hexr, Es, d)) for d in hexr.STEPS if d != hvA)
            win += g > o
            mine.append(g); oth.append(o)
        mine.sort(); oth.sort()
        print("   %-8s   %d     %2d/10        %4d            %4d"
              % (name, T, win, mine[5], oth[5]))
