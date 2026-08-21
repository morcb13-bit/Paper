#  検定C1-0  空試験 ── 網膜が単独で作る成分は巻き数のどこに出るか
#      対象：Z[ζ10] の輪郭を N 点で標本し、各点を最も近い五角形中心へ着地させたときの
#            変位だけを巻き取り、巻き数 k ごとの取り分を見る。図形そのものは巻き取らない
#      測る量：r = Σ_{2≤|k|≤4} |c_k|²  ÷  Σ_{2≤|k|≤N/2} |c_k|²
#              （k=0 は重心、|k|=1 は大きさと向き。どちらも形の歪みではないので外す）
#              一様な散り方の期待値は 6/57 = 0.1053
#      OK なら：着地の変位は散る（r ≤ 0.21）── 歪みは低い巻き数、網膜の署名は散り、
#               二つは別の場所に出る。検定C1 の本体に進んでよい
#      NG なら：着地の変位が低い巻き数に集中する（r ≥ 0.32）── 網膜が単独で作る成分が
#               歪みと同じ場所に出る。巻き取りでは分けられない
#      判定不能：0.21 < r < 0.32
#
#      【事前登録の欠陥・走らせたあとに判明】この閾値は N=60 のときの一様値 0.1053 を
#      基準にしている。N を変えると一様値も動く（N=15 では 0.5455）ので、掃引には
#      そのまま使えない。下の「掃引」は一様値で割った値と乱数の帯で見ているが、
#      これは事前登録された判定ではない。次に測るときは帯のほうを先に登録すること
#      対照（同時に走らせる）：
#          雑音 …… 同じ二乗平均の乱数変位。散るはず（散らなければ物差しが壊れている）
#          歪み …… 同じ二乗平均のアフィン変形 z → z + ε·conj(z)。
#                   集中するはず（集中しなければ物差しが壊れている）
#      必ず落ちる設定：N=3 では巻き取りが三成分で尽き、復元が完全一致する（v227 §2.1）
#
#  厳密性の内訳
#      着地の判定（どのセルが最も近いか）…… Q(φ) の符号だけ。浮動小数を使わない
#      候補の絞り込み …… 画面座標で半径3以内を拾うだけ（間隔1.618に対し誤差1e-15）
#      巻き取りの評価 …… 50桁の十進。変位そのものは Q(ζ10) で厳密に持っている

from fractions import Fraction as F
from decimal import Decimal as D, getcontext
from collections import defaultdict
import random, sys

import b13_chain_units as U

getcontext().prec = 50
PI = D("3.14159265358979323846264338327950288419716939937510582097494459230781")

# ── 50桁の三角関数（Taylor・浮動小数を通さない） ─────────────
def dsin(x):
    t = x; s = x; n = 1
    while abs(t) > D(10) ** -46:
        n += 2
        t = -t * x * x / (n * (n - 1))
        s += t
    return s

def dcos(x):
    t = D(1); s = D(1); n = 0
    while abs(t) > D(10) ** -46:
        n += 2
        t = -t * x * x / (n * (n - 1))
        s += t
    return s

# ζ10^k の 50桁表示
Z10 = [(dcos(PI * D(k) / 5), dsin(PI * D(k) / 5)) for k in range(10)]
ROT18 = (dcos(-PI / 10), dsin(-PI / 10))          # b13_chain_units の描画回転 −18°

def to_xy(a):
    """Q(ζ10) の元 (a0,a1,a2,a3) → 50桁の (x,y)。b13_chain_units.xy と同じ向き"""
    x = D(0); y = D(0)
    for k in range(4):
        c = D(a[k].numerator) / D(a[k].denominator) if isinstance(a[k], F) else D(a[k])
        x += c * Z10[k][0]; y += c * Z10[k][1]
    return (x * ROT18[0] - y * ROT18[1], x * ROT18[1] + y * ROT18[0])

# ── 担体 ─────────────────────────────────────────────
rows, place, offs = U.build_stack()
CELLS = list(U.fits(sum(place, [])))
GRID = defaultdict(list)
for q in CELLS:
    p = U.xy(q); GRID[(int(p[0] // 3), int(p[1] // 3))].append(q)

def land(p):
    """p（Q(ζ10)）に最も近い五角形中心。決定は Q(φ) の符号だけで行う"""
    fx, fy = U.xy(tuple(float(t) for t in p))
    gx, gy = int(fx // 3), int(fy // 3)
    cand = []
    for dx in (-2, -1, 0, 1, 2):
        for dy in (-2, -1, 0, 1, 2):
            cand += GRID.get((gx + dx, gy + dy), ())
    if not cand:
        raise RuntimeError("網膜の外")
    best = None; bn = None
    for q in cand:
        n = U.norm2(U.zsub(tuple(F(t) for t in q), p))
        if bn is None or U.phi_lt(n, bn):
            best, bn = q, n
    return best

# ── 輪郭・変位場・巻き取り ──────────────────────────────
import math
CX, CY = 0.0, 38.0
CTR = min(CELLS, key=lambda q: (U.xy(q)[0] - CX) ** 2 + (U.xy(q)[1] - CY) ** 2)

def find_offset(tx, ty, B=6):
    best = None; bd = None
    for a0 in range(-B, B + 1):
        for a1 in range(-B, B + 1):
            for a2 in range(-B, B + 1):
                for a3 in range(-B, B + 1):
                    x, y = U.xy((a0, a1, a2, a3))
                    d = (x - tx) ** 2 + (y - ty) ** 2
                    if bd is None or d < bd:
                        bd, best = d, (a0, a1, a2, a3)
    return best, bd ** 0.5

def to_xyD(a):
    return to_xy(a)

def ms(ds, N):
    s = D(0)
    for d in ds:
        x, y = to_xy(d); s += x * x + y * y
    return s / N

def wind(ds, N, W):
    P = [to_xy(d) for d in ds]
    out = []
    for k in range(N):
        re = D(0); im = D(0)
        for j in range(N):
            c, s = W[(-j * k) % N]
            re += P[j][0] * c - P[j][1] * s
            im += P[j][0] * s + P[j][1] * c
        out.append(re * re + im * im)
    return out

def ratio(pw, N):
    lo = sum(pw[k] + pw[N - k] for k in (2, 3, 4))
    al = sum(pw[k] + pw[N - k] for k in range(2, N // 2)) + pw[N // 2]
    return lo / al

def run(R, M=20, seed=13, verbose=False):
    V = []
    for k in range(3):
        th = math.radians(90 + 120 * k)
        o, _ = find_offset(R * math.cos(th), R * math.sin(th))
        V.append(U.zadd(CTR, o))
    N = 3 * M
    BASE = []
    for k in range(3):
        a, b = V[k], V[(k + 1) % 3]
        for i in range(M):
            t = F(i, M)
            BASE.append(tuple(F(a[j]) + t * (b[j] - a[j]) for j in range(4)))
    CEN = tuple(sum(p[j] for p in BASE) / N for j in range(4))
    W = [(dcos(2 * PI * D(m) / D(N)), dsin(2 * PI * D(m) / D(N))) for m in range(N)]

    DA = [U.zsub(tuple(F(t) for t in land(p)), p) for p in BASE]
    msA = ms(DA, N)
    random.seed(seed)
    DB0 = [tuple(F(random.randint(-1000, 1000), 1000) for _ in range(4)) for _ in range(N)]
    kB = F(int((msA / ms(DB0, N)).sqrt() * 10 ** 12), 10 ** 12)
    DB = [tuple(kB * c for c in d) for d in DB0]
    DC0 = [U.zconj(U.zsub(p, CEN)) for p in BASE]
    kC = F(int((msA / ms(DC0, N)).sqrt() * 10 ** 12), 10 ** 12)
    DC = [tuple(kC * c for c in d) for d in DC0]

    pw = [wind(d, N, W) for d in (DA, DB, DC)]
    rs = [float(ratio(p, N)) for p in pw]
    side = math.hypot(*(lambda a, b: (U.xy(b)[0] - U.xy(a)[0], U.xy(b)[1] - U.xy(a)[1]))(V[0], V[1]))
    return rs, float(msA.sqrt()), side, pw[0], N

def verdict(r):
    return "集中" if r >= 0.32 else ("散る" if r <= 0.21 else "判定不能")



# ══ 1. 事前登録された判定（辺長19・N=60・M=20固定） ══════════
rs, rms, side, pwA, N = run(11.0, M=20)
print("【事前登録された判定】辺長 %.2f・標本 %d 点・変位rms %.4f" % (side, N, rms))
print("  一様な散り方の期待値 r = 6/57 = 0.1053／集中 r≥0.32／散る r≤0.21")
for nm, r in zip(("着地（網膜が単独で作る成分）", "雑音（対照）", "歪み（対照・アフィン変形）"), rs):
    print("    %-28s r = %.4f  %s" % (nm, r, verdict(r)))
print("  判定：%s" % ("OK ── 網膜が単独で作る成分は散る" if rs[0] <= 0.21
                    else ("NG ── 低い巻き数に集中する" if rs[0] >= 0.32 else "判定不能")))
tot = sum(pwA[k] + pwA[N - k] for k in range(2, N // 2)) + pwA[N // 2]
prof = sorted(((float((pwA[k] + pwA[N - k]) / tot), k) for k in range(2, N // 2)), reverse=True)
print("  着地の成分が乗る巻き数（上位5本／一様なら 1/57=0.0175）：",
      "  ".join("|k|=%d %.3f" % (k, v) for v, k in prof[:5]))

import statistics
def run_multi(R, nseed=40):
    _, _, side0, _, _ = run(R, M=6)
    M = max(3, int(round(side0 / 1.618034)))
    V = []
    for k in range(3):
        th = math.radians(90 + 120 * k)
        o, _ = find_offset(R * math.cos(th), R * math.sin(th))
        V.append(U.zadd(CTR, o))
    N = 3 * M
    BASE = []
    for k in range(3):
        a, b = V[k], V[(k + 1) % 3]
        for i in range(M):
            t = F(i, M)
            BASE.append(tuple(F(a[j]) + t * (b[j] - a[j]) for j in range(4)))
    CEN = tuple(sum(p[j] for p in BASE) / N for j in range(4))
    W = [(dcos(2 * PI * D(m) / D(N)), dsin(2 * PI * D(m) / D(N))) for m in range(N)]
    DA = [U.zsub(tuple(F(t) for t in land(p)), p) for p in BASE]
    msA = ms(DA, N)
    rA = float(ratio(wind(DA, N, W), N))
    DC0 = [U.zconj(U.zsub(p, CEN)) for p in BASE]
    kC = F(int((msA / ms(DC0, N)).sqrt() * 10 ** 12), 10 ** 12)
    rC = float(ratio(wind([tuple(kC * c for c in d) for d in DC0], N, W), N))
    rB = []
    for sd in range(nseed):
        random.seed(1000 + sd)
        DB0 = [tuple(F(random.randint(-1000, 1000), 1000) for _ in range(4)) for _ in range(N)]
        kB = F(int((msA / ms(DB0, N)).sqrt() * 10 ** 12), 10 ** 12)
        rB.append(float(ratio(wind([tuple(kB * c for c in d) for d in DB0], N, W), N)))
    side = math.hypot(U.xy(V[1])[0]-U.xy(V[0])[0], U.xy(V[1])[1]-U.xy(V[0])[1])
    flat = 6.0 / (2 * (N // 2 - 2) + 1)
    rB.sort()
    return side, N, flat, rA, rC, rB

print("\n【掃引・事前登録なし】標本の間隔は五角形の間隔（1.618）に合わせる。取り分は一様値で割って示す")
print("雑音は乱数の種40本の分布（5%〜95%）\n")
print("辺長    N   着地      雑音の帯（40本）      歪み")
for R in (5.0, 8.0, 11.0, 16.0, 22.0):
    side, N, flat, rA, rC, rB = run_multi(R)
    lo, md, hi = rB[2], statistics.median(rB), rB[-3]
    inside = "帯の中" if lo <= rA <= hi else ("帯の上" if rA > hi else "帯の下")
    print("%6.2f %4d  %5.2f    %5.2f 〜 %5.2f (中央%5.2f)  %5.2f   着地は%s"
          % (side, N, rA/flat, lo/flat, hi/flat, md/flat, rC/flat, inside))
