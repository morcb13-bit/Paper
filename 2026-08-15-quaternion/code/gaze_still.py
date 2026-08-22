#  静止した形を、視線を動かして読む
#
#  対象を一行で（判別法9）：世界に固定された輪郭を、担体を格子ベクトル g だけ
#  平行移動させて（＝視線を移して）読んだとき、読み値が g で変わるか。
#
#  機構（B13 の言葉）：一つの視線＝一つの窓 W。視線を d 動かした担体は
#      σ(v) ∈ W + σ(d) の場所。したがって
#      同じ場所を両方で見られる ⟺ σ(v) ∈ W ∩ (W+σ(d))   … shift_window.py の SW1（既確定）
#      重ねて標本する場所      ⟺ σ(v) ∈ W ∪ (W+σ(d))   … こちらが本件
#      残渣が大きい一歩（＝実長の短い一歩）ほど窓は重ならず、新しい場所が乗る。
#
#  検定GS0 一歩ごとの新規率は残渣で決まるか（対照・shift_window の裏取り）
#      OK なら：新規率は |σ(d)| の単調関数で、実長では決まらない。
#               σ(d) が W+(−W) の外接半径 0.832544 を超える一歩では新規率 1.00
#      NG なら：実長で決まる／どちらとも無関係
#      必ず落ちる設定：d=0 の新規率が 0 であること
#
#  検定GS1 一つの視線で読めない形が、別の視線で読めるか
#      対象：v229 §2.8 の水の読み（網の縁から、点いていないセルだけを通って
#            近接2.7 で伝播させ、届かないセルを内部とする）。読みは一切変えない。
#      OK なら：まず既知の表（線0・二本0・三角の輪1・丸R=5.0/4.0/3.0 は1・
#               C字0・Y字0・二重の輪2・丸R=6.5 は漏れる）を視線0で再現し、
#               そのうえで丸R=6.5 が視線によっては 1 になる（＝盲点の(a)型・ずらせば入る）
#      NG なら：どの視線でも漏れる（＝担体の密度が決めており、視線では補えない）
#               あるいは既知の表が再現しない（＝組み直した水が別物）
#      必ず落ちる設定：対象のいない画面で内部 0。同じ視線を繰り返しても読みは変わらない
#
#  検定GS2 重ねると輪郭の取りこぼしが減るか（閾値を持たない量で測る）
#      測る量：輪郭上の各点から最も近い「乗った標本点」までの距離。その最大と平均。
#      OK なら：K（視線の数）とともに単調に減る。残渣の大きい一歩を使ったときだけ減り、
#               残渣の小さい一歩（φ³・単数）ではほとんど減らない
#      NG なら：減らない／残渣の大小で差が出ない
#      必ず落ちる設定：d=0 を K 回重ねても一切変わらないこと
#
#  厳密性の内訳：担体の生成と番地は Q(φ)（b13_chain_units）。着地の判定は
#      mg_lib と同じ浮動小数の最近傍（候補は空間ハッシュ）。距離そのものを測る検定なので
#      符号判定ではなく距離を使う。残渣は 50桁ではなく倍精度（表示に足りる桁だけ）。

import math
from collections import deque
import b13_chain_units as U

rows, place, offs = U.build_stack()
CELLS = list(U.fits(sum(place, [])))
XY = {q: U.xy(q) for q in CELLS}
CX, CY = 0.0, 38.0
CTR = min(CELLS, key=lambda q: (XY[q][0] - CX) ** 2 + (XY[q][1] - CY) ** 2)
OX, OY = XY[CTR]
print("担体 %d 環 / セル %d 枚" % (len(place), len(CELLS)))

# 空間ハッシュ（着地用）
H = {}
for q in CELLS:
    x, y = XY[q]
    H.setdefault((int(x // 3), int(y // 3)), []).append(q)

def land(px, py):
    gx, gy = int(px // 3), int(py // 3)
    best, bd = None, None
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for q in H.get((gx + dx, gy + dy), ()):
                x, y = XY[q]
                d = (x - px) ** 2 + (y - py) ** 2
                if bd is None or d < bd:
                    bd, best = d, q
    return best

# ── 残渣（ガロア共役 ζ10 → ζ10^3） ──────────────────────
def resid(v):
    re = im = 0.0
    for k in range(4):
        th = math.pi * (3 * k) / 5.0
        re += v[k] * math.cos(th); im += v[k] * math.sin(th)
    return math.hypot(re, im)

def realen(v):
    x, y = U.xy(v)
    return math.hypot(x, y)

# ── 近接2.7 の隣（担体は平行移動しても内部の幾何は同じなので一度だけ） ──
NB = {q: [] for q in CELLS}
L = list(CELLS)
for i in range(len(L)):
    xi, yi = XY[L[i]]
    for j in range(i + 1, len(L)):
        xj, yj = XY[L[j]]
        if (xi - xj) ** 2 + (yi - yj) ** 2 <= 2.7 ** 2:
            NB[L[i]].append(L[j]); NB[L[j]].append(L[i])
print("近接2.7 の隣：平均 %.1f 枚" % (sum(len(v) for v in NB.values()) / len(CELLS)))

# ── 図形（世界座標・原点は CTR の位置） ───────────────────
def arc(R, a0, a1, n):
    return [(OX + R * math.cos(a0 + (a1 - a0) * i / n),
             OY + R * math.sin(a0 + (a1 - a0) * i / n)) for i in range(n + 1)]

def seg(p, q, n):
    return [(p[0] + (q[0] - p[0]) * i / n, p[1] + (q[1] - p[1]) * i / n) for i in range(n + 1)]

def poly(R, k, n=40, turn=0.0):
    V = [(OX + R * math.cos(math.radians(90 + turn + 360 * i / k)),
          OY + R * math.sin(math.radians(90 + turn + 360 * i / k))) for i in range(k)]
    out = []
    for i in range(k):
        out += seg(V[i], V[(i + 1) % k], n)
    return out

SHAPES = {
    "一本の線":   seg((OX - 9, OY), (OX + 9, OY), 90),
    "離れた二本": seg((OX - 9, OY + 3), (OX + 9, OY + 3), 90) + seg((OX - 9, OY - 3), (OX + 9, OY - 3), 90),
    "三角の輪":   poly(6.0, 3),
    "丸 R=3.0":   arc(3.0, 0, 2 * math.pi, 120),
    "丸 R=4.0":   arc(4.0, 0, 2 * math.pi, 160),
    "丸 R=5.0":   arc(5.0, 0, 2 * math.pi, 200),
    "丸 R=6.5":   arc(6.5, 0, 2 * math.pi, 260),
    "C字":        arc(5.0, 0.45, 2 * math.pi - 0.45, 200),
    "Y字":        seg((OX, OY), (OX, OY + 7), 70) + seg((OX, OY), (OX - 6, OY - 4), 70) + seg((OX, OY), (OX + 6, OY - 4), 70),
    "二重の輪":   arc(3.0, 0, 2 * math.pi, 120) + arc(7.0, 0, 2 * math.pi, 280),
    "対象なし":   [],
}

# ── 視線 g（格子ベクトル）で読む ─────────────────────────
def lit(shape, g):
    """視線 g の担体で点いたセル（担体の番地で返す）"""
    gx, gy = U.xy(g)
    return set(land(px - gx, py - gy) for px, py in shape)

def water(E):
    """網の縁から点いていないセルだけを通す。届かないセル＝内部。かたまりの数を返す"""
    if not E:
        return 0, 0
    rmax = max(math.hypot(XY[q][0] - OX, XY[q][1] - OY) for q in E)
    seed = [q for q in CELLS if q not in E
            and math.hypot(XY[q][0] - OX, XY[q][1] - OY) > rmax]
    seen = set(seed); dq = deque(seed)
    while dq:
        q = dq.popleft()
        for r in NB[q]:
            if r not in seen and r not in E:
                seen.add(r); dq.append(r)
    inner = [q for q in CELLS if q not in E and q not in seen]
    # 内部のかたまりの数
    iset = set(inner); c = 0; done = set()
    for q in inner:
        if q in done:
            continue
        c += 1; dq = deque([q]); done.add(q)
        while dq:
            p = dq.popleft()
            for r in NB[p]:
                if r in iset and r not in done:
                    done.add(r); dq.append(r)
    return c, len(inner)

def miss(shape, gazes):
    """輪郭の各点から、乗った標本点（世界座標）までの最短距離。最大と平均"""
    pts = []
    for g in gazes:
        gx, gy = U.xy(g)
        for q in lit(shape, g):
            pts.append((XY[q][0] + gx, XY[q][1] + gy))
    pts = list(set(pts))
    mx = 0.0; s = 0.0
    for px, py in shape:
        d = min((px - x) ** 2 + (py - y) ** 2 for x, y in pts) ** 0.5
        mx = max(mx, d); s += d
    return mx, s / len(shape), len(pts)

# ── 一歩の候補 ───────────────────────────────────────
ZERO = (0, 0, 0, 0)
UNITS = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1),
         (-1, 1, -1, 1), (-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0),
         (0, 0, 0, -1), (1, -1, 1, -1)]                     # ζ10^k・実長1.0
PHI3 = None
best = None
for a0 in range(-3, 4):
    for a1 in range(-3, 4):
        for a2 in range(-3, 4):
            for a3 in range(-3, 4):
                v = (a0, a1, a2, a3)
                if v == ZERO:
                    continue
                r = resid(v)
                if r < 0.24 and (best is None or realen(v) < best[0]):
                    best = (realen(v), v)
PHI3 = best[1]

print("\n検定GS0 一歩ごとの新規率（担体を平行移動したとき、元の担体に無い場所の割合）")
print("   %-16s %8s %8s %8s" % ("一歩", "実長", "残渣", "新規率"))
S = set(CELLS)
for name, v in [("0（対照）", ZERO), ("ζ10（十方向の1歩）", UNITS[0]), ("ζ10^2", UNITS[2]),
                ("φ³ の一歩 %s" % (PHI3,), PHI3), ("連続接続", (1, 1, 1, 1))]:
    new = sum(1 for q in CELLS if U.zadd(q, v) not in S)
    print("   %-16s %8.4f %8.4f %8.3f" % (name, realen(v), resid(v), new / len(CELLS)))

print("\n検定GS1 視線ごとの読み（水の内部のかたまり／内部のセル数）")
GAZ = [ZERO] + UNITS
print("   %-12s %s" % ("図形", "視線0  " + "  ".join("g%d" % i for i in range(1, 11))))
for nm, sh in SHAPES.items():
    out = []
    for g in GAZ:
        c, n = water(lit(sh, g)) if sh else (0, 0)
        out.append(c)
    print("   %-12s %s" % (nm, "  ".join("%3d" % c for c in out)))

print("\n   同じ視線を繰り返す（対照）：丸 R=6.5 を 0 で5回 → %s"
      % [water(lit(SHAPES["丸 R=6.5"], ZERO))[0] for _ in range(5)])

print("\n検定GS2 重ねると輪郭の取りこぼしが減るか（取りこぼし＝輪郭の点から最も近い標本点までの距離）")
for nm in ("丸 R=6.5", "三角の輪"):
    sh = SHAPES[nm]
    print("   %s" % nm)
    print("      %-22s %6s %6s %8s" % ("視線の重ね方", "最大", "平均", "標本点"))
    for label, gs in [("0 を1回", [ZERO]),
                      ("0 を5回（対照）", [ZERO] * 5),
                      ("ζ10 の一歩 2回", GAZ[:2]),
                      ("ζ10 の一歩 3回", GAZ[:3]),
                      ("ζ10 の一歩 5回", GAZ[:5]),
                      ("ζ10 の一歩 11回", GAZ),
                      ("φ³ の一歩 2回", [ZERO, PHI3]),
                      ("φ³ の一歩 3回", [ZERO, PHI3, U.zadd(PHI3, PHI3)])]:
        mx, av, n = miss(sh, gs)
        print("      %-22s %6.3f %6.3f %8d" % (label, mx, av, n))
