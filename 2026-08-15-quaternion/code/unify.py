#  検定UN1  静止した形と動く形を、一つの輪で扱えるか
#
#      対象を一行で（判別法9）：毎枚 ①沈みの深いところへ視線を寄せ ②水で読む、
#          という一つの輪を回したとき、対象が止まっていても動いていても働くか。
#
#      輪（新しい定数を持ち込まない）
#          ① 沈み（sink.py・半径2.7）の深いところの重心に最も近いセルへ視線を移す
#          ② 寄る先が前の枚と同じセルなら（＝視線が止まった）、装置が持つ10方向へ
#             ζ₁₀ の一歩だけ振る。順に回すだけで、向きは選ばない
#          ③ 水（water.py・網の縁から近接2.7）で内部を読む。直近3枚の多数決が
#             2枚続けて同じなら「固まった」として止める（3枚は GS1b から）
#
#      OK なら：止まっている丸は「閉じている」で固まり、C字は「開いている」で固まり、
#               動く丸は視線が離れずに追随し、対象なしでは固まらない（または開いたまま）
#      NG なら：静止で固まらない／動く対象を見失う／C字が閉じで固まる／対象なしが固まる
#      必ず落ちる設定：対象なしの画面で「閉じている」が出ないこと

exec(open('water.py').read().split('print("■ 検定TP3')[0])
exec(open('sink.py').read().split('print("■ 沈んだところへ移る')[0])

import math

X0, Y0 = O
UNITS = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (-1, 1, -1, 1),
         (-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1), (1, -1, 1, -1)]

HS = {}
for q in CELLS:
    x, y = XYC[q]
    HS.setdefault((int(x // 3), int(y // 3)), []).append(q)

def near(px, py):
    gx, gy = int(px // 3), int(py // 3)
    best = None; bd = None
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for q in HS.get((gx + dx, gy + dy), ()):
                d = (XYC[q][0] - px) ** 2 + (XYC[q][1] - py) ** 2
                if bd is None or d < bd:
                    bd = d; best = q
    return best

def paint2(path):
    return set(near(x, y) for x, y in path)

def circle(R, n, cx=0.0, cy=0.0):
    return [(X0 + cx + R * math.cos(2 * math.pi * i / n),
             Y0 + cy + R * math.sin(2 * math.pi * i / n)) for i in range(n)]

def carc(R, n, a0, a1, cx=0.0, cy=0.0):
    return [(X0 + cx + R * math.cos(a0 + (a1 - a0) * i / n),
             Y0 + cy + R * math.sin(a0 + (a1 - a0) * i / n)) for i in range(n + 1)]

NETS = set(NET)

def run(shape_at, frames=30, back=0, seed=13):
    """shape_at(t) → その枚の世界座標の輪郭。back>0 ならその枚数だけ背景を混ぜる"""
    import random as _r
    Gv = (0, 0, 0, 0)
    prev_q = None
    reads = []; verdicts = []
    settled = None; dist = []
    bg = set(_r.Random(seed).sample(CELLS, back)) if back else set()
    k = 0
    for t in range(frames):
        path = shape_at(t)
        gx, gy = U.xy(Gv)
        if not path:
            E = set(); dist.append(float('nan'))
            Eb = E | bg
            c, n = flood(Eb & NETS) if (Eb & NETS) else (0, 0)
            reads.append(1 if n >= 1 else 0)
            if len(reads) >= 3:
                verdicts.append(1 if sum(reads[-3:]) >= 2 else 0)
            if settled is None and len(verdicts) >= 2 and verdicts[-1] == verdicts[-2]:
                settled = (t + 1, verdicts[-1])
            w, ctr, nd = sink(Eb, (X0, Y0), RNET)
            q = near(*ctr) if ctr else None
            if q is not None and q != prev_q:
                Gv = U.zadd(Gv, tuple(a - b for a, b in zip(q, CTR))); prev_q = q
            else:
                Gv = U.zadd(Gv, UNITS[k % 10]); k += 1; prev_q = None
            continue
        E = paint2([(x - gx, y - gy) for x, y in path])
        Eb = E | bg
        # 対象の中心（採点用・装置は使わない）
        cx = sum(x for x, y in path) / len(path) - gx
        cy = sum(y for x, y in path) / len(path) - gy
        dist.append(math.hypot(cx - X0, cy - Y0))
        # ③ 読む
        c, n = flood(Eb & NETS) if (Eb & NETS) else (0, 0)
        reads.append(1 if n >= 1 else 0)
        if len(reads) >= 3:
            verdicts.append(1 if sum(reads[-3:]) >= 2 else 0)
        if settled is None and len(verdicts) >= 2 and verdicts[-1] == verdicts[-2]:
            settled = (t + 1, verdicts[-1])
        # ① 深いところへ寄る
        w, ctr, nd = sink(Eb, (X0, Y0), RNET)
        q = near(*ctr) if ctr else None
        if q is not None and q != prev_q:
            Gv = U.zadd(Gv, tuple(a - b for a, b in zip(q, CTR)))
            prev_q = q
        else:
            # ② 止まったので振る
            Gv = U.zadd(Gv, UNITS[k % 10]); k += 1
            prev_q = None
    return settled, sum(dist[-10:]) / 10, reads

CASES = [
    ("止まった丸 R=6.5",        lambda t: circle(6.5, 260), 0),
    ("止まった三角の輪",         lambda t: FIG["閉じた輪(三角)"], 0),
    ("止まったC字",             lambda t: FIG["開いた輪(C字)"], 0),
    ("止まった一本の線",         lambda t: FIG["一本の線"], 0),
    ("動く丸（一歩1.0）",        lambda t: circle(6.5, 260, cx=t * 1.0), 0),
    ("動く丸（一歩3.0）",        lambda t: circle(6.5, 260, cx=t * 1.5, cy=t * 2.598), 0),
    ("対象なし",                lambda t: [], 0),
    ("止まった丸＋背景100枚",     lambda t: circle(6.5, 260), 100),
    ("止まったC字＋背景100枚",    lambda t: FIG["開いた輪(C字)"], 100),
    ("対象なし＋背景100枚",       lambda t: [], 100),
]

print("  %-22s %10s %8s %10s" % ("場面", "固まった枚", "判定", "対象との距離"))
for nm, f, bk in CASES:
    st, d, reads = run(f, 30, bk)
    print("  %-22s %10s %8s %10.2f"
          % (nm, ("%d枚" % st[0]) if st else "固まらず",
             ("閉じている" if st[1] else "開いている") if st else "—",
             float('nan') if d != d else d))

print("\n  ※ 対象との距離は採点用（装置は対象の真の中心を使っていない）")
