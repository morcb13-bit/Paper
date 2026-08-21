#  検定SR3  中心を与えずに、表を丸ごと引けるか
#
#      検定SR2 で、探索の中心に真の位置を与えると 60枚全部を追い、ずれ 0.44 だった。
#      落ちていたのは表ではなく中心の与え方である。ならば中心を与えなければよい。
#      毎枚 150ビットを読み、表 1222 件の全部に対して「符牒の 1 が全部立っているか」を
#      問う。予測も速さも使わない。力学を持たない。
#
#      対象：検定TB と同じ表・同じ場面。追跡の状態を一切持たない読み
#      測る量：乗った項目の件数／それらの位置の広がり／1が最も多い項目と真の位置の距離
#      OK なら：乗った項目が一箇所に集まり（広がり ≤ 1.618）、距離が CH3 のずれ
#              （1.55〜3.07）より小さい ── 追跡は力学なしの表引きだけで足りる
#      NG なら：離れた場所の項目が同時に乗る ── 中心が要る。要る理由は曖昧さであって
#              予測の精度ではない
#      必ず落ちる設定：対象のいない画面では 0 件（検定SR1 で確認済み・ここでも再掲）

exec(open('wind_core.py').read())

import math, random

XYC = {q: U.xy(q) for q in CELLS}
FGRID = {}
for q in CELLS:
    x, y = XYC[q]
    FGRID.setdefault((int(x // 3), int(y // 3)), []).append(q)

def fast_land(p):
    fx, fy = U.xy(tuple(float(t) for t in p))
    gx, gy = int(fx // 3), int(fy // 3)
    best = None; bd = None
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for q in FGRID.get((gx + dx, gy + dy), ()):
                x, y = XYC[q]
                d = (x - fx) ** 2 + (y - fy) ** 2
                if bd is None or d < bd:
                    bd = d; best = q
    return best

B1 = U.xy((1, 0, 0, 0)); B2 = U.xy((0, 1, 0, 0))
DET = B1[0] * B2[1] - B1[1] * B2[0]
DEN = 10 ** 12
def to_field(x, y):
    a = (x * B2[1] - y * B2[0]) / DET
    b = (B1[0] * y - B1[1] * x) / DET
    return (F(int(a * DEN), DEN), F(int(b * DEN), DEN), F(0), F(0))

GAP = 1.618034
O = U.xy(CTR)

rows, place, offs = U.build_stack()
faces = U.gaps(U.fits(sum(place, [])))
STARS = []
for area, cyc in faces:
    if abs(area - 2.9389) < 0.01:
        P = [U.xy(p) for p in cyc]
        STARS.append((sum(a for a, _ in P) / len(P), sum(b for _, b in P) / len(P)))
AROUND = {}
for g in STARS:
    ds = sorted(CELLS, key=lambda q: (XYC[q][0] - g[0]) ** 2 + (XYC[q][1] - g[1]) ** 2)[:5]
    AROUND[g] = tuple(sorted(ds, key=lambda q: math.atan2(XYC[q][1] - g[1], XYC[q][0] - g[0])))

TV = []
for k in range(3):
    th = math.radians(90 + 120 * k)
    o, _ = find_offset(6.0 * math.cos(th), 6.0 * math.sin(th))
    TV.append(o)
PARTS = []
for k in range(3):
    a, b = TV[k], TV[(k + 1) % 3]
    for i in range(8):
        t = F(i, 8)
        PARTS.append(tuple(F(a[j]) + t * (b[j] - a[j]) for j in range(4)))
C0 = tuple(sum(p[j] for p in PARTS) / len(PARTS) for j in range(4))
PARTS = [U.zsub(p, C0) for p in PARTS]

def put_at(x, y):
    base = to_field(x, y)
    out = set()
    for p in PARTS:
        q = fast_land(tuple(a + b for a, b in zip(base, p)))
        if q is not None:
            out.add(q)
    return out

def on_bits(E):
    on = []
    for g in STARS:
        for q in AROUND[g]:
            if q in E:
                on.append(q)
    return on

def traj(T, k, speed=0.45, turn=60.0, seed=1):
    rnd = random.Random(seed)
    x, y = O[0] - 3.0, O[1] - 2.0
    ang = rnd.uniform(0, 2 * math.pi)
    out = []
    for t in range(T):
        if k and t % k == 0 and t:
            ang += math.radians(rnd.uniform(-turn, turn))
        x += speed * math.cos(ang); y += speed * math.sin(ang)
        if math.hypot(x - O[0], y - O[1]) > 5.0:
            ang = math.atan2(O[1] - y, O[0] - x) + rnd.uniform(-0.4, 0.4)
        out.append((x, y))
    return out

def scenes(pts, seed, nback=100, nflick=8):
    back = set(random.Random(seed).sample(CELLS, nback))
    return [back | set(random.Random(50 + 7 * t).sample(CELLS, nflick)) | put_at(x, y)
            for t, (x, y) in enumerate(pts)]

def bake(minbits, step=0.30, n=44):
    tbl = []
    for i in range(n):
        for j in range(n):
            x, y = O[0] - 6.5 + i * step, O[1] - 6.5 + j * step
            on = on_bits(put_at(x, y))
            if len(on) < minbits:
                continue
            tbl.append(((x, y), tuple(on), len(on)))
    return tbl

def pull_all(E, tbl, allow):
    out = []
    for pos, on, n1 in tbl:
        miss = 0
        for q in on:
            if q not in E:
                miss += 1
                if miss > allow:
                    break
        if miss <= allow:
            out.append((pos, n1))
    return out

for MB in (6, 8):
    tbl = bake(MB)
    print("■ 表：1 が %d 本以上 → %d 件" % (MB, len(tbl)))
    print("  欠け  乗った件数  広がり(最大)  1が最多の項目のずれ  重心のずれ  空の枚")
    for allow in (0, 1):
        cnt = []; spread = []; dbest = []; dcent = []; empty = 0; n = 0
        for sd in range(5):
            pts = traj(60, 4, seed=sd + 1)
            Es = scenes(pts, sd + 1)
            for t, (x, y) in enumerate(pts):
                hits = pull_all(Es[t], tbl, allow)
                n += 1
                if not hits:
                    empty += 1
                    continue
                cnt.append(len(hits))
                xs = [p[0] for p, _ in hits]; ys = [p[1] for p, _ in hits]
                sp = 0.0
                for i in range(len(hits)):
                    for j in range(i + 1, len(hits)):
                        sp = max(sp, math.hypot(xs[i] - xs[j], ys[i] - ys[j]))
                spread.append(sp)
                bp = max(hits, key=lambda h: h[1])[0]
                dbest.append(math.hypot(bp[0] - x, bp[1] - y))
                cx, cy = sum(xs) / len(xs), sum(ys) / len(ys)
                dcent.append(math.hypot(cx - x, cy - y))
        m = len(cnt)
        print("   %d    %6.1f     %6.2f        %6.2f            %6.2f      %.3f"
              % (allow, sum(cnt) / m, sum(spread) / m, sum(dbest) / m,
                 sum(dcent) / m, empty / n))
    # 必ず落ちる設定
    hit = 0; n = 0
    for sd in range(5):
        back = set(random.Random(sd + 1).sample(CELLS, 100))
        for t in range(60):
            E = back | set(random.Random(50 + 7 * t).sample(CELLS, 8))
            if pull_all(E, tbl, 1):
                hit += 1
            n += 1
    print("  必ず落ちる設定（対象なし・欠け1本まで許す）：乗った枚 %d / %d\n" % (hit, n))
