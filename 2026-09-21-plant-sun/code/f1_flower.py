# 検定F1：花びら 3・4・5・6 枚の花を、どの向き・どの位置でも読む（ハチ型ロボットの目と脳）
# 床：②の床（扇10枚の担体＋舟・細ひし形・五芒星）。網膜を中心の半径 12 に広げる。出口 4 つ（縁の 45°・135°・225°・315° に最も近い五角形）。
# 花：中心から伸びる n 枚の楕円の花びら（長さ 6）。花びらの幅は 5/n に比例させ、花全体の面積を n によらずそろえる。
#     向きは 0〜360° の任意、花の中心は網膜の中心から半径 3 以内のでたらめな位置。
# 学習：K1a と同じ積み重ね（240 枚ずつ 4 回、前回の札を出発点に、これまでに見たもの全部で続きを学ぶ）。規則・読み（最初に着いた出口）は②と同じ。
# 基準（走らせる前に決めたもの）
#   合格  見せていない花 4000 枚で、当たりが偶然の 25% と点の数だけの読みを、どちらも 10 ポイントを越えて上回る
#   あわせて測る：花びらの数ごとの当たり（5 枚だけ良くても合格とは言わない）
#   負の対照  ラベルを混ぜて 25% 前後
#   NG  上のどちらかを越えない
import sys; sys.argv = [sys.argv[0], '0', '0']
import math, random, json
from collections import Counter
import b13_amoeba10 as A0
import b13_amoeba10_gap as G
A = G.A
G._base = lambda: A0.make_floor(ret=12.0)
NS = (3, 4, 5, 6); L = 4

def make_floor():
    fl = G.make_floor()
    rim = [n for n in range(fl['N']) if fl['kind'][n] == '五角形' and math.hypot(*fl['xy'][n]) > 17]
    ex = []
    for ang in (45, 135, 225, 315):
        ux, uy = math.cos(math.radians(ang)), math.sin(math.radians(ang))
        ex.append(max(rim, key=lambda n: fl['xy'][n][0] * ux + fl['xy'][n][1] * uy))
    fl['ex'] = ex
    return fl
fl = make_floor(); DS = [A.scent(fl, e) for e in fl['ex']]
P = len(fl['pix']); N = fl['N']

def flower(ci, rng, R=6.0, B0=1.3):
    n = NS[ci]; th = rng.uniform(0, 2 * math.pi)
    rr = 3 * math.sqrt(rng.random()); ph = rng.uniform(0, 2 * math.pi); cx, cy = rr * math.cos(ph), rr * math.sin(ph)
    a = R / 2; b = B0 * 5 / n; img = []
    for k in fl['pix']:
        x, y = fl['xy'][k]; dx, dy = x - cx, y - cy; on = 0
        for i in range(n):
            t = th + 2 * math.pi * i / n; c, s = math.cos(t), math.sin(t)
            u = dx * c + dy * s - a; v = -dx * s + dy * c
            if (u / a) ** 2 + (v / b) ** 2 <= 1: on = 1; break
        img.append(on)
    return img
def make_data(n, rng):
    X, y = [], []
    for _ in range(n):
        ci = rng.randrange(L); X.append(flower(ci, rng)); y.append(ci)
    return X, y

def learn_from(X, y, tokens, rng, start, swap, passes=8):
    W = set(); cur = A.score(fl, DS, X, y, start, swap, tokens, W)
    for _ in range(passes):
        moved = False
        moves = [("s", k, v) for k in range(P) for v in range(L)] + ([("w", i, 0) for i in sorted(W)] if tokens else [])
        rng.shuffle(moves)
        for kind, k, v in moves:
            s2, w2 = start, swap
            if kind == "s":
                if start[k] == v: continue
                s2 = start.copy(); s2[k] = v
            else: w2 = swap.copy(); w2[k] ^= 1
            W2 = set(); sc = A.score(fl, DS, X, y, s2, w2, tokens, W2)
            if sc > cur: start, swap, cur, moved = s2, w2, sc, True; W |= W2
        if not moved: break
    return start, swap

if __name__ == "__main__" and len(sys.argv) > 0:
    pass

# ---- 速い学習：札を一枚裏返したとき、結果が変わりうる画像だけを走らせ直す（結果は全部走らせ直すのと同じ）----
def one(img, yy, st, sw, tokens):
    W = set(); arr, f = A.simulate(fl, DS, img, st, sw, tokens, waits=W)
    other = max(arr[j] for j in range(L) if j != yy)
    return int(f == yy), arr[yy] - other, W

def learn_fast(X, y, tokens, rng, start, swap, passes=8):
    lit = [[i for i, img in enumerate(X) if img[k]] for k in range(P)]
    res = [one(img, yy, start, swap, tokens) for img, yy in zip(X, y)]
    tot = lambda R: (sum(r[0] for r in R), sum(r[1] for r in R))
    cur = tot(res)
    for _ in range(passes):
        moved = False
        waitcells = set().union(*(r[2] for r in res)) if tokens else set()
        moves = [("s", k, v) for k in range(P) for v in range(L)] + [("w", i, 0) for i in sorted(waitcells)]
        rng.shuffle(moves)
        for kind, k, v in moves:
            s2, w2 = start, swap
            if kind == "s":
                if start[k] == v: continue
                s2 = start.copy(); s2[k] = v; idx = lit[k]
            else:
                w2 = swap.copy(); w2[k] ^= 1; idx = [i for i, r in enumerate(res) if k in r[2]]
            new = {i: one(X[i], y[i], s2, w2, tokens) for i in idx}
            sc = (cur[0] + sum(new[i][0] - res[i][0] for i in idx), cur[1] + sum(new[i][1] - res[i][1] for i in idx))
            if sc > cur:
                start, swap, cur, moved = s2, w2, sc, True
                for i in idx: res[i] = new[i]
        if not moved: break
    return start, swap

def test(st, sw, Xte, yte, tokens=True):
    per = Counter(); hit = 0
    for img, yy in zip(Xte, yte):
        _, f = A.simulate(fl, DS, img, st, sw, tokens)
        if f == yy: hit += 1; per[yy] += 1
    return hit, per

def job(s):
    import time
    Xte, yte = make_data(4000, random.Random(12345))
    rng = random.Random(8000 + s)
    batches = [make_data(240, rng) for _ in range(4)]
    init = [rng.randrange(L) for _ in range(P)]
    out = []
    ch = {'tok': (init[:], [0] * N), 'ctl': (init[:], [0] * N)}
    sx, sy, sysh = [], [], []
    for r, (Xb, yb) in enumerate(batches):
        t0 = time.time()
        sx += Xb; sy += yb; ysh = list(yb); rng.shuffle(ysh); sysh += ysh
        row = {}
        for k in ch:
            st, sw = learn_fast(sx, sy if k == 'tok' else sysh, True, rng, ch[k][0][:], ch[k][1][:])
            h, per = test(st, sw, Xte, yte)
            row[k] = h; row[k + '_花びら別'] = [per[c] for c in range(L)]; ch[k] = (st, sw)
        row['点の数'] = A.naive(sx, sy, Xte, yte)
        row['各数の試験枚数'] = [yte.count(c) for c in range(L)]
        out.append(row); print(s, r + 1, row, f"{int(time.time()-t0)}秒", flush=True)
    return s, out, ch['tok']

if __name__ == "__main__":
    from multiprocessing import Pool
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] != '0' else 5
    with Pool(2) as p: res = p.map(job, range(seeds), chunksize=1)
    json.dump(res, open('f1_result.json', 'w'))
    for r in range(4):
        tok = sum(o[r]['tok'] for _, o, _ in res); ctl = sum(o[r]['ctl'] for _, o, _ in res); nv = sum(o[r]['点の数'] for _, o, _ in res)
        per = [sum(o[r]['tok_花びら別'][c] for _, o, _ in res) for c in range(L)]
        cnt = [sum(o[r]['各数の試験枚数'][c] for _, o, _ in res) for c in range(L)]
        print(f"■ {r+1} 回目（見せた {240*(r+1)} 枚、試験 4000×{len(res)}）：通票あり {tok}  点の数 {nv}  対照 {ctl}  ｜ 花びら 3/4/5/6 の当たり {per}（各 {cnt}）")
