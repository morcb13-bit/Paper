# b13_retina_ai.py — 床の左側を網膜にして棒の画像を載せ、縦か横かを学ばせる
#
# 床：ペンタグリッドの菱形の床（半径 R）。左側のセルが網膜、右端のセルが出力。
# 網膜のセルは、自分の中心が棒の上にあれば +1（扇0・表）、なければ黙る。左端の一つは源（いつも +1）。
# 信号は出力からの段数が一つ減る向きにだけ流れる（一方向）。セルは伸び縮みする（b13_ricci.py と同じ）。
# 学習：位相 ±1／表裏の入れ替え／周波数 ±1／閉じたまま止める。外れが減る手を採る（止める手は同点でも採る）。
# 読み：出力が表で着地したら「横」。
#
# 画像を作るところ（棒の位置・向きとセルの中心の比較）だけ浮動小数を使う。これは題材のデータで、装置の外。
#
# 基準（走らせる前に決めたもの）
#   I1  見せていない位置の棒で、試験の当たりが 7 割以上
#   対照 札を混ぜると半分前後
#   素朴な読み：点いたセルの数だけで縦横を当てる最良の閾値（学習用で選ぶ）の試験の当たり。
#       I1 がこれを上回らなければ、網は数を数えただけで形を読んでいない
#   NG  試験が半分前後、または素朴な読み以下
import math, random, sys, functools
import numpy as np
from b13_iconet import build, fpos, xpos
from b13_icolearn import sign as zsign, sub as zsub
import b13_spinnet as S, b13_ricci as R
R.CLOSE['on'] = True

def make_retina_net(Rad=5.0, frac=0.4):
    verts, edges = build(6)
    verts = [v for v in verts if math.hypot(*fpos(v)) < Rad]
    vs = set(verts)
    adj = {v: set() for v in verts}
    for u, v in edges:
        if u in vs and v in vs: adj[u].add(v); adj[v].add(u)
    order = sorted(verts, key=functools.cmp_to_key(lambda a, b: zsign(zsub(xpos(a), xpos(b)))))
    out = order[-1]
    d = {out: 0}; fr = [out]
    while fr:
        nf = []
        for u in fr:
            for w in adj[u]:
                if w not in d: d[w] = d[u] + 1; nf.append(w)
        fr = nf
    verts = [v for v in order if v in d]
    idx = {v: i for i, v in enumerate(verts)}
    xs = [fpos(v)[0] for v in verts]; x0, x1 = min(xs), max(xs)
    retina = [v for v in verts if fpos(v)[0] < x0 + frac * (x1 - x0)]
    src = retina[0]                              # 左端＝源
    pix = retina[1:]
    fwd = [(idx[u], idx[w]) for u in verts for w in adj[u] if d[w] == d[u] - 1]
    net = dict(N=len(verts), inputs=[idx[v] for v in pix] + [idx[src]], out=idx[out],
               L=max(d[v] for v in retina), fwd=fwd, nb=[(a, b) for a, b in fwd])
    return net, [fpos(v) for v in pix], verts, d

def draw_bar(cells, horiz, cx, cy, length=2.6, width=0.9):
    a, b = (length / 2, width / 2) if horiz else (width / 2, length / 2)
    return np.array([1 if abs(x - cx) <= a and abs(y - cy) <= b else 0 for x, y in cells], dtype=np.int64)

def make_data(cells, n, rng):
    xs = [c[0] for c in cells]; ys = [c[1] for c in cells]
    X, y = [], []
    while len(X) < n:
        h = rng.randrange(2)
        cx = rng.uniform(min(xs) + 0.6, max(xs) - 0.6); cy = rng.uniform(min(ys) + 0.6, max(ys) - 0.6)
        img = draw_bar(cells, h, cx, cy)
        if img.sum() < 2: continue                # ほとんど網膜の外の棒は捨てる
        X.append(np.append(img, 1)); y.append(h) # 最後の列＝源
    return np.array(X), np.array(y, dtype=np.int64)

def naive(Xtr, ytr, Xte, yte):
    ctr = Xtr[:, :-1].sum(1); cte = Xte[:, :-1].sum(1); best = (0, 0, 1)
    for th in range(0, int(ctr.max()) + 2):
        for sgn in (1, -1):
            h = int((((ctr >= th) if sgn > 0 else (ctr < th)).astype(int) == ytr).sum())
            if h > best[0]: best = (h, th, sgn)
    _, th, sgn = best
    return int((((cte >= th) if sgn > 0 else (cte < th)).astype(int) == yte).sum())

if __name__ == "__main__":
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    NTE = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    OFF = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    PASSES = int(sys.argv[5]) if len(sys.argv) > 5 else 15
    net, cells, verts, d = make_retina_net()
    T = net['L'] + 3
    print(f"床：セル {net['N']} 個、網膜 {len(cells)} 個＋源 1、出力までの段数 最大 {net['L']}、一方向の辺 {len(net['fwd'])}", flush=True)
    for s in range(OFF, OFF + seeds):
        rng = random.Random(6000 + s)
        Xtr, ytr = make_data(cells, NTR, rng); Xte, yte = make_data(cells, NTE, rng)
        p, w = R.learn(net, net['fwd'], Xtr, ytr, T, rng, passes=PASSES)
        a_tr = S.score(R.run(net, net['fwd'], p, w, Xtr, T), ytr)
        a_te = S.score(R.run(net, net['fwd'], p, w, Xte, T), yte)
        nv = naive(Xtr, ytr, Xte, yte)
        ysh = ytr.copy(); rng.shuffle(ysh)
        pc, wc = R.learn(net, net['fwd'], Xtr, ysh, T, rng, passes=PASSES)
        yte_sh = np.array([rng.randrange(2) for _ in yte])
        ctl = S.score(R.run(net, net['fwd'], pc, wc, Xte, T), yte_sh)
        closed = int(((w == 0) & (p % 2 == 0)).sum())
        print(f"  {s}：学習用 {a_tr}/{NTR}  試験 {a_te}/{NTE}  ｜ 素朴な読み {nv}/{NTE}  ｜ 対照 {ctl}/{NTE}  ｜ 閉じたセル {closed}/{net['N']}", flush=True)
