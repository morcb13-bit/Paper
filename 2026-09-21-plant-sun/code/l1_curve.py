# 検定L1：②の模型（隙間を足した床・最初に着いた出口）を離散の学習器として、学習枚数 240/480/960 で固定の試験 4000 枚に当てる。
# 基準（走らせる前に決めたもの）
#   L1a  通票ありの当たりが 240 < 480 < 960 と上がる
#   L1b  各枚数で 通票あり − 通票なし ＞ 2 ポイント
#   負の対照  ラベルを混ぜる → 33% 前後
#   NG  枚数を増やしても通票ありが上がらない、または点の数だけの読み以下
#   参考（合否に使わない）：最近傍（食い違いの数）、パーセプトロン（外れたら ±1）
import sys; sys.argv = [sys.argv[0], '0', '0']
import random, json
from multiprocessing import Pool
import b13_amoeba10_gap as G
A = G.A
fl = A.make_floor(); DS = [A.scent(fl, e) for e in fl['ex']]
Xte, yte = A.make_data(fl, 4000, random.Random(12345))

def nn(Xtr, ytr, X):
    out = []
    for img in X:
        best = None
        for t, yy in zip(Xtr, ytr):
            d = sum(a != b for a, b in zip(img, t))
            if best is None or d < best[0]: best = (d, yy)
        out.append(best[1])
    return out

def perceptron(Xtr, ytr, X, rng, epochs=30):
    P = len(Xtr[0]); W = [[0] * (P + 1) for _ in range(3)]
    sc = lambda w, img: w[P] + sum(w[k] for k in range(P) if img[k])
    order = list(range(len(Xtr)))
    for _ in range(epochs):
        rng.shuffle(order)
        for i in order:
            s = [sc(W[c], Xtr[i]) for c in range(3)]; p = s.index(max(s))
            if p != ytr[i]:
                for k in range(P):
                    if Xtr[i][k]: W[ytr[i]][k] += 1; W[p][k] -= 1
                W[ytr[i]][P] += 1; W[p][P] -= 1
    return [max(range(3), key=lambda c: sc(W[c], img)) for img in X]

def job(args):
    ntr, s = args
    rng = random.Random(8000 + s)
    Xtr, ytr = A.make_data(fl, ntr, rng)
    r = {}
    for tok in (True, False):
        st, sw = A.learn(fl, DS, Xtr, ytr, tok, rng)
        r['tok' if tok else 'free'] = A.score(fl, DS, Xte, yte, st, sw, tok)[0]
        if tok: r['model'] = (st, sw)
    ysh = list(ytr); rng.shuffle(ysh)
    stc, swc = A.learn(fl, DS, Xtr, ysh, True, rng)
    r['ctl'] = A.score(fl, DS, Xte, yte, stc, swc, True)[0]
    r['nv'] = A.naive(Xtr, ytr, Xte, yte)
    r['nn'] = sum(int(p == t) for p, t in zip(nn(Xtr, ytr, Xte), yte))
    r['pc'] = sum(int(p == t) for p, t in zip(perceptron(Xtr, ytr, Xte, random.Random(s)), yte))
    print(ntr, s, {k: v for k, v in r.items() if k != 'model'}, flush=True)
    return ntr, s, r

if __name__ == "__main__":
    jobs = [(n, s) for n in (240, 480, 960) for s in range(10)]
    with Pool(2) as p: res = p.map(job, jobs, chunksize=1)
    json.dump([[n, s, r] for n, s, r in res], open('l1_result.json', 'w'))
    for n in (240, 480, 960):
        rs = [r for nn_, s, r in res if nn_ == n]
        tot = {k: sum(r[k] for r in rs) for k in ('tok', 'free', 'nv', 'ctl', 'nn', 'pc')}
        print(f"学習用 {n}（試験 4000×10＝40000）：" + "  ".join(f"{k} {v}" for k, v in tot.items()))
