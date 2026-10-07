# 検定K1：積み重ねる学習。240 枚ずつ 4 回に分けて見せ、前回の札を出発点に続きを学ぶ（白紙に戻さない）。
#   K1a：これまでに見たもの全部で続きを学ぶ   K1b：新しく見せた 240 枚だけで続きを学ぶ
# 基準（走らせる前に決めたもの）
#   K1-1  回を重ねるごとに当たりが上がる（1 < 2 < 3 < 4 回目）
#   K1-2  4 回目の当たりが、白紙から 960 枚で学んだ L1（19635/40000 = 49.1%）に並ぶか上回る
#   K1-3  通票あり − 通票なし が各回で 2 ポイントを越える
#   負の対照  ラベルを混ぜて積み重ねる（K1a の形）→ 33% 前後
#   NG  回を重ねても上がらない、または K1b で前に覚えたことが崩れて下がる
import sys; sys.argv = [sys.argv[0], '0', '0']
import random, json
from multiprocessing import Pool
import b13_amoeba10_gap as G
A = G.A
fl = A.make_floor(); DS = [A.scent(fl, e) for e in fl['ex']]
Xte, yte = A.make_data(fl, 4000, random.Random(12345)); pixset = set(fl['pix'])
P = len(fl['pix']); N = fl['N']

def learn_from(X, y, tokens, rng, start, swap, passes=8):
    W = set(); cur = A.score(fl, DS, X, y, start, swap, tokens, W)
    for _ in range(passes):
        moved = False
        moves = [("s", k, v) for k in range(P) for v in range(3)] + ([("w", i, 0) for i in sorted(W)] if tokens else [])
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

def test(st, sw, tok): return A.score(fl, DS, Xte, yte, st, sw, tok)[0]

def job(s):
    rng = random.Random(8000 + s)
    batches = [A.make_data(fl, 240, rng) for _ in range(4)]
    init = [rng.randrange(3) for _ in range(P)]
    chains = {k: (init[:], [0] * N) for k in ('a_tok', 'a_free', 'b_tok', 'b_free', 'ctl')}
    out = []
    seenX, seenY, seenYsh = [], [], []
    for r, (Xb, yb) in enumerate(batches):
        seenX += Xb; seenY += yb
        ysh = list(yb); rng.shuffle(ysh); seenYsh += ysh
        row = {}
        for k in chains:
            tok = k.endswith('tok') or k == 'ctl'
            X, y = (seenX, seenY) if k.startswith('a') else (Xb, yb) if k.startswith('b') else (seenX, seenYsh)
            old = chains[k]
            st, sw = learn_from(X, y, tok, rng, old[0][:], old[1][:])
            row[k] = test(st, sw, tok)
            row[k + '_変わった札'] = sum(a != b for a, b in zip(st, old[0])) + sum(a != b for a, b in zip(sw, old[1]))
            if k in ('a_tok', 'b_tok'):
                row[k + '_網膜を外す'] = test(st, [0 if i in pixset else v for i, v in enumerate(sw)], True)
                row[k + '_床を外す'] = test(st, [0 if i not in pixset else v for i, v in enumerate(sw)], True)
                row[k + '_乗り換え札'] = sum(sw)
            chains[k] = (st, sw)
        out.append(row); print(s, r + 1, row, flush=True)
    return s, out

if __name__ == "__main__":
    with Pool(2) as p: res = p.map(job, range(10), chunksize=1)
    json.dump(res, open('k1_result.json', 'w'))
    keys = list(res[0][1][0].keys())
    for r in range(4):
        tot = {k: sum(o[r][k] for _, o in res) for k in keys}
        print(f"■ {r+1} 回目（見せた {240*(r+1)} 枚、試験 40000）")
        print("   " + "  ".join(f"{k} {v}" for k, v in tot.items()))
