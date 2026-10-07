# 検定②：扇10枚の担体に、隙間（舟・細ひし形・五芒星）を番地として足す。各隙間は囲む五角形とだけ繋ぐ。
# 正十角形は足さない。規則・読み（最初に着いた出口）・種・枚数・判定は b13_amoeba10_first.py と同じ。
# 基準（走らせる前に決めたもの）
#   合格  通票あり − 通票なし ＞ 2 ポイント、かつ通票ありが点の数だけの読みを上回る
#   負の対照  ラベルを混ぜる → 33% 前後
#   NG  差が 2 ポイント以内、または点の数だけの読み以下
import sys, pickle, random
from collections import Counter
import b13_amoeba10_first as A

KINDS = tuple(sys.argv[3].split(',')) if len(sys.argv) > 3 else ("舟", "細ひし形", "五芒星")
_base = A.make_floor
def make_floor():
    fl = _base()
    gl = pickle.load(open('/home/claude/icolearn/gaps_floor.pkl', 'rb'))
    adj = [set(a) for a in fl['adj']]; xy = list(fl['xy']); fan = list(fl['fan']); kind = ['五角形'] * fl['N']
    for name, bd in gl:
        if name not in KINDS: continue
        g = len(adj); adj.append(set(bd))
        for b in bd: adj[b].add(g)
        xy.append((sum(xy[b][0] for b in bd) / len(bd), sum(xy[b][1] for b in bd) / len(bd)))
        fan.append(sorted(set(sum((fan[b] for b in bd), [])))); kind.append(name)
    fl.update(N=len(adj), adj=[sorted(s) for s in adj], xy=xy, fan=fan, kind=kind)
    return fl
A.make_floor = make_floor

if __name__ == "__main__":
    seeds = int(sys.argv[1]); NTR = int(sys.argv[2]); NTE = 400
    fl = A.make_floor(); DS = [A.scent(fl, e) for e in fl['ex']]
    print(f"床：番地 {fl['N']}（{dict(Counter(fl['kind']))}）、網膜 {len(fl['pix'])}、出口 {fl['ex']}", flush=True)
    tot = Counter()
    for s in range(seeds):
        rng = random.Random(8000 + s)
        Xtr, ytr = A.make_data(fl, NTR, rng); Xte, yte = A.make_data(fl, NTE, rng)
        res = []
        for tok in (True, False):
            st, sw = A.learn(fl, DS, Xtr, ytr, tok, rng)
            res.append(A.score(fl, DS, Xte, yte, st, sw, tok)[0])
        ysh = list(ytr); rng.shuffle(ysh)
        stc, swc = A.learn(fl, DS, Xtr, ysh, True, rng)
        ctl = A.score(fl, DS, Xte, yte, stc, swc, True)[0]
        nv = A.naive(Xtr, ytr, Xte, yte)
        tot['tok'] += res[0]; tot['free'] += res[1]; tot['nv'] += nv; tot['ctl'] += ctl
    n = NTE * seeds
    print(f"学習用 {NTR}（試験 {n}）：通票あり {tot['tok']}  通票なし {tot['free']}  点の数 {tot['nv']}  対照 {tot['ctl']}")
