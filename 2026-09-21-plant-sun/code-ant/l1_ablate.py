# 検定AB：L1 の通票ありの模型から、乗り換えの札を場所ごとに外して、固定の試験 4000 枚の当たりを測る（学習はやり直さない）
# 基準（走らせる前に決めたもの）：外して当たりが 2 ポイントを越えて下がれば、その場所の札は読みに使われている。2 ポイント以内なら効いていない。
import sys; sys.argv = [sys.argv[0], '0', '0']
import json, random
from collections import Counter
import b13_amoeba10_gap as G
A = G.A
fl = A.make_floor(); DS = [A.scent(fl, e) for e in fl['ex']]
Xte, yte = A.make_data(fl, 4000, random.Random(12345))
res = json.load(open('l1_result.json')); pixset = set(fl['pix'])
def kinds(st, sw):
    C = Counter()
    for img, yy in zip(Xte, yte):
        arr, f = A.simulate(fl, DS, img, st, sw, True)
        if f == yy: C['当たり'] += 1
        elif f is None: C['誰も着かない'] += 1
        elif f == -1: C['同時着']  += 1
        elif arr[yy] == 0: C['正解の出口に誰も着かない'] += 1
        else: C['正解にも着くが別が先'] += 1
    return C
for n in (240, 960):
    tot = {k: Counter() for k in ('そのまま', '(a)床を外す', '(b)網膜を外す', '(c)両方を外す')}
    for nn, s, r in res:
        if nn != n: continue
        st, sw = r['model']
        variants = {'そのまま': sw,
                    '(a)床を外す': [0 if i not in pixset else v for i, v in enumerate(sw)],
                    '(b)網膜を外す': [0 if i in pixset else v for i, v in enumerate(sw)],
                    '(c)両方を外す': [0] * len(sw)}
        for k, w in variants.items(): tot[k] += kinds(st, w)
    print(f"■ 学習用 {n}（試験 4000×10＝40000）")
    for k, C in tot.items():
        print(f"  {k:10s} 当たり {C['当たり']}  ｜ 正解の出口に誰も着かない {C['正解の出口に誰も着かない']}  正解にも着くが別が先 {C['正解にも着くが別が先']}  同時着 {C['同時着']}  誰も着かない {C['誰も着かない']}", flush=True)
