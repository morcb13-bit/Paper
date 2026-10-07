# 測定（判定なし）：種ごとに、通票ありの模型で「違う層のロボットに道を塞がれた待機」と「同じ層に塞がれた待機」を層ごとに数える。
# 使い方：python3 seed_probe.py warp|warp3b
import sys
which = sys.argv[1]; sys.argv = [sys.argv[0], '0', '0']
import random
from collections import Counter
if which == 'warp':
    import b13_amoeba10_warp as M; A = M.A
    fl = A.make_floor(); DS = [A.scent(fl, e) for e in fl['ex']]
else:
    import b13_amoeba10_warp3b as M; A = M.A; import b13_amoeba10_warp3 as W3
    fl = A.make_floor(); DS = [W3.scent_l(fl, l) for l in range(3)]
L = 3; own = fl.get('own'); twin = fl.get('twin')
def g_of(p):
    k = fl['kind'][p]; return k[:2] if k.startswith('層') else '層0'
def probe(img, start, swap, C, T=120):
    rob = [[fl['pix'][k], start[k]] for k in range(len(img)) if img[k]]
    exset = set(fl['ex']); first = None; arrived = [0] * L
    for t in range(T):
        if not rob: break
        occ = {r[0]: r[1] for r in rob}; claimed = {}; nxt = []
        for pos, lay in rob:
            D = DS[lay]; nb = sorted(fl['adj'][pos], key=lambda j: -D[j]); pos2 = None; who = None
            for c in nb[:2]:
                if D[c] <= D[pos]: break
                if c in claimed: who = who if who is not None else claimed[c]; continue
                if c in occ and c not in exset: who = who if who is not None else occ[c]; continue
                claimed[c] = lay; pos2 = c; break
            if pos2 is None:
                pos2 = pos
                if who is not None: C[(g_of(pos), '違う層' if who != lay else '同じ層')] += 1
                if swap[pos]:
                    lay = (lay + 1) % L
                    if own is not None and own[pos] is not None: pos2 = twin[pos][lay]
            if pos2 in exset: arrived[fl['ex'].index(pos2)] += 1
            else: nxt.append([pos2, lay])
        if first is None:
            got = {j for j in range(L) if arrived[j]}
            if got: first = got.pop() if len(got) == 1 else -1
        rob = nxt
    return first
print(which)
for s in range(10):
    rng = random.Random(8000 + s)
    Xtr, ytr = A.make_data(fl, 240, rng); Xte, yte = A.make_data(fl, 400, rng)
    st, sw = A.learn(fl, DS, Xtr, ytr, True, rng)
    C = Counter(); hit = Counter(); tie = 0
    for img, yy in zip(Xte, yte):
        f = probe(img, st, sw, C); hit[yy] += int(f == yy); tie += int(f == -1)
    sc = Counter(st)
    row = "  ".join(f"{g} 違{C[(g,'違う層')]}/同{C[(g,'同じ層')]}" for g in ('層0', '層1', '層2'))
    print(f"種{s}：当たり {sum(hit.values())}（A{hit[0]} B{hit[1]} C{hit[2]}）同時 {tie}  最初の層 A{sc[0]} B{sc[1]} C{sc[2]}  乗り換え札 {sum(sw)}  ｜ {row}", flush=True)
