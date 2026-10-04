# 記事・動く図の素材：学習済みの札、待ち・乗り換えの数え上げ、図用の床と匂いの順位を書き出す
import json, random, numpy as np
import b13_amoeba_ai as M, b13_retina_ai as RA
from b13_iconet import fpos

def sim_trace(fl, DS, img, start, swap, tokens=True, T=40):
    """simulate と同じ規則。待った回数・待避・乗り換えをセルごとに数える"""
    rob = [[fl['pix'][k], start[k]] for k in range(len(img)) if img[k]]
    arrived = [0, 0]; exits = (fl['A'], fl['B'])
    waits = {}; sides = 0; swaps = {}
    for t in range(T):
        if not rob: break
        occ = {r[0]: n for n, r in enumerate(rob)} if tokens else {}
        claimed = set(); nxt = []
        for n, (pos, lay) in enumerate(rob):
            D = DS[lay]; nb = sorted(fl['adj'][pos], key=lambda j: -D[j]); moved = False
            for r_, c in enumerate(nb[:2]):
                if D[c] <= D[pos]: break
                if tokens and (c in claimed or (c in occ and c not in exits)): continue
                claimed.add(c); pos2 = c; moved = True; sides += (r_ == 1); break
            if not moved:
                pos2 = pos; waits[pos] = waits.get(pos, 0) + 1
                if tokens and swap[pos]: lay = 1 - lay; swaps[pos] = swaps.get(pos, 0) + 1
            if pos2 in exits: arrived[0 if pos2 == exits[0] else 1] += 1
            else: nxt.append([pos2, lay])
        rob = nxt
    return arrived, waits, sides, swaps

fl, cells = M.make_floor()
DS = [M.scent(fl, fl['A']), M.scent(fl, fl['B'])]
rng = random.Random(7000)
Xtr, ytr = RA.make_data(cells, 240, rng); Xte, yte = RA.make_data(cells, 400, rng)
Xtr = Xtr[:, :-1]; Xte = Xte[:, :-1]
st, sw = M.learn(fl, DS, Xtr, ytr, True, rng)
te = M.score(M.predict(fl, DS, Xte, st, sw, True), yte)[0]
print("学習用240・試験", te, "/400")
# 縦横ごとの待ち・待避・乗り換え
agg = {0: [0, 0, 0, 0], 1: [0, 0, 0, 0]}
swapcell = {}
for img, yy in zip(Xte, yte):
    arr, waits, sides, swaps = sim_trace(fl, DS, img, st, sw)
    a = agg[int(yy)]; a[0] += 1; a[1] += sum(waits.values()); a[2] += sides; a[3] += sum(swaps.values())
    for c, n in swaps.items(): swapcell.setdefault(c, [0, 0])[int(yy)] += n
for yy, name in ((1, "横"), (0, "縦")):
    n, w, s, x = agg[yy]
    print(f"{name}の棒 {n} 枚：待った回数 {w}（1枚あたり {w/n:.2f}）  待避 {s}  乗り換え {x}")
print("乗り換えの札が立ったセル", [i for i in range(fl['N']) if sw[i]])
print("札ごとの乗り換え回数 [縦, 横]", {c: v for c, v in sorted(swapcell.items())})
# 図用の書き出し
pos = [fpos(v) for v in fl['verts']]
up = [[[j for j in sorted(fl['adj'][i], key=lambda j: -D[j]) if D[j] > D[i]] for i in range(fl['N'])] for D in DS]
pred = M.predict(fl, DS, Xte[:40], st, sw, True)
pred0 = M.predict(fl, DS, Xte[:40], st, sw, False)
out = dict(pos=[[round(x, 5), round(y, 5)] for x, y in pos],
           edges=sorted({(min(i, j), max(i, j)) for i in range(fl['N']) for j in fl['adj'][i]}),
           pix=fl['pix'], A=fl['A'], B=fl['B'], up=up, start=st, swap=sw,
           check=dict(imgs=Xte[:40].tolist(), y=yte[:40].tolist(), tok=pred, notok=pred0),
           acc=dict(test=te, n=400))
json.dump(out, open("amoeba_fig.json", "w"))
print("書き出し amoeba_fig.json")
