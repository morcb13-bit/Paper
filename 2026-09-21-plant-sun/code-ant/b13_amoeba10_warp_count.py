# 測定（判定なし）：③の床で、層ごとに入ったロボット・一体一刻の手・待機を数える。通票あり／なしの両方。
# 模型は③と同じ種・同じ学習（学習用 240）。見せていない試験 400 枚 × 10 回。
import sys; sys.argv = [sys.argv[0], '0', '0']
import random
from collections import Counter, defaultdict
import b13_amoeba10_warp as W
A = W.A
fl = A.make_floor(); DS = [A.scent(fl, e) for e in fl['ex']]; L = 3
lay = ['層0' if k in ('五角形', '舟', '細ひし形', '五芒星') else k for k in fl['kind']]
C = {t: defaultdict(Counter) for t in (True, False)}
def run(img, start, swap, tok, T=120):
    rob = [[fl['pix'][k], start[k], set()] for k in range(len(img)) if img[k]]
    exset = set(fl['ex'])
    for t in range(T):
        if not rob: break
        occ = set(r[0] for r in rob) if tok else set(); claimed = set(); nxt = []
        for pos, ly, seen in rob:
            g = lay[pos]; c = C[tok][g]; c['手'] += 1
            if g not in seen: seen.add(g); c['入った体'] += 1
            D = DS[ly]; nb = sorted(fl['adj'][pos], key=lambda j: -D[j]); pos2 = None
            for n, cc in enumerate(nb[:2]):
                if D[cc] <= D[pos]: break
                if tok and (cc in claimed or (cc in occ and cc not in exset)): continue
                claimed.add(cc); pos2 = cc
                if n == 1: c['待避'] += 1
                break
            if pos2 is None:
                pos2 = pos; c['待機'] += 1
                if tok and swap[pos]: ly = (ly + 1) % L; c['乗り換え'] += 1
            if pos2 in exset: C[tok]['全体']['着いた'] += 1
            else: nxt.append([pos2, ly, seen])
        rob = nxt
    C[tok]['全体']['着かない'] += len(rob)
for s in range(10):
    rng = random.Random(8000 + s)
    Xtr, ytr = A.make_data(fl, 240, rng); Xte, yte = A.make_data(fl, 400, rng)
    models = [A.learn(fl, DS, Xtr, ytr, tok, rng) for tok in (True, False)]
    for tok, (st, sw) in zip((True, False), models):
        for img in Xte:
            C[tok]['全体']['体'] += sum(img); run(img, st, sw, tok)
nlay = Counter(lay)
for tok in (True, False):
    g = C[tok]['全体']
    print(f"通票{'あり' if tok else 'なし'}：ロボット {g['体']}、着いた {g['着いた']}、120刻で着かない {g['着かない']}")
    for name in ('層0', '層1', '層2'):
        c = C[tok][name]
        print(f"  {name}（{nlay[name]}番地）：入った体 {c['入った体']}  手 {c['手']}  待機 {c['待機']}  待避 {c['待避']}  乗り換え {c['乗り換え']}")
