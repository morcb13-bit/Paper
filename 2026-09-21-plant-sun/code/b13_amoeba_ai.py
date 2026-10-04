# b13_amoeba_ai.py — アメーバ型ロボットで網膜の棒（縦か横か）を読ませる
#
# 床：b13_retina_ai.py と同じ（菱形の床の頂点 90 個、左側 32 個が網膜）。辺は両向き。
# 出口：右側に二つ（上の出口 A＝「横」、下の出口 B＝「縦」）。
# 匂いの層：出口ごとに一枚。D_t(i) = Σ 隣の D_{t−1} ＋（出口なら k^t）を整数で回す（第19章と同じ。k は隣の数の最大＋1）。
# ロボット：棒で点いた網膜のセルに一体ずつ置く（点くのは網膜の 1 割前後）。
#   各刻、自分の層の匂いがいちばん濃い隣へ一歩進む。
#   通票：一つのセルに入れるのは一刻に一体だけ。先に札を取った一体が入る（番号の若い順）。
#   待避：いちばん濃い隣が埋まっていたら、二番目に濃い隣が今より濃く、空いていればそこへ（駅の二本目）。
#         どちらも塞がっていたら、その場で一刻待つ（待機所）。
#   乗り換え：待たされたセルに「乗り換え」の札があれば、そのロボットは層を替える（A↔B）。
#   出口に着いたロボットは床から降りる。
# 読み：出口 A に着いた数 ＞ 出口 B に着いた数 なら「横」。
# 学ぶもの（すべて 0/1 の札）：
#   ・網膜のセルごとに、そこから出るロボットの最初の層（A か B か）
#   ・床のセルごとに、待たされたら層を乗り換えるか
#   外れが減る札だけ裏返す。外れが同じなら、着いた数の差（当たりの側に足し、外れの側に引く）で比べる。
#
# 基準（走らせる前に決めたもの）
#   M1  通票あり：見せていない棒の当たりが、前回の網（約 69%）を上回ること
#   M0  通票なし（ロボットが重なってすり抜ける）＝セルごとの票を数えるだけの多数決。M1 がこれを上回れば、
#       渋滞のさばき方（待つ・待避・乗り換え）が読みに効いている
#   素朴な読み：点いたセルの数だけで当てる最良の閾値
#   対照 札を混ぜると半分前後
#   NG  M1 が 69% 以下、または M0 と差がない
import math, random, sys, functools
import numpy as np
import b13_retina_ai as RA
from b13_iconet import build, fpos, xpos
from b13_icolearn import sign as zsign, sub as zsub

def make_floor():
    net, cells, verts, d = RA.make_retina_net()
    idx = {v: i for i, v in enumerate(verts)}
    N = len(verts)
    adj = [set() for _ in range(N)]
    for a, b in net['fwd']: adj[a].add(b); adj[b].add(a)
    adj = [sorted(s) for s in adj]
    pix = net['inputs'][:-1]                         # 網膜（源は使わない）
    # 出口：右側 1/4 のうち、y がいちばん大きいセル（A）といちばん小さいセル（B）
    xs = [fpos(v)[0] for v in verts]; x0, x1 = min(xs), max(xs)
    right = [i for i, v in enumerate(verts) if fpos(v)[0] > x1 - 0.25 * (x1 - x0)]
    A = max(right, key=lambda i: fpos(verts[i])[1]); B = min(right, key=lambda i: fpos(verts[i])[1])
    return dict(N=N, adj=adj, pix=pix, A=A, B=B, verts=verts), cells

def scent(fl, ex, T=40):
    """床は二部グラフなので、一刻ごとに濃さが偶奇で入れ替わる。続く二刻を足して偶奇をならす"""
    N = fl['N']; D = [0] * N; prev = D
    K = max(len(a) for a in fl['adj']) + 1        # 第19章：k が隣の数の最大以上なら最短の段数どおりに着く
    for t in range(T + 1):
        nd = [sum(D[j] for j in fl['adj'][i]) for i in range(N)]
        nd[ex] += K ** t
        prev, D = D, nd
    return [a + b for a, b in zip(D, prev)]

def simulate(fl, DS, img, start, swap, tokens=True, T=40):
    """img：網膜の 0/1。start[k]：網膜 k 番の最初の層。swap[i]：セル i で待たされたら層を替えるか"""
    rob = [[fl['pix'][k], start[k]] for k in range(len(img)) if img[k]]   # [位置, 層 0=A 1=B]
    arrived = [0, 0]; exits = (fl['A'], fl['B'])
    for t in range(T):
        if not rob: break
        occ = {} if not tokens else {r[0]: n for n, r in enumerate(rob)}
        claimed = set()
        nxt = []
        for n, (pos, lay) in enumerate(rob):
            D = DS[lay]
            nb = sorted(fl['adj'][pos], key=lambda j: -D[j])
            moved = False
            for c in nb[:2]:
                if D[c] <= D[pos]: break
                if tokens and (c in claimed or (c in occ and c not in exits)): continue
                claimed.add(c); pos2 = c; moved = True; break
            if not moved:
                pos2 = pos
                if tokens and swap[pos]: lay = 1 - lay
            if pos2 == exits[lay] or pos2 in exits:
                arrived[0 if pos2 == exits[0] else 1] += 1
            else:
                nxt.append([pos2, lay])
        rob = nxt
    return arrived

def predict(fl, DS, X, start, swap, tokens):
    out = []
    for img in X:
        a, b = simulate(fl, DS, img, start, swap, tokens)
        out.append((a, b))
    return out

def score(out, y):
    hit = 0; marg = 0
    for (a, b), yy in zip(out, y):
        pred = 1 if a > b else 0
        hit += int(pred == yy)
        d = (a - b) if yy == 1 else (b - a)
        marg += d
    return hit, marg

def learn(fl, DS, X, y, tokens, rng, passes=8):
    P = len(fl['pix']); N = fl['N']
    start = [rng.randrange(2) for _ in range(P)]
    swap = [0] * N
    cur = score(predict(fl, DS, X, start, swap, tokens), y)
    for _ in range(passes):
        moved = False
        moves = [("s", k) for k in range(P)] + ([("w", i) for i in range(N)] if tokens else [])
        rng.shuffle(moves)
        for kind, k in moves:
            s2, w2 = start, swap
            if kind == "s": s2 = start.copy(); s2[k] ^= 1
            else: w2 = swap.copy(); w2[k] ^= 1
            sc = score(predict(fl, DS, X, s2, w2, tokens), y)
            if sc > cur: start, swap, cur, moved = s2, w2, sc, True
        if not moved: break
    return start, swap

if __name__ == "__main__":
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    NTE = 400
    fl, cells = make_floor()
    DS = [scent(fl, fl['A']), scent(fl, fl['B'])]
    print(f"床：セル {fl['N']}、網膜 {len(fl['pix'])}、出口 A={fl['A']} B={fl['B']}", flush=True)
    tot = [0, 0, 0, 0]
    for s in range(seeds):
        rng = random.Random(7000 + s)
        Xtr, ytr = RA.make_data(cells, NTR, rng); Xte, yte = RA.make_data(cells, NTE, rng)
        Xtr = Xtr[:, :-1]; Xte = Xte[:, :-1]
        res = []
        for tok in (True, False):
            st, sw = learn(fl, DS, Xtr, ytr, tok, rng)
            a_tr = score(predict(fl, DS, Xtr, st, sw, tok), ytr)[0]
            a_te = score(predict(fl, DS, Xte, st, sw, tok), yte)[0]
            res.append((a_tr, a_te, sum(sw)))
        ysh = list(ytr); rng.shuffle(ysh); ysh = np.array(ysh)
        stc, swc = learn(fl, DS, Xtr, ysh, True, rng)
        yte_sh = np.array([rng.randrange(2) for _ in yte])
        ctl = score(predict(fl, DS, Xte, stc, swc, True), yte_sh)[0]
        nv = RA.naive(np.hstack([Xtr, np.ones((len(Xtr), 1), int)]), ytr, np.hstack([Xte, np.ones((len(Xte), 1), int)]), yte)
        lit = round(float(Xte.sum(1).mean()), 1)
        tot[0] += res[0][1]; tot[1] += res[1][1]; tot[2] += nv; tot[3] += ctl
        print(f"  {s}：通票あり 学習用 {res[0][0]}/{NTR} 試験 {res[0][1]}/{NTE}（乗り換えの札 {res[0][2]}）  ｜ "
              f"通票なし 学習用 {res[1][0]}/{NTR} 試験 {res[1][1]}/{NTE}  ｜ 素朴 {nv}/{NTE}  ｜ 対照 {ctl}/{NTE}  ｜ 点いたセルの平均 {lit}", flush=True)
    print(f"合計（試験 {NTE*seeds}）：通票あり {tot[0]}  通票なし {tot[1]}  素朴 {tot[2]}  対照 {tot[3]}")
