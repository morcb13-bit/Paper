# b13_bondnet.py — ペンタゴン・リッチフローAIプロセッサ（最初の形）
#
# セル（床の番地）は表裏つきの指し（20 で一周）を持ち、届いた指しを扇に着地させる（b13_spinnet と同じ）。
# 隣どうしのセルは「結合」したときだけ信号を通す（面を共有して仕切りを取り除いた＝手術のあと）。
# 学習の一手は三つ：
#   ・セルの向きを扇一枚回す（±1）／表裏を入れ替える（+10）
#   ・結合を一つ作る／一つ外す
#   外れが減るなら採る。外れが同じなら「結合を外す」手だけ採る（使っていない繋がりは縮んで消える）。
#
# 基準（走らせる前に決めたもの）
#   B1  見せていない入力の当たりが、全部開いた網（前回 NG、半分前後）を明らかに上回ること（6 割以上）
#   対照 札を混ぜると半分前後に落ちること
#   B2  「同符号」を学んだあと、題に関係ない入力（3〜6 番目）から出力へ、結合を辿って届く道が残っているか数える
#       （残っていれば、当たりが上がっても理由は切り離しではない）
#   NG  B1 で半分前後なら、結合を選んでも足りない
import itertools, random, sys
import numpy as np
import b13_fannet as F
import b13_spinnet as S

def learn(net, pairs, directed_mode, X, y, T, rng, passes=10):
    N = net['N']; P = len(pairs)
    r = np.array([rng.randrange(20) for _ in range(N)], dtype=np.int64)
    bond = np.ones(P, dtype=bool)                     # 始めは全部結合（前回の網）
    def edges_of(b):
        sel = [pairs[i] for i in range(P) if b[i]]
        return sel if directed_mode else F.directed(sel)
    def ev(r, b): return S.score(S.run(net, edges_of(b), r, X, T), y)
    cur = ev(r, bond)
    for _ in range(passes):
        moved = False
        moves = [("r", v, d) for v in range(N) for d in (1, -1, 10)] + [("b", i, 0) for i in range(P)]
        rng.shuffle(moves)
        for kind, a, d in moves:
            if kind == "r":
                r2 = r.copy(); r2[a] = (r2[a] + d) % 20; b2 = bond
            else:
                r2 = r; b2 = bond.copy(); b2[a] = not b2[a]
            s2 = ev(r2, b2)
            if s2 > cur or (s2 == cur and kind == "b" and not b2[a]):
                r, bond, moved = r2, b2, (moved or s2 > cur or kind == "b")
                cur = s2
        if not moved: break
    return r, bond, edges_of(bond)

def reaches(net, edges, src):
    adj = {}
    for a, b in edges: adj.setdefault(a, []).append(b)
    seen = {src}; fr = [src]
    while fr:
        nf = []
        for u in fr:
            for w in adj.get(u, []):
                if w not in seen: seen.add(w); nf.append(w)
        fr = nf
    return net['out'] in seen

if __name__ == "__main__":
    net = F.make_floor(); m = len(net['inputs'])
    ALL0 = np.array([s for s in itertools.product((-1, 0, 1), repeat=m) if any(s)])
    NETS = {"A 一方向": (net['fwd'], True, net['L']),
            "C 隣＋φ²": (net['nb'] + net['phi2'], False, net['L'] + 4)}
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    only = sys.argv[3] if len(sys.argv) > 3 else None
    for tname, (f, sel) in S.TASKS.items():
        if only and tname != only: continue
        A = sel(ALL0); y_all = f(A)
        print(f"題「{tname}」 入力 {len(A)} 通り（はい {int(y_all.sum())}）、学習用 {NTR}、{seeds} 回の合計", flush=True)
        for nname, (pairs, dmode, T) in NETS.items():
            trh = teh = ctl = tot = 0; bonds = []; leak = []
            for s in range(seeds):
                rng = random.Random(3000 + s); idx = list(range(len(A))); rng.shuffle(idx)
                tr, te = np.array(idx[:NTR]), np.array(idx[NTR:])
                r, b, E = learn(net, pairs, dmode, A[tr], y_all[tr], T, rng)
                trh += S.score(S.run(net, E, r, A[tr], T), y_all[tr])
                teh += S.score(S.run(net, E, r, A[te], T), y_all[te])
                bonds.append(int(b.sum()))
                leak.append([reaches(net, E, net['inputs'][i]) for i in range(m)])
                perm = list(range(len(A))); rng.shuffle(perm); ysh = y_all[perm]
                rc, bc, Ec = learn(net, pairs, dmode, A[tr], ysh[tr], T, rng)
                ctl += S.score(S.run(net, Ec, rc, A[te], T), ysh[te]); tot += len(te)
            print(f"  {nname}：学習用 {trh}/{NTR*seeds}  試験 {teh}/{tot}  ｜ 対照 {ctl}/{tot}  ｜ "
                  f"残った結合 {bonds}/{len(pairs)}  入力1〜6から出力へ届く道 "
                  + " ".join("".join("○" if x else "・" for x in L) for L in leak), flush=True)
