# b13_spinnet.py — 表裏つきの指し（20 で一周）で、扇の輪の鏡の罠が外れるか
#
# 指し j ∈ Z20：扇は j mod 10、表裏は j div 10（扇を 10 枚回ると表裏が入れ替わり、20 で元に戻る）。
# 一歩：届いた指しの扇を数えて、いちばん近い扇に着地（並んだら黙る）。
#       表裏は、届いた指しのうち裏の数の偶奇（裏が奇数なら裏）。
#       自分の向き r ∈ Z20 を足して隣へ渡す。
# 入力：+1 → j = 0（扇 0・表）、−1 → j = 10（扇 0・裏）、0 → 黙る。
#       −1 は「半周（扇 5 枚）」ではなく「表裏の入れ替え」で表す。
# 読み：出力の番地が表で着地したら「はい」、裏か黙りなら「いいえ」。
# 学習：外れが減るときだけ、番地の向きを ±1（扇一枚）か +10（表裏の入れ替え）だけ変える。
#
# 基準（走らせる前に決めたもの）
#   S0  入力を全部裏返したとき、10 扇の網は出力も必ず裏返る（前回の罠）。表裏つきの網ではそうならないこと
#   S1  「同符号」の試験の当たりが、10 扇の網の上限（半分）を明らかに上回ること（6 割以上）
#   対照 札を混ぜると、試験の当たりは半分前後に落ちること
#   NG  S1 で試験が半分前後なら、鏡の罠のほかに（繋がりの側に）原因がある
import itertools, random, sys
import numpy as np
import b13_fannet as F

def run(net, edges, r, X, T):
    N = net['N']; Ex = X.shape[0]
    src = np.array([a for a, b in edges], dtype=np.int64)
    dst = np.array([b for a, b in edges], dtype=np.int64)
    inp = np.array(net['inputs'])
    xin = np.where(X > 0, 0, np.where(X < 0, 10, -1))
    def clamp(P):
        P[:, inp] = np.where(xin >= 0, (xin + r[inp][None, :]) % 20, -1)
    P = -np.ones((Ex, N), dtype=np.int64); clamp(P)
    for t in range(T):
        C = np.zeros((Ex, N, 10), dtype=np.int64)
        Z = np.zeros((Ex, N), dtype=np.int64)            # 届いた裏の数
        on = P[:, src] >= 0
        e_idx, k_idx = np.nonzero(on)
        pj = P[e_idx, src[k_idx]]
        np.add.at(C, (e_idx, dst[k_idx], pj % 10), 1)
        np.add.at(Z, (e_idx, dst[k_idx]), pj // 10)
        K = F.land(C)
        J = np.where(K >= 0, K + 10 * (Z % 2), -1)       # 着地した指し（向きを足す前）
        P = np.where(J >= 0, (J + r[None, :]) % 20, -1); clamp(P)
    return J[:, net['out']]

def score(out, y):
    pred = ((out >= 0) & (out < 10)).astype(np.int64)    # 表で着地＝はい
    return int((pred == y).sum())

def learn(net, edges, X, y, T, rng, passes=12):
    N = net['N']
    r = np.array([rng.randrange(20) for _ in range(N)], dtype=np.int64)
    cur = score(run(net, edges, r, X, T), y)
    for _ in range(passes):
        moved = False
        order = list(range(N)); rng.shuffle(order)
        for v in order:
            for d in (1, -1, 10):
                r2 = r.copy(); r2[v] = (r2[v] + d) % 20
                s2 = score(run(net, edges, r2, X, T), y)
                if s2 > cur: r, cur, moved = r2, s2, True; break
        if not moved: break
    return r

TASKS = {
    "同符号": (lambda X: (X[:, 0] * X[:, 1] == 1).astype(np.int64), lambda A: A[(A[:, 0] != 0) & (A[:, 1] != 0)]),
    "全体の積が正": (lambda X: (np.prod(np.where(X == 0, 1, X), axis=1) == 1).astype(np.int64), lambda A: A),
    "一つ目が正": (lambda X: (X[:, 0] == 1).astype(np.int64), lambda A: A[(A[:, 0] != 0) & (A[:, 1] != 0)]),
}

if __name__ == "__main__":
    net = F.make_floor(); m = len(net['inputs'])
    ALL0 = np.array([s for s in itertools.product((-1, 0, 1), repeat=m) if any(s)])
    NETS = {"A 一方向": (net['fwd'], net['L']), "C 隣＋φ²": (F.directed(net['nb'] + net['phi2']), net['L'] + 4)}
    # S0：入力を全部裏返したときの出力
    print("S0 入力を全部裏返したとき（無作為の向き 3 通り × 網 2 通り、入力 728 通り）")
    for nname, (edges, T) in NETS.items():
        for s in range(3):
            rng = random.Random(s)
            r10 = np.array([rng.randrange(10) for _ in range(net['N'])])
            r20 = np.array([rng.randrange(20) for _ in range(net['N'])])
            o10a = F.run(net, edges, r10, ALL0, T); o10b = F.run(net, edges, r10, -ALL0, T)
            o20a = run(net, edges, r20, ALL0, T); o20b = run(net, edges, r20, -ALL0, T)
            f10 = int(((o10a >= 0) & (o10b == (o10a + 5) % 10)).sum() + ((o10a < 0) & (o10b < 0)).sum())
            p20a = (o20a >= 0) & (o20a < 10); p20b = (o20b >= 0) & (o20b < 10)
            same20 = int((p20a == p20b).sum())
            print(f"  {nname} 種{s}：10扇で出力が半周して裏返った（または両方黙る） {f10}/728  ｜ "
                  f"表裏つきで読みが変わらなかった {same20}/728")
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    for tname, (f, sel) in TASKS.items():
        A = sel(ALL0); y_all = f(A)
        print(f"\n題「{tname}」 入力 {len(A)} 通り（はい {int(y_all.sum())}）、学習用 {NTR}、{seeds} 回の合計", flush=True)
        for nname, (edges, T) in NETS.items():
            trh = teh = ctl = tot = 0
            for s in range(seeds):
                rng = random.Random(2000 + s); idx = list(range(len(A))); rng.shuffle(idx)
                tr, te = np.array(idx[:NTR]), np.array(idx[NTR:])
                r = learn(net, edges, A[tr], y_all[tr], T, rng)
                trh += score(run(net, edges, r, A[tr], T), y_all[tr])
                teh += score(run(net, edges, r, A[te], T), y_all[te])
                perm = list(range(len(A))); rng.shuffle(perm); ysh = y_all[perm]
                rc = learn(net, edges, A[tr], ysh[tr], T, rng)
                ctl += score(run(net, edges, rc, A[te], T), ysh[te]); tot += len(te)
            maj = max(int(y_all.sum()), len(A) - int(y_all.sum()))
            print(f"  {nname}：学習用 {trh}/{NTR*seeds}  試験 {teh}/{tot}  ｜ 対照 {ctl}/{tot}"
                  f"（多数決の割合 {maj}/{len(A)}）", flush=True)
